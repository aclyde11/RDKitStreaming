#include <iostream>
#include <fstream>
#include <bitset>
#include <cinttypes>
#include <string>
#include <chrono>
#include <ctime>

#include <parallel_hashmap/phmap.h>

#include <nop/serializer.h>
#include <nop/utility/die.h>
#include <nop/utility/stream_reader.h>
#include <nop/utility/stream_writer.h>


#include <StreamingMoleculeReader/parallelDatabase.h>


using parallel_smile_set = phmap::parallel_flat_hash_set<std::string, std::hash<std::string>, std::equal_to<std::string>, std::allocator<std::string>, 8, std::mutex>;
using namespace SMR;

// for ceralizing
namespace {
    auto Die() { return nop::Die(std::cerr); }
}

void task(std::string const& item, MutexCounter *total_counter, MutexCounter *valid_counter, MutexCounter *unqiue_counter, parallel_smile_set *meset) {
    boost::optional<std::string> value = getCannonicalSmileFromSmile(item);
    total_counter->increment();
    if (value.has_value()) {
        valid_counter->increment();
        if (std::get<1>(meset->insert(value.value())))
            unqiue_counter->increment();
    }
}

void monitarThread(std::vector<MutexCounter> * valid_counters,
                   std::vector<MutexCounter> * unique_counters,
                   std::vector<MutexCounter> * total_counters,  int monitar_freq, std::atomic<bool> *stop, moodycamel::ConcurrentQueue<std::string> *q) {

    std::cout <<"time,queuesize,total,valid,unique"<<std::endl;
    while(!(stop->load(std::memory_order_acquire))) {
        std::this_thread::sleep_for(std::chrono::seconds(monitar_freq));
        int total = 0;
        int valid = 0;
        int unique = 0;

        for(size_t i = 0; i < valid_counters->size(); i++) {
            total += (*total_counters)[i].view();
            valid += (*valid_counters)[i].view();
            unique += (*unique_counters)[i].view();
        }

        std::cout << std::time(nullptr) << "," << q->size_approx() << "," << total << "," << valid << "," << unique << std::endl;
    }
}

// number threads, load, file_saved, file_read
int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "check your args. Exiting." << std::endl;
        return 1;
    }

    bool LOAD = (bool)(atoi(argv[2]));
    size_t n_threads = (size_t)(atoi(argv[1]));

    // parallel test.
    moodycamel::ConcurrentQueue<std::string> q;
    std::vector<std::thread> threads;

    std::atomic<int> doneConsumers(0);

    std::vector<MutexCounter> total_counters(n_threads);
    std::vector<MutexCounter> valid_counters(n_threads);
    std::vector<MutexCounter> unique_counters(n_threads);

    std::unordered_map<std::string, bool> initial_set;
    parallel_smile_set dbase;

    if (LOAD) {
        if (argc < 5) {
            std::cout << "check your args. Exiting." << std::endl;
            return 1;
        }

        std::cout << "Starting here." << std::endl;
        initial_set = getInitialSetFromFile(n_threads, 1000000000);
        std::cout << "Got set." << std::endl;
        using Writer = nop::StreamWriter<std::ofstream>;
        nop::Serializer<Writer> serializer{argv[3]};
        serializer.Write(initial_set) || Die();
        serializer.writer().stream().close();
        std::cout << "wrote out to file " << argv[3] << std::endl;
        return 0;
    } else {
        using Reader = nop::StreamReader<std::ifstream>;
        std::cout << "reading in " << argv[3] << " as initial databse of SMILES." << std::endl;
        nop::Deserializer<Reader> deserializer{argv[3]};
        deserializer.Read(&initial_set) || Die();
        std::cout << initial_set.size() << std::endl;
    }

    //convert STL to phmap
    for (std::pair<std::string, bool>  element : initial_set)
    {
        dbase.insert(element.first);
    }

    std::cout << "new dbase has initial size " << dbase.size() << " Moving on to your problem." <<std::endl;

    //Producer Thread
    std::atomic<bool> doneProducer(true);
    threads[0] = std::thread([&](std::atomic<bool> *stop) {
        for (std::string line; std::getline(std::cin, line);) {
            q.enqueue(line);
        }
        *stop = false;
    }, &doneProducer);


    //Monitar Thread
    std::atomic<bool> stopMonitar(false);
    std::thread monitar(
            [&](std::vector<MutexCounter> * valid_counters,
                std::vector<MutexCounter> * unique_counters,
                std::vector<MutexCounter> * total_counters,  std::atomic<bool> *stop) {

                monitarThread(valid_counters, unique_counters, total_counters, 5, stop, &q);
            }, &valid_counters, &unique_counters, &total_counters, &stopMonitar);

    // Consumers
    for (size_t i = 1; i != n_threads; ++i) {
        threads.emplace_back(
                [&](MutexCounter *total_counter, MutexCounter *valid_counter, MutexCounter *unqiue_counter,
                        std::atomic<bool> *stop, parallel_smile_set *meset) {
                    std::string item;
                    boost::optional<std::string> value;
                    bool itemsLeft;
                    do {
                        // It's important to fence (if the producers have finished) *before* dequeueing
                        itemsLeft = stop->load(std::memory_order_acquire);
                        while (q.try_dequeue(item)) {
                            itemsLeft = true;
                            task(item, total_counter, valid_counter, unqiue_counter, meset);
                        }
                    } while (itemsLeft || doneConsumers.fetch_add(1, std::memory_order_acq_rel) + 1 == (n_threads - 1));
                }, &total_counters[i],
                   &valid_counters[i],
                   &unique_counters[i],
                   &doneProducer, &dbase);
    }

    // Wait for all threads
    for (size_t i = 0; i != n_threads; ++i) {
        threads[i].join();
    }

    std::cout << "ok everyone else finsihed. Waiting for monitar" << std::endl;
    stopMonitar = true;
    monitar.join();

    // Collect any leftovers (could be some if e.g. consumers finish before producers)
    std::string item;
    boost::optional<std::string> value;
    while (q.try_dequeue(item)) {
        task(item, &total_counters[0], &valid_counters[0], &unique_counters[0], &dbase);
    }

    int unique = 0;
    int total = 0;
    int valid = 0;
    for (size_t i = 0; i < n_threads; i++) {
        unique += unique_counters[i].view();
        total += total_counters[i].view();
        valid += valid_counters[i].view();

    }

    std::cout << "total,unique,valid" << std::endl;
    std::cout << total << ", " << unique << ", " << valid << std::endl;
    std::cout << dbase.size() << std::endl;
}
