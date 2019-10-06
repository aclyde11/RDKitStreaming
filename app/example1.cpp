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


using parallel_smile_set = phmap::parallel_flat_hash_set<std::string, std::hash<std::string>, std::equal_to<std::string>, std::allocator<std::string>, 8>;
using namespace SMR;

// for ceralizing
namespace {
    auto Die() { return nop::Die(std::cerr); }
}

void task(std::string const& item, MutexCounter *total_counter, MutexCounter *valid_counter, moodycamel::ConcurrentQueue<std::string> *qout) {
    boost::optional<std::string> value = getCannonicalSmileFromSmile(item);
    total_counter->increment();
    if (value.has_value()) {
        valid_counter->increment();
        qout->enqueue(value.value());
    }
}

void monitarThread(std::vector<MutexCounter> * valid_counters,
                   std::vector<MutexCounter> * unique_counters,
                   std::vector<MutexCounter> * total_counters,
                   std::vector<MutexCounter> * enamine_counter,
                   int monitar_freq, std::atomic<bool> *stop, moodycamel::ConcurrentQueue<std::string> *q, moodycamel::ConcurrentQueue<std::string> *q_out) {

    std::cout <<"time,queuesize,writersize,total,valid,unique,uniqueenamine"<<std::endl;
    while(!(stop->load(std::memory_order_acquire))) {
        std::this_thread::sleep_for(std::chrono::seconds(monitar_freq));
        int total = 0;
        int valid = 0;
        int unique = 0;
        int enamine = 0;

        for(size_t i = 0; i < valid_counters->size(); i++) {
            total += (*total_counters)[i].view();
            valid += (*valid_counters)[i].view();
            unique += (*unique_counters)[i].view();
            enamine += (*enamine_counter)[i].view();
        }

        std::cout << std::time(nullptr) << "," << q->size_approx() << "," << q_out->size_approx() << "," << total << "," << valid << "," << unique <<"," << enamine << std::endl;
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
    moodycamel::ConcurrentQueue<std::string> q_out;

    std::vector<std::thread> threads;

    std::atomic<int> doneConsumers(0);

    std::vector<MutexCounter> total_counters(n_threads);
    std::vector<MutexCounter> valid_counters(n_threads);
    std::vector<MutexCounter> unique_counters(n_threads);
    std::vector<MutexCounter> enamine_unique_counters(n_threads);

    std::unordered_map<std::string, bool> initial_set, big_initial_set;
    parallel_smile_set dbase, dbase_enamine;

    if (LOAD) {
        std::cout << "Starting here." << std::endl;
        if (argc == 5) {
            std::cout << "Reading from file." << std::endl;
            initial_set = getInitialSetFromFile(argv[4], n_threads, 1000000000);
        } else {
            initial_set = getInitialSetFromFile(n_threads, 1000000000);
        }

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

     {
        using Reader = nop::StreamReader<std::ifstream>;
        std::cout << "reading in " << argv[4] << " as ENAMINE databse of SMILES." << std::endl;
        nop::Deserializer<Reader> deserializer{argv[4]};
        big_initial_set.reserve(10000000000);
        deserializer.Read(&big_initial_set) || Die();
        std::cout << initial_set.size() << std::endl;
    }

    dbase_enamine.reserve(big_initial_set.size());
    int couttester=0;
    //convert STL to phmap
    for (std::pair<std::string, bool>  element : big_initial_set)
    {
        dbase_enamine.insert(element.first);
        couttester++;
        if (couttester % 1000000 == 0) {
            std::cout << "si" << std::endl;
        }
    }
    big_initial_set.clear();
    std::cout << "loadded huge map into struct." << std::endl;

    //convert STL to phmap
    for (std::pair<std::string, bool>  element : initial_set)
    {
        dbase.insert(element.first);
    }

    std::cout << "Cleaning..." << std::endl;
    initial_set.clear();

    std::cout << "new dbase has initial size " << dbase.size() << " Moving on to your problem." <<std::endl;
    std::cout << "new dbase_enaine has initial size " << dbase.size() << " Moving on to your problem." <<std::endl;

    //Producer Thread
    std::atomic<bool> doneProducer(true);
    threads.emplace_back([&](std::atomic<bool> *stop) {
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
                std::vector<MutexCounter> * total_counters,
                std::vector<MutexCounter> * enamine_counters,
                std::atomic<bool> *stop) {

                monitarThread(valid_counters, unique_counters, total_counters, enamine_counters, 5, stop, &q, &q_out);
            }, &valid_counters, &unique_counters, &total_counters, &enamine_unique_counters, &stopMonitar);


    //Writer Thread
    std::atomic<bool> stopWriter(false);
    std::thread writer(
            [&]( std::atomic<bool> *stop, parallel_smile_set *meset, parallel_smile_set *enamine) {
                while (!(stop->load(std::memory_order_acquire))) {
                    std::string item;
                    while (q_out.try_dequeue(item)) {
                        if (std::get<1>(meset->insert(item)))
                            unique_counters[0].increment();
                        if (std::get<1>(enamine->insert(item)))
                            enamine_unique_counters[0].increment();
                    }
                }
            }, &stopWriter, &dbase, &dbase_enamine);

    // Consumers
    for (size_t i = 1; i != n_threads; ++i) {
        threads.emplace_back(
                [&](MutexCounter *total_counter, MutexCounter *valid_counter,
                        std::atomic<bool> *stop) {
                    std::string item;
                    boost::optional<std::string> value;
                    bool itemsLeft;
                    do {
                        // It's important to fence (if the producers have finished) *before* dequeueing
                        itemsLeft = stop->load(std::memory_order_acquire);
                        while (q.try_dequeue(item)) {
                            itemsLeft = true;
                            task(item, total_counter, valid_counter, &q_out);
                        }
                    } while (itemsLeft || doneConsumers.fetch_add(1, std::memory_order_acq_rel) + 1 == (n_threads - 1));
                }, &total_counters[i],
                   &valid_counters[i],
                   &doneProducer);
    }

    // Wait for all threads
    for (size_t i = 0; i != n_threads; ++i) {
        threads[i].join();
    }

    std::cout << "ok everyone else finsihed. Waiting for monitar" << std::endl;
    stopMonitar = true;
    monitar.join();

    std::cout << "ok everyone else finsihed. Waiting for Writer" << std::endl;
    stopWriter = true;
    writer.join();

    // Collect any leftovers (could be some if e.g. consumers finish before producers)
    std::string item;
    boost::optional<std::string> value;
    while (q.try_dequeue(item)) {
        task(item, &total_counters[0], &valid_counters[0], &q_out);
    }

    while (q_out.try_dequeue(item)) {
        if (std::get<1>(dbase.insert(item)))
            unique_counters[0].increment();
        if (std::get<1>(dbase_enamine.insert(item)))
            enamine_unique_counters[0].increment();
    }

    int unique = 0;
    int total = 0;
    int valid = 0;
    int enamine =0;
    for (size_t i = 0; i < n_threads; i++) {
        unique += unique_counters[i].view();
        total += total_counters[i].view();
        valid += valid_counters[i].view();
        enamine += enamine_unique_counters[i].view();
    }

    std::cout << "total,unique,valid,emiane" << std::endl;
    std::cout << total << ", " << unique << ", " << valid << "," << enamine << std::endl;
    std::cout << dbase.size() << std::endl;
    std::cout << dbase_enamine.size() << std::endl;

}
