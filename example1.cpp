#include <iostream>
#include <fstream>
#include <iostream>
#include <bitset>
#include <cinttypes>
#include <string>
#include <atomic>
#include <chrono>
#include <unordered_map>

#include <boost/optional.hpp>

#include "concurrentqueue.h"

#include "parallel_hashmap/phmap.h"

#include "nop/serializer.h"
#include "nop/utility/die.h"
#include "nop/utility/stream_reader.h"
#include "nop/utility/stream_writer.h"

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>

#define THREADS 6

using parallel_smile_set = phmap::parallel_flat_hash_map<std::string, int>;

namespace {

// Sends fatal errors to std::cerr.
    auto Die() { return nop::Die(std::cerr); }

}

struct MutexCounter {
    std::mutex lock;
    int value;

    void increment() {
//        lock.lock();
        value += 1;
//        lock.unlock();
    }

    int view() {
        int tmp;
        lock.lock();
        tmp = value;
        lock.unlock();
        return tmp;
    }
};

boost::optional<std::string> getSmileCannon(std::string const &smi) {
    RDKit::ROMol *mol1;
    try {
        mol1 = RDKit::SmilesToMol(smi);
    } catch (...) {
        return {};
    }

    if (mol1 != nullptr) {
        auto tmp = RDKit::MolToSmiles(*mol1);
        delete mol1;
        return tmp;
    } else {
        return {};
    }
}


boost::optional<std::string> consumeItem(std::string const &item) {
    return getSmileCannon(item);
}

void monitarThread(MutexCounter *valid_counters,
                   MutexCounter *unique_counters,
                   MutexCounter *total_counters,  int monitar_freq, std::atomic<bool> *stop) {

    while(!(stop->load(std::memory_order_acquire))) {
        std::this_thread::sleep_for(std::chrono::seconds(5));
        int total = 0;
        int valid = 0;
        int unique = 0;

        for(int i = 0; i < THREADS; i++) {
            total += total_counters[i].view();
            valid += valid_counters[i].view();
            unique += unique_counters[i].view();
        }

        std::cout << total << ", " << valid << ", " << unique << std::endl;
    }
}



std::unordered_map<std::string, int> getInitalSet(std::string const &filename) {
    std::unordered_map<std::string, int> email{2000000};

    moodycamel::ConcurrentQueue<std::string> q;
    std::thread threads[THREADS];

    std::atomic<bool> doneProducer(true);
    std::atomic<int> doneConsumers(0);

    std::cout << "startin size, max size" << email.bucket_count() << " " << std::endl;

    std::cout << "Reading in train.txt" << std::endl;

    // Producers
    // This guy is reading in from std in and add it to the queues.
    threads[0] = std::thread([&](std::atomic<bool> *stop) {
        std::ifstream ifs(filename);
        int counter =0;
        for (std::string line; std::getline(ifs, line);) {
            q.enqueue(line);
            counter++;
            if (counter >= 100000) {
                *stop = false;
                return;
            }
        }
        *stop = false;
    }, &doneProducer);

    std::mutex lock;
    for (int i = 1; i != THREADS; ++i) {
        threads[i] = std::thread(
                [&](
                        std::atomic<bool> *stop, std::unordered_map<std::string, int> *meset, std::mutex *lock) {
                    std::string item;
                    boost::optional<std::string> value;
                    bool itemsLeft;
                    do {
                        // It's important to fence (if the producers have finished) *before* dequeueing
                        itemsLeft = stop->load(std::memory_order_acquire);
                        while (q.try_dequeue(item)) {
                            itemsLeft = true;
                            value = consumeItem(item);
                            if (value.has_value()) {
                                lock->lock();
                                meset->insert({value.value(),1});
                                lock->unlock();
                            }
                        }
                    } while (itemsLeft || doneConsumers.fetch_add(1, std::memory_order_acq_rel) + 1 == THREADS
                                                                                                       - 1);
                },
                &doneProducer, &email, &lock);
    }

    // Wait for all threads
    for (int i = 0; i != THREADS;
         ++i) {
        threads[i].join();
    }

    // Collect any leftovers (could be some if e.g. consumers finish before producers)
    std::string item;
    boost::optional<std::string> value;
    while (q.try_dequeue(item)) {
        value = consumeItem(item);
        if (value.has_value()) {
            email.insert({value.value(),1});
        }
    }
    std::cout << "done." << std::endl;

    return email;
}


int main(int argc, char **argv) {
    bool LOAD = (bool)(atoi(argv[1]));
    // read in smiles from worker queue

    // parallel test.
    moodycamel::ConcurrentQueue<std::string> q;
    std::thread threads[THREADS];

    std::atomic<bool> doneProducer(true);
    std::atomic<int> doneConsumers(0);

    MutexCounter total_counters[THREADS] = {};
    MutexCounter valid_counters[THREADS] = {};
    MutexCounter unique_counters[THREADS] = {};

    auto start = std::chrono::high_resolution_clock::now();
    std::unordered_map<std::string, int> email_;
    parallel_smile_set email;

    if (LOAD) {
        email_ = getInitalSet("/Users/austin/train.txt");
        using Writer = nop::StreamWriter<std::ofstream>;
        nop::Serializer<Writer> serializer{"example.txt"};
        serializer.Write(email_) || Die();
        std::cout << "wrote" << std::endl;
        serializer.writer().stream().close();
    } else {
        using Reader = nop::StreamReader<std::ifstream>;
        nop::Deserializer<Reader> deserializer{"example.txt"};
        deserializer.Read(&email_) || Die();
        std::cout << email_.size() << std::endl;
    }

     //convert email_ to what I want
    for (std::pair<std::string, int> element : email_)
    {
        std::string st = element.first;
        int a = element.second;
        email.insert({st, a});
    }

    std::cout << "new email " << email.size() <<std::endl;

    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    std::cout << "inital loading took " << microseconds << std::endl;

    std::cout << email.size() << std::endl;
    std::cout << "Moving on to your problem." << std::endl;
    //producers
    // This guy is reading in from std in and add it to the queues.
    threads[0] = std::thread([&](std::atomic<bool> *stop) {
        for (std::string line; std::getline(std::cin, line);) {
            q.enqueue(line);
        }
        *stop = false;
    }, &doneProducer);

    std::atomic<bool> stopMonitar(false);
    //monitar

    std::thread monitar(
            [&](MutexCounter *valid_counters,
                MutexCounter *unique_counters,
                MutexCounter *total_counters,  std::atomic<bool> *stop) {

                monitarThread(valid_counters, unique_counters, total_counters, 10, stop);
            }, valid_counters, unique_counters, total_counters, &stopMonitar);

    // Consumers
    for (int i = 1; i != THREADS; ++i) {
        threads[i] = std::thread(
                [&](MutexCounter *total_counter, //std::mutex *total_counter_lock,
                         MutexCounter *valid_counter, //std::mutex *valid_coutner_lock,
                         MutexCounter *unqiue_counter, //std::mutex *unqiue_counter_lock,
                         std::atomic<bool> *stop, parallel_smile_set *meset) {
                    std::string item;
                    boost::optional<std::string> value;
                    bool itemsLeft;
                    do {
                        // It's important to fence (if the producers have finished) *before* dequeueing
                        itemsLeft = stop->load(std::memory_order_acquire);
                        while (q.try_dequeue(item)) {
                            itemsLeft = true;
                            value = consumeItem(item);

                            total_counter->increment();

                            if (value.has_value()) {
                                valid_counter->increment();
                                if (std::get<1>(meset->insert({value.value(), 1})))
                                    unqiue_counter->increment();
                            }
                        }
                    } while (itemsLeft || doneConsumers.fetch_add(1, std::memory_order_acq_rel) + 1 == THREADS - 1);
                }, &total_counters[i], //&total_counter_locks[i],
                   &valid_counters[i], //&valid_counter_locks[i],
                   &unique_counters[i], //&unique_counter_locks[i],
                   &doneProducer, &email);
    }

    // Wait for all threads
    for (int i = 0; i != THREADS;
         ++i) {
        threads[i].join();
    }

    std::cout << "ok everyone else finsihed. Waiting for monitar" << std::endl;
    stopMonitar = true;
    monitar.join();

    // Collect any leftovers (could be some if e.g. consumers finish before producers)
    std::string item;
    boost::optional<std::string> value;
    while (q.try_dequeue(item)) {
        value = consumeItem(item);
        total_counters[0].increment();
        if (value.has_value()) {
            if (std::get<1>(email.insert({value.value(),1})))
                valid_counters[0].increment();
        }
    }

    int unique = 0;
    int total = 0;
    int valid = 0;
    for (int i = 1; i < THREADS; i++) {
        unique += unique_counters[i].view();
        total += total_counters[i].view();
        valid += valid_counters[i].view();

    }

    std::cout << total << ", " << unique << ", " << valid << std::endl;
    std::cout << email.size() << std::endl;
}