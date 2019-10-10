//
// Created by Austin Clyde on 2019-10-06.
//

#ifndef RDKITSV_PARALLELDATABASE_H
#define RDKITSV_PARALLELDATABASE_H

#include "rdkit_interface.h"

#include <thread>
#include <vector>
#include <string>
#include <atomic>
#include <unordered_map>
#include <boost/optional.hpp>
#include <concurrentqueue/concurrentqueue.h>
#include <city.h>

namespace SMR {

    /**
     * This struct is only used by a single writer and a single reader. Locking is not needed for writing.
     */
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

    struct CityHasher {
        std::size_t operator()(const std::string& k) const {
            return CityHash64(k.data(), k.size());
        }
    };

    using myStdMap = std::unordered_map<std::string, bool, CityHasher>;

    /**
     * Generic WorkQueue Function
     */
    myStdMap getInitialSetFromFile( size_t n_threads = 4,
                                                                size_t init_size = 2000000) {
        myStdMap email{init_size};
        std::mutex lock;

        moodycamel::ConcurrentQueue<std::string> q;
        moodycamel::ConcurrentQueue<std::string> q_out;

        std::vector<std::thread> threads;
        std::atomic<bool> doneProducer(true);
        std::atomic<size_t> doneConsumers(0);

        // Producers
        threads.emplace_back([&](std::atomic<bool> *stop) {
            for (std::string line; std::getline(std::cin, line);) {
                if (q.size_approx() < 500000) {
                    q.enqueue(line);
                } else {
                    std::this_thread::sleep_for(std::chrono::seconds(1));
                    q.enqueue(line);
                }
            }
            *stop = false;
        }, &doneProducer);

        for (size_t i = 1; i < n_threads; ++i) {
            threads.emplace_back(
                    [&](std::atomic<bool> *stop, myStdMap *meset, std::mutex *lock) {
                        bool itemsLeft;
                        do {
                            itemsLeft = stop->load(std::memory_order_acquire);
                            std::string item;
                            while (q.try_dequeue(item)) {
                                itemsLeft = true;
                                //work here:
                                boost::optional<std::string> value = getCannonicalSmileFromSmile(item);
                                if (value.has_value()) {
                                    q_out.enqueue(value.value());
                                }
                            }
                        } while (itemsLeft ||
                                 doneConsumers.fetch_add(1, std::memory_order_acq_rel) + 1 == n_threads - 1);
                    }, &doneProducer, &email, &lock);
        }

        //Monitar Thread
        std::atomic<bool> stopMonitar(false);
        std::thread monitar(
                [&](moodycamel::ConcurrentQueue<std::string> *q, std::atomic<bool> *stop) {

                    std::cout << "time,queuesize,total,valid,unique" << std::endl;
                    auto start_time = std::time(nullptr);
                    while (!(stop->load(std::memory_order_acquire))) {
                        std::this_thread::sleep_for(std::chrono::seconds(10));

                        std::cout << std::time(nullptr) - start_time << "," << q->size_approx() << ", " << q_out.size_approx() << ", " << email.size() << std::endl;
                    }
                }, &q, &stopMonitar);

        //Writer Thread
        std::atomic<bool> stopWriter(false);
        std::thread writer(
                [&]( std::atomic<bool> *stop, myStdMap *meset) {

                    std::cout << "time,queuesize,total,valid,unique" << std::endl;
                    while (!(stop->load(std::memory_order_acquire))) {
                        std::string item;
                        while (q_out.try_dequeue(item)) {
                            meset->insert({item, true});
                        }
                    }
                }, &stopWriter, &email);

        // Wait for all threads
        for (size_t i = 0; i != n_threads; ++i) {
            threads[i].join();
        }

        stopMonitar = true;
        monitar.join();

        stopWriter = true;
        writer.join();

        // Collect any leftovers (could be some if e.g. consumers finish before producers)
        std::string item;
        boost::optional<std::string> value;
        while (q.try_dequeue(item)) {
            value = getCannonicalSmileFromSmile(item);
            if (value.has_value()) {
                email.insert({value.value(), true});
            }
        }

        return email;
    }


    /**
     * Generic WorkQueue Function
     */
    myStdMap getInitialSetFromFile(std::string const &filename, size_t n_threads = 4,
                                                                size_t init_size = 2000000) {
        myStdMap email{init_size};
        std::mutex lock;

        moodycamel::ConcurrentQueue<std::string> q;
        moodycamel::ConcurrentQueue<std::string> q_out;

        std::vector<std::thread> threads;
        std::atomic<bool> doneProducer(true);
        std::atomic<size_t> doneConsumers(0);

        // Producers
        threads.emplace_back([&](std::atomic<bool> *stop) {
            std::ifstream ifs(filename);
            for (std::string line; std::getline(ifs, line);) {
                if (q.size_approx() < 500000) {
                    q.enqueue(line);
                } else {
                    std::this_thread::sleep_for(std::chrono::seconds(1));
                    q.enqueue(line);
                }
            }
            *stop = false;
        }, &doneProducer);

        for (size_t i = 1; i < n_threads; ++i) {
            threads.emplace_back(
                    [&](std::atomic<bool> *stop, myStdMap *meset, std::mutex *lock) {
                        bool itemsLeft;
                        do {
                            itemsLeft = stop->load(std::memory_order_acquire);
                            std::string item;
                            while (q.try_dequeue(item)) {
                                itemsLeft = true;
                                //work here:
                                boost::optional<std::string> value = getCannonicalSmileFromSmile(item);
                                if (value.has_value()) {
                                    q_out.enqueue(value.value());
                                }
                            }
                        } while (itemsLeft ||
                                 doneConsumers.fetch_add(1, std::memory_order_acq_rel) + 1 == n_threads - 1);
                    }, &doneProducer, &email, &lock);
        }

        //Monitar Thread
        std::atomic<bool> stopMonitar(false);
        std::thread monitar(
                [&](moodycamel::ConcurrentQueue<std::string> *q, std::atomic<bool> *stop) {

                    std::cout << "time,queuesize,total,valid,unique" << std::endl;
                    auto start_time = std::time(nullptr);
                    while (!(stop->load(std::memory_order_acquire))) {
                        std::this_thread::sleep_for(std::chrono::seconds(10));

                        std::cout << std::time(nullptr) - start_time << "," << q->size_approx() << ", " << q_out.size_approx() << ", " << email.size() << std::endl;
                    }
                }, &q, &stopMonitar);

        //Writer Thread
        std::atomic<bool> stopWriter(false);
        std::thread writer(
                [&]( std::atomic<bool> *stop, myStdMap *meset) {

                    std::cout << "time,queuesize,total,valid,unique" << std::endl;
                    while (!(stop->load(std::memory_order_acquire))) {
                        std::string item;
                        while (q_out.try_dequeue(item)) {
                            meset->insert({item, true});
                        }
                    }
                }, &stopWriter, &email);

        // Wait for all threads
        for (size_t i = 0; i != n_threads; ++i) {
            threads[i].join();
        }

        stopMonitar = true;
        monitar.join();

        stopWriter = true;
        writer.join();

        // Collect any leftovers (could be some if e.g. consumers finish before producers)
        std::string item;
        boost::optional<std::string> value;
        while (q.try_dequeue(item)) {
            value = getCannonicalSmileFromSmile(item);
            if (value.has_value()) {
                email.insert({value.value(), true});
            }
        }

        return email;
    }

}

#endif //RDKITSV_PARALLELDATABASE_H
