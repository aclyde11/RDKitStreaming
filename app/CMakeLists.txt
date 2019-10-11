

#set(Boost_NO_SYSTEM_PATHS ON)
#set(Boost_NO_BOOST_CMAKE TRUE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" )
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")


#set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_MULTITHREADED OFF)
#set(Boost_USE_STATIC_RUNTIME ON)
find_package( Boost COMPONENTS iostreams filesystem system)


find_package( MyRDKit COMPONENTS DataStructs Fingerprints SmilesParse Depictor
        RDGeneral GraphMol )


find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)

set( LIBS ${RDKIT_LIBRARIES} Boost::iostreams ${RDKit_THREAD_LIBS} z  )

include_directories(${RDKIT_INCLUDE_DIR})

add_library(nop INTERFACE)
target_include_directories(nop INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/app>
        $<INSTALL_INTERFACE:include/nop>
        )


add_library(parallel_hashmap INTERFACE)
target_include_directories(parallel_hashmap INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/app>
        $<INSTALL_INTERFACE:include/parallel_hashmap>
        )

add_library(concurrentqueue INTERFACE)
target_include_directories(concurrentqueue INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/app>
        $<INSTALL_INTERFACE:include/concurrentqueue>
        )

add_library(StreamingMolecularReader INTERFACE)
target_include_directories(StreamingMolecularReader INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/app>
        $<INSTALL_INTERFACE:include/StreamingMolecularReader>
        )

find_package(cityhash)
include_directories(${CITYHASH_INCLUDE_DIR})


include_directories("${CMAKE_SOURCE_DIR}/include/")

add_executable(example1 example1.cpp)
target_link_libraries(example1 nop)
target_link_libraries(example1 parallel_hashmap)
target_link_libraries(example1 concurrentqueue)
target_link_libraries(example1 ${LIBS} )
target_link_libraries(example1 ${CITYHASH_LIBRARIES})
