cmake_minimum_required( VERSION 3.1 )

set (CMAKE_CXX_STANDARD 17)

project(RDKitSV)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(Boost_NO_SYSTEM_PATHS ON)
#set(BOOST_ROOT:PATHNAME /usr/local/Cellar/boost/1.71.0/)
set(Boost_NO_BOOST_CMAKE TRUE)

set(BOOST_ROOT /usr/local/Cellar/boost/1.71.0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" )

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
find_package( Boost COMPONENTS iostreams filesystem system)


find_package( MyRDKit COMPONENTS ChemReactions FileParsers Fingerprints SmilesParse Depictor RDGeometryLib
        RDGeneral SubstructMatch Subgraphs MolDraw2D GraphMol DistGeometry DistGeomHelpers
        MolAlign Optimizer ForceField ForceFieldHelpers Alignment
        ForceField  MolTransforms EigenSolvers)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})


find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)

set( LIBS ${RDKIT_LIBRARIES} Boost::iostreams ${RDKit_THREAD_LIBS} z  )


add_library(cereal INTERFACE)
target_include_directories(cereal INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        )


include_directories(${RDKIT_INCLUDE_DIR})
foreach( exnum RANGE 1 1 )
    message( "Example : ${exnum}" )
    set( exe "example${exnum}" )
    set( src "${exe}.cpp" )
    add_executable( ${exe} ${src} concurrentqueue.h)
    target_include_directories(${exe} PUBLIC .)
    target_include_directories(${exe} PRIVATE "${CMAKE_SOURCE_DIR}/parallel_hashmap")
    target_link_libraries( ${exe} ${LIBS} )
endforeach( exnum )

