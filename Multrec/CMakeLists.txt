cmake_minimum_required(VERSION 2.8)

project(Multrec)



set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++0x")


set(SOURCES
    main.cpp
        trees/genespeciestreeutil.cpp
        trees/newicklex.cpp
        trees/node.cpp
        trees/treeinfo.cpp
        trees/treeiterator.cpp
        multigenereconciler.cpp
)


include_directories(.)


add_executable(Multrec ${SOURCES} )
