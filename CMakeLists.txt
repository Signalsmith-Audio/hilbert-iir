cmake_minimum_required(VERSION 3.24)

add_library(signalsmith-hilbert INTERFACE)
target_include_directories(signalsmith-hilbert INTERFACE ${CMAKE_CURRENT_LIST_DIR}/include)
