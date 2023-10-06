# cmake_minimum_required(VERSION 3.12)
project(pmetParallel)

# Set compilation parameters
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

# set out dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_executable(pmetParallel
  main.cpp
  Output.cpp
  motif.cpp
  motifComparison.cpp
  utils.cpp
  Output.hpp
  motif.hpp
  motifComparison.hpp
  utils.hpp
)

# Add performance testing tools
# target_link_libraries(${PROJECT_NAME} profiler unwind)
target_link_libraries(${PROJECT_NAME})