# General parameters
cmake_minimum_required(VERSION 2.8.12)
project (splines)

# Add boost dependencies
find_package(Boost)

# Add eigen dependencies
find_path(EIGEN3_INCLUDE_DIRS "Eigen" PATHS "${SDK_DIR}/include" PATH_SUFFIXES "eigen3")
if (EIGEN3_INCLUDE_DIRS)
  include_directories ("${EIGEN3_INCLUDE_DIRS}")
endif (EIGEN3_INCLUDE_DIRS)

add_definitions(-std=c++11)
include_directories(include ${Boost_INCLUDE_DIRS})

file (GLOB SOURCE "src/*.cpp")
add_library (splines SHARED ${SOURCE})
target_link_libraries(splines ${Boost_LIBRARIES})