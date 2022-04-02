cmake_minimum_required (VERSION 3.1)

project (SoapFilm3D)
set (CMAKE_CXX_STANDARD 11)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# add directory with macros
set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# add cmake file
include (common)

# set options

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

# Eigen library is required
find_package (Eigen3 REQUIRED)
include_directories (${EIGEN3_INCLUDE_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/LosTopos3D)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/newsparse)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/tunicate)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Apps/SoapFilm3D/fmmtl)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory (LosTopos)

add_subdirectory (Apps)

