cmake_minimum_required (VERSION 3.1)

project (SoapFilm3D)
set (CMAKE_CXX_STANDARD 11)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# add directory with macros
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# add cmake file
include (common)

# set options

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

# Eigen library is required
find_package (Eigen3 REQUIRED)
include_directories (${EIGEN3_INCLUDE_DIR})

option (USE_PNG "Builds in support for exporting to png" ON)
if (USE_PNG)
  find_package (PNG)
  if (PNG_FOUND)
    add_definitions (-DHAVE_PNG)
    add_definitions (${PNG_DEFINITIONS})
    include_directories (${PNG_INCLUDE_DIR})
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${PNG_LIBRARIES})
  endif (PNG_FOUND)
endif (USE_PNG)

find_package (LAPACK REQUIRED)
if (LAPACK_FOUND)
add_definitions (-DHAVE_LAPACK)
include_directories (${LAPACK_INCLUDE_DIR})
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${LAPACK_LIBRARIES})
endif (LAPACK_FOUND)

find_package (BLAS REQUIRED)
add_definitions (-DHAVE_BLAS)
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${BLAS_LIBRARIES})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/LosTopos3D)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/newsparse)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/tunicate)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Apps/SoapFilm3D/fmmtl)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory (LosTopos)

add_subdirectory (Apps)

