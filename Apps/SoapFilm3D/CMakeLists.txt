append_files (Headers "h" .)
append_files (Sources "cpp" .)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.45.0) 

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS}) 
  set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${Boost_LIBRARIES})
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Locate GLUT
find_package (GLUT REQUIRED glut)
if (GLUT_FOUND)
  include_directories (${GLUT_INCLUDE_DIR})
  set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${GLUT_glut_LIBRARY})
else (GLUT_FOUND)
  message (SEND_ERROR "Unable to locate GLUT")
endif (GLUT_FOUND)

find_package(GLEW REQUIRED)
if (GLEW_FOUND)
  include_directories(${GLEW_INCLUDE_DIRS})
  set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${GLEW_LIBRARIES})
else (GLEW_FOUND)
  message (SEND_ERROR "Unable to locate GLEW")
endif()

option(NO_SHADER "Build Without Shader" OFF)

if(NO_SHADER)
add_definitions(-DNO_SHADER)
endif()

include_directories (${EIGEN_INCLUDE_DIR})

add_executable (SoapFilm3D ${Headers} ${Sources})
target_link_libraries (SoapFilm3D LosTopos ${DEFAULT_LIBRARIES})

