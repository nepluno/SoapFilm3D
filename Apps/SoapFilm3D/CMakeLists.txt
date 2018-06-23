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

include_directories (${EIGEN_INCLUDE_DIR})

add_executable (SoapFilm3D ${Headers} ${Sources})
target_link_libraries (SoapFilm3D LosTopos ${DEFAULT_LIBRARIES})

