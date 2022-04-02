add_definitions(-DUSE_FORTRAN_BLAS)

set (Directories
 .
 common
 common/tunicate
 common/newsparse
 LosTopos3D
)

append_files (Headers "h" ${Directories})
append_files (Sources "cpp" ${Directories})

find_package (OpenGL REQUIRED)
include_directories (${OPENGL_INCLUDE_DIR})

find_package (GLUT REQUIRED glut)
include_directories (${GLUT_INCLUDE_DIR})

find_package (Eigen3 REQUIRED)
include_directories (${EIGEN3_INCLUDE_DIR})

add_library(LosTopos STATIC ${Headers} ${Sources})
target_link_libraries (LosTopos ${DEFAULT_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})


