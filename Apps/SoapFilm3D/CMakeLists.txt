append_files (Headers "h" .)
append_files (Sources "cpp" .)

include_directories (${EIGEN_INCLUDE_DIR})

add_executable (SoapFilm3D ${Headers} ${Sources})
target_link_libraries (SoapFilm3D LosTopos)

