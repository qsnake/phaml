include_directories(../../src)
add_executable(simple master.f90 pde.f90)
target_link_libraries(simple phaml arpack lapack blas)
