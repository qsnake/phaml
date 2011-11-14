include_directories(../../src)
add_executable(elliptic master.f90 usermod.f90 pde.f90)
target_link_libraries(elliptic phaml arpack lapack blas)
