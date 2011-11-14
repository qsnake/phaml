# compile elliptic example
include_directories(../../src)
add_executable(master master.f90 usermod.f90 pde.f90)
target_link_libraries(master phaml arpack lapack blas)
