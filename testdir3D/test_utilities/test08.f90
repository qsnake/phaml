program phaml_master
use phaml
use global
implicit none

type(phaml_solution_type) :: soln

call phaml_create(soln,nproc=1)

call phaml_solve_pde(soln,max_vert=200,sequential_vert=10)

open(unit=11,file="savemsh.msh")
call phaml_store_grid(soln,11,GRIDFILE_MSH)
close(11)

call phaml_destroy(soln)

end program phaml_master
