program phaml_master
use phaml
use global
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=1,triangle_files="polyhole", &
                  max_blen=.25_my_real,id=3)
call phaml_solve_pde(soln,max_vert=2000,sequential_vert=10)
open(unit=11,file="savepolyhole.poly")
call phaml_store_grid(soln,11,GRIDFILE_POLY)
close(11)
call phaml_destroy(soln)
end program phaml_master
