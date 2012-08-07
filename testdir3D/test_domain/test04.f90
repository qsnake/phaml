program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=1,triangle_files="brickcyl.geo")
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     print_error_when=FINAL,print_error_what=ENERGY_LINF_ERR, &
                     print_error_who=MASTER,max_vert=1000)
call phaml_destroy(soln)
end program phaml_master
