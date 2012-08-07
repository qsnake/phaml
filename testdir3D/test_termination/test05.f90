program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=FINAL,print_grid_who=MASTER, &
                     print_error_when=FINAL,print_error_who=MASTER, &
                     print_errest_what=ENERGY_ERREST, &
                     term_energy_err=7.0e-1_my_real)
call phaml_destroy(soln)
end program phaml_master
