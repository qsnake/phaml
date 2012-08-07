program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     print_error_when=PHASES,print_error_what=ENERGY_L2_ERR, &
                     print_error_who=MASTER, &
                     max_eq=5000,reftype=HP_ADAPTIVE, &
                     hp_strategy=HP_COEF_ROOT)
call phaml_destroy(soln)
end program phaml_master
