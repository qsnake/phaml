program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,                   &
                     max_vert=500,           &
                     mg_cycles=5,            &
                     print_warnings=.false., &
                     print_grid_when=PHASES, &
                     print_grid_who=MASTER  ,&
                     print_linsys_when=PHASES, &
                     print_linsys_who=MASTER, &
                     print_solver_who=MASTER, &
                     print_solver_when=PHASES, &
                     print_error_when=TOO_MUCH  , &
                     print_error_what=ENERGY_ERR, &
                     print_errest_what=ENERGY_ERREST, &
                     print_error_who=MASTER)
call phaml_destroy(soln)
end program phaml_master
