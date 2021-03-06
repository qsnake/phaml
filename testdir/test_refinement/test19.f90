program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     max_eq=500,mg_cycles=5)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     max_refsolveloop=1, mg_cycles=5, &
                     refterm=KEEP_NEQ_SMOOTH)
call phaml_destroy(soln)
end program phaml_master
