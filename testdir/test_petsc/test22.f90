program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     print_error_when=PHASES,print_error_who=MASTER, &
                     print_error_what=ENERGY_LINF_L2_ERR,max_vert=500, &
                     petsc_matrix_free=.true., &
                     petsc_rtol=1.0e-10_my_real, &
                     solver=PETSC_GMRES_SOLVER, &
                     preconditioner=MG_PRECONDITION,mg_cycles=1)
call phaml_destroy(soln)
end program phaml_master
