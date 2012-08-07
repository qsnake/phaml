program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: pde1
call phaml_create(pde1,nproc=4,eq_type=EIGENVALUE)
call phaml_solve_pde(pde1,print_grid_when=PHASES,print_grid_who=MASTER, &
                     print_warnings=.false.,print_error_who=MASTER, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_eval_when=PHASES,print_eval_who=MASTER, &
                     print_solver_when=PHASES,print_solver_who=MASTER, &
                     max_vert=500, &
                     solver=PETSC_GMRES_SOLVER,  &
                     preconditioner=PETSC_SOR_PRECONDITION, &
                     eigensolver=SLEPC_JACOBI_DAVIDSON, &
                     num_eval=1)
call phaml_destroy(pde1)

end program phaml_master
