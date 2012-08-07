!---------------------------------------------------------------------!
!                                PHAML                                !
!                                                                     !
! The Parallel Hierarchical Adaptive MultiLevel code for solving      !
! linear elliptic partial differential equations of the form          !
! (PUx)x + (QUy)y + RU = F on 2D polygonal domains with mixed         !
! boundary conditions, and eigenvalue problems where F is lambda*U.   !
!                                                                     !
! PHAML is public domain software.  It was produced as part of work   !
! done by the U.S. Government, and is not subject to copyright in     !
! the United States.                                                  !
!                                                                     !
!     William F. Mitchell                                             !
!     Applied and Computational Mathematics Division                  !
!     National Institute of Standards and Technology                  !
!     william.mitchell@nist.gov                                       !
!     http://math.nist.gov/phaml                                      !
!                                                                     !
!---------------------------------------------------------------------!

!----------------------------------------------------
! This file contains the main program supplied by the user.
! This is a simple example that just solves one linear elliptic pde.
!----------------------------------------------------

!       ------------
program phaml_master
!       ------------

!----------------------------------------------------
! This is the main program.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use global
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Local variables

type(phaml_solution_type) :: soln
character(len=32) :: tfiles
!----------------------------------------------------
! Begin executable code

! pick which domain to use by uncommenting one of these lines

! tfiles = "tet.geo"
! tfiles = "unitcube.geo"
 tfiles = "unitcube.msh"
! tfiles = "parpipe.geo"
! tfiles = "cylinder.geo"
! tfiles = "cylinder.geomsh"
! tfiles = "sphere.geo"
! tfiles = "brickcyl.geo"

call phaml_create(soln,nproc=1,draw_grid_who=MASTER, &
!                  spawn_form=DEBUG_SLAVE,debug_command="idbc", &
                  triangle_files=tfiles)

! for the tetrahedron domain, refinement stalls unless we do a few uniform
! refinements first

if (tfiles == "tet.geo") then
   call phaml_solve_pde(soln,reftype=H_UNIFORM,max_refsolveloop=4)
endif

call phaml_solve_pde(soln,                                  &
                     max_vert=10000,                        &
                     draw_grid_when=PHASES    ,             &
                     pause_after_phases=.true.,             &
                     pause_at_start=.true.,                 &
                     petsc_rtol = 1.0e-4_my_real,           &
                     print_grid_when=PHASES,                &
                     print_grid_who=MASTER,                 &
                     print_solver_who=MASTER,               &
                     print_solver_when=PHASES,              &
                     errtype=RELATIVE_ERROR,                &
                     print_error_when=PHASES,               &
                     print_error_what=ENERGY_LINF_ERR,      &
                     print_errest_what=ENERGY_LINF_ERREST,  &
                     print_error_who=MASTER)

call phaml_destroy(soln)

end program phaml_master
