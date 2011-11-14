program phaml_master
use phaml
use phaml_user_mod
implicit none
type(phaml_solution_type) :: soln
real(my_real), pointer :: x(:), y(:), u(:)
integer :: i
call phaml_create(soln,nproc=1)
call update_usermod(soln)
call phaml_solve_pde(soln,max_refsolveloop=1,mg_cycles=20)
call phaml_get_grid_soln(soln,x,y,u)
do i=1,size(x)
   write(6,"(SS,1P,3E19.12E2)") x(i),y(i),u(i)
end do
deallocate(x,y,u)
call phaml_destroy(soln)
end program phaml_master
