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
!     Mathematical and Computational Sciences Division                !
!     National Institute of Standards and Technology                  !
!     william.mitchell@nist.gov                                       !
!     http://math.nist.gov/phaml                                      !
!                                                                     !
!---------------------------------------------------------------------!

module stack_mod

!----------------------------------------------------
! This module implements stacks.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use message_passing, only: fatal
!----------------------------------------------------

implicit none
private
public create_stack, destroy_stack, push_stack, pop_stack, peek_stack, &
       size_stack, END_OF_STACK, stack_integer, stack_my_real

!----------------------------------------------------
! The following types are defined:

type stack_integer
   integer :: p = 0
   integer, pointer :: val(:) => NULL()
end type stack_integer

type stack_my_real
   integer :: p = 0
   real(my_real), pointer :: val(:) => NULL()
end type stack_my_real

!----------------------------------------------------
! The following defined constants are defined:

integer, parameter :: END_OF_STACK = -huge(0)
!----------------------------------------------------
! The following generic interfaces are defined:

interface create_stack
   module procedure create_stack_integer, create_stack_my_real
end interface

interface destroy_stack
   module procedure destroy_stack_integer, destroy_stack_my_real
end interface

interface push_stack
   module procedure push_stack_integer, push_stack_my_real
end interface

interface pop_stack
   module procedure pop_stack_integer, pop_stack_my_real
end interface

interface peek_stack
   module procedure peek_stack_integer, peek_stack_my_real
end interface

interface size_stack
   module procedure size_stack_integer, size_stack_my_real
end interface

!----------------------------------------------------

contains

subroutine create_stack_integer(stack)
type(stack_integer), intent(inout) :: stack
if (associated(stack%val)) call destroy_stack_integer(stack)
allocate(stack%val(8))
stack%p = 0
end subroutine create_stack_integer

subroutine create_stack_my_real(stack)
type(stack_my_real), intent(inout) :: stack
if (associated(stack%val)) call destroy_stack_my_real(stack)
allocate(stack%val(8))
stack%p = 0
end subroutine create_stack_my_real

subroutine destroy_stack_integer(stack)
type(stack_integer), intent(inout) :: stack
if (associated(stack%val)) deallocate(stack%val)
end subroutine destroy_stack_integer

subroutine destroy_stack_my_real(stack)
type(stack_my_real), intent(inout) :: stack
if (associated(stack%val)) deallocate(stack%val)
end subroutine destroy_stack_my_real

subroutine push_stack_integer(stack,val)
type(stack_integer), intent(inout) :: stack
integer, intent(in) :: val
integer, allocatable :: temp(:)
integer :: astat
if (.not. associated(stack%val)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("push_stack: stack has not been created, or has been destroyed")
   stop
endif
if (stack%p == size(stack%val)) then
   allocate(temp(size(stack%val)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("push_stack: allocation failed")
      stop
   endif
   temp = stack%val
   deallocate(stack%val)
   allocate(stack%val(2*size(temp)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("push_stack: allocation failed")
      stop
   endif
   stack%val(1:size(temp)) = temp
   deallocate(temp)
endif
stack%p = stack%p + 1
stack%val(stack%p) = val
end subroutine push_stack_integer

subroutine push_stack_my_real(stack,val)
type(stack_my_real), intent(inout) :: stack
real(my_real), intent(in) :: val
real(my_real), allocatable :: temp(:)
integer :: astat
if (.not. associated(stack%val)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("push_stack: stack has not been created, or has been destroyed")
   stop
endif
if (stack%p == size(stack%val)) then
   allocate(temp(size(stack%val)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("push_stack: allocation failed")
      stop
   endif
   temp = stack%val
   deallocate(stack%val)
   allocate(stack%val(2*size(temp)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("push_stack: allocation failed")
      stop
   endif
   stack%val(1:size(temp)) = temp
   deallocate(temp)
endif
stack%p = stack%p + 1
stack%val(stack%p) = val
end subroutine push_stack_my_real

function pop_stack_integer(stack)
type(stack_integer), intent(inout) :: stack
integer :: pop_stack_integer
if (.not. associated(stack%val)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("pop_stack: stack has not been created, or has been destroyed")
   stop
endif
if (stack%p == 0) then
   pop_stack_integer = END_OF_STACK
else
   pop_stack_integer = stack%val(stack%p)
   stack%p = stack%p - 1
endif
end function pop_stack_integer

function pop_stack_my_real(stack)
type(stack_my_real), intent(inout) :: stack
real(my_real) :: pop_stack_my_real
if (.not. associated(stack%val)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("pop_stack: stack has not been created, or has been destroyed")
   stop
endif
if (stack%p == 0) then
   pop_stack_my_real = END_OF_STACK
else
   pop_stack_my_real = stack%val(stack%p)
   stack%p = stack%p - 1
endif
end function pop_stack_my_real

function peek_stack_integer(stack)
type(stack_integer), intent(in) :: stack
integer :: peek_stack_integer
if (.not. associated(stack%val)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("peek_stack: stack has not been created, or has been destroyed")
   stop
endif
if (stack%p == 0) then
   peek_stack_integer = END_OF_STACK
else
   peek_stack_integer = stack%val(stack%p)
endif
end function peek_stack_integer

function peek_stack_my_real(stack)
type(stack_my_real), intent(in) :: stack
real(my_real) :: peek_stack_my_real
if (.not. associated(stack%val)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("peek_stack: stack has not been created, or has been destroyed")
   stop
endif
if (stack%p == 0) then
   peek_stack_my_real = END_OF_STACK
else
   peek_stack_my_real = stack%val(stack%p)
endif
end function peek_stack_my_real

function size_stack_integer(stack)
type(stack_integer), intent(in) :: stack
integer :: size_stack_integer
size_stack_integer = stack%p
end function size_stack_integer

function size_stack_my_real(stack)
type(stack_my_real), intent(in) :: stack
integer :: size_stack_my_real
size_stack_my_real = stack%p
end function size_stack_my_real

end module stack_mod
