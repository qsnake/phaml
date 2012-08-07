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

module refine_adapt_mod

!----------------------------------------------------
! This module contains routines for the basic control of adaptive refinement.
!
! communication tags in this module are of the form 33xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use gridtype_mod
use grid_util
use linsystype_mod
use hp_strategies
use refine_elements
use message_passing
use hash_mod
use linsys_io
use load_balance
use error_estimators
use omp_lib
use zoltan_interf
!----------------------------------------------------

implicit none
private
public refine_adaptive, use_old_refinement, reconcile

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!        ------------------
function use_old_refinement(refine_control,procs,predictive)
!        ------------------

!----------------------------------------------------
! This routine whether we can use the new approach to the overall refinement
! scheme, or if this is one of the cases where we still need the old approach.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(refine_options), intent(in) :: refine_control
type(proc_info), intent(in) :: procs
logical, intent(in) :: predictive
logical :: use_old_refinement
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! Assume we can use the new code.

use_old_refinement = .false.

! Uniform refinement has it's own code.

if (refine_control%reftype == P_UNIFORM .or. &
    refine_control%reftype == H_UNIFORM) use_old_refinement = .true.

! The reference solution and NLP strategies have their own refinement code.

if (refine_control%reftype == HP_ADAPTIVE .and. &
    (refine_control%hp_strategy == HP_REFSOLN_EDGE .or. &
     refine_control%hp_strategy == HP_REFSOLN_ELEM .or. &
     refine_control%hp_strategy == HP_NLP)) use_old_refinement = .true.

! TEMP methods that may perform more than one refinememnt/derefinement at
!      a time.  Allowing these requires modification to mark_for_refinement.

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_T3S) use_old_refinement = .true.

! TEMP STEEPEST_SLOPE can change p by more than 1, and do both h and change p

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_STEEPEST_SLOPE) use_old_refinement = .true.

end function use_old_refinement

!          ---------------
subroutine refine_adaptive(grid,procs,refine_control,solver_control,io_control,&
                           lb,still_sequential,init_nvert,init_nelem,init_dof, &
                           loop,partition_method,balance_what,predictive, &
                           stalled,no_time)
!          ---------------

!----------------------------------------------------
! This routine performs adaptive refinement of the grid.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), target, intent(in) :: procs
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
type(io_options), intent(in) :: io_control
type(Zoltan_Struct), pointer :: lb
integer, intent(in) :: init_nvert, init_nelem, init_dof, loop, &
                       partition_method, balance_what
logical, intent(in) :: still_sequential, predictive
logical, intent(out) :: stalled
logical, intent(in), optional :: no_time
!----------------------------------------------------
! Local variables:

logical :: any_changes, any_p_coarsened, any_h_coarsened, any_p_refined, &
           any_h_refined
real(my_real) :: target
integer :: nelem, loop_count, astat, itemp
integer, pointer :: desired_level(:),desired_degree(:), &
                    leaf_elements(:)
!----------------------------------------------------
! Begin executable code


if (use_old_refinement(refine_control,procs,predictive)) then
   call old_refine_adaptive(grid,procs,refine_control,solver_control, &
                            io_control,still_sequential,init_nvert,init_nelem, &
                            init_dof,loop,balance_what,predictive,no_time)
   stalled = .false.
   return
endif

any_h_coarsened = .false.
any_p_coarsened = .false.
any_h_refined   = .false.
any_p_refined   = .false.

! Allocate memory.

allocate(desired_level(size(grid%element)),desired_degree(size(grid%element)), &
         leaf_elements(grid%biggest_elem),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in refine_adaptive")
   stop
endif

! Recompute error indicators, if they're stale

if (.not. grid%errind_up2date) then
   call all_error_indicators(grid,refine_control%error_estimator)
endif

! Set the target that determines when enough refinement has occured.

target = set_termination_target(grid,procs,refine_control,init_nvert, &
                                init_nelem,init_dof,loop,still_sequential)

! Repeat until enough refinement has occured.

loop_count = 0
do
   loop_count = loop_count + 1

! Make a list of leaf elements.

   if (size(leaf_elements) /= grid%biggest_elem) then
      deallocate(leaf_elements)
      allocate(leaf_elements(grid%biggest_elem),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in refine_adaptive")
         stop
      endif
   endif

   call list_elements(grid,leaf_elements,nelem,leaf=.true.)

! Mark elements for refinement.  This sets, for each leaf element, the new
! h-level and p in the conforming adapted grid.

   call mark_for_refinement(grid,procs,refine_control,still_sequential, &
                            leaf_elements,nelem,desired_level,desired_degree, &
                            any_changes,loop_count<4)

! If there are no changes on any processor, quit.

   if (.not. any_changes) exit

! Perform predictive load balancing.

   if (predictive .and. .not. still_sequential) then
      call pred_load_balance(grid,procs,refine_control,solver_control,lb, &
                             partition_method,still_sequential,balance_what, &
                             desired_level,desired_degree,any_changes)

! Reallocate if predictive load balance reallocated grid%element or performed
! refinement that changed the biggest element lid

      if (size(desired_level) /= size(grid%element)) then
         deallocate(desired_level,desired_degree)
         allocate(desired_level(size(grid%element)), &
                  desired_degree(size(grid%element)), stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in refine_adaptive")
            stop
         endif
      endif

      if (size(leaf_elements) /= grid%biggest_elem) then
         deallocate(leaf_elements)
         allocate(leaf_elements(grid%biggest_elem),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in refine_adaptive")
            stop
         endif
      endif

! redo mark for refinement because load balance changed grid

      if (any_changes) then
         call list_elements(grid,leaf_elements,nelem,leaf=.true.)
         call mark_for_refinement(grid,procs,refine_control,still_sequential, &
                                  leaf_elements,nelem,desired_level, &
                                  desired_degree,any_changes,loop_count<4)
      endif
   endif

! Stop doing coarsening after 3 loops (arbitrary 3) to avoid an infinite
! loop of coarsen/refine.

   if (loop_count < 4) then

      if (refine_control%derefine) then

! Perform coarsenings, first p then h.

         if (refine_control%reftype == P_UNIFORM .or. &
             refine_control%reftype == P_ADAPTIVE .or. &
             refine_control%reftype == HP_ADAPTIVE) then
            call p_coarsen_grid(grid,refine_control,desired_degree, &
                                leaf_elements,nelem,any_p_coarsened)
         else
            any_p_coarsened = .false.
         endif
         if (refine_control%reftype == H_UNIFORM .or. &
             refine_control%reftype == H_ADAPTIVE .or. &
             refine_control%reftype == HP_ADAPTIVE) then
            call h_coarsen_grid(grid,leaf_elements,nelem,desired_level, &
                                desired_degree,refine_control,any_h_coarsened)
         else
            any_h_coarsened = .false.
         endif

! If not increasing the size of the grid, check to see if we're done after
! the coarsening phase.

         if (refine_control%refterm == KEEP_NVERT        .or. &
             refine_control%refterm == KEEP_NELEM        .or. &
             refine_control%refterm == KEEP_NEQ          .or. &
             refine_control%refterm == KEEP_ERREST) then
            if (adapt_done(grid,procs,refine_control,still_sequential, &
                           target,.true.)) then
               any_p_refined = .false.
               any_h_refined = .false.
               exit
            endif
         endif

      else

         any_p_coarsened = .false.
         any_h_coarsened = .false.

      endif

   else
      any_p_coarsened = .false.
      any_h_coarsened = .false.

   endif

! Remake the list of leaf elements because h coarsening invalidates it.

   if (any_h_coarsened) then
      call list_elements(grid,leaf_elements,nelem,leaf=.true.)
   endif

! Perform refinements, first p then h.

   if (refine_control%reftype == P_UNIFORM .or. &
       refine_control%reftype == P_ADAPTIVE .or. &
       refine_control%reftype == HP_ADAPTIVE) then
      call p_refine_grid(grid,refine_control,desired_degree,leaf_elements, &
                         nelem,any_p_refined)
   else
      any_p_refined = .false.
   endif
   if (refine_control%reftype == H_UNIFORM .or. &
       refine_control%reftype == H_ADAPTIVE .or. &
       refine_control%reftype == HP_ADAPTIVE) then
      call h_refine_grid(grid,refine_control,solver_control,desired_level, &
                         desired_degree,leaf_elements,nelem,any_h_refined)
   else
      any_h_refined = .false.
   endif

! See if enough refinement has occured.

   if (adapt_done(grid,procs,refine_control,still_sequential,target, &
                  any_p_coarsened .or. any_h_coarsened .or. &
                  any_p_refined .or. any_h_refined)) exit

end do

! check for stalling

if (loop_count == 1) then
   if (any_p_coarsened .or. any_h_coarsened .or. any_p_refined .or. &
       any_h_refined) then
      itemp = 1
   else
      itemp = 0
   endif
   itemp = phaml_global_max(procs,itemp,3307)
   if (itemp == 0) then
      stalled = .true.
   else
      stalled = .false.
   endif
else
   stalled = .false.
endif

! Free memory.

deallocate(desired_level,desired_degree,leaf_elements)

! Enforce the overlap conditions around the partition boundary.

if (num_proc(procs) > 1 .and. .not. still_sequential) then
   call enforce_overlap(grid,refine_control,solver_control,procs)
endif

end subroutine refine_adaptive

!        ----------------------
function set_termination_target(grid,procs,refine_control,init_nvert, &
                                init_nelem,init_dof,loop,still_sequential) &
                                result(target)
!        ----------------------

!----------------------------------------------------
! This function determines the value to use as the target for terminating
! refinement.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: init_nvert, init_nelem, init_dof, loop
logical, intent(in) :: still_sequential
real(my_real) :: target
!----------------------------------------------------
! Local variables:

integer :: itarget
!----------------------------------------------------
! Begin executable code

! special cases

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_T3S) then
   if (refine_control%t3s_reftype == H_UNIFORM) then
      target = 0
      return
   endif
endif

! Set the target value depending on what type of termination is used.

select case (refine_control%refterm)

case (DOUBLE_NVERT)

   target = init_nvert*refine_control%inc_factor**loop
   if (target > refine_control%max_vert) target = refine_control%max_vert

case (DOUBLE_NELEM)

   target = init_nelem*refine_control%inc_factor**loop
   if (target > refine_control%max_elem) target = refine_control%max_elem

case (DOUBLE_NEQ)

   target = init_dof*refine_control%inc_factor**loop
   if (target > refine_control%max_dof) target = refine_control%max_dof

case (HALVE_ERREST)

   if (.not. grid%errind_up2date) then
      call all_error_indicators(grid,refine_control%error_estimator)
   endif
   target = compute_global_max_errind(grid,procs,still_sequential)/2

case (KEEP_NVERT)

   if (refine_control%max_vert /= huge(0)) then
      target = refine_control%max_vert
   else
      call get_grid_info(grid,procs,still_sequential,3301,total_nvert=itarget)
      target = itarget
      if (target > refine_control%max_vert) target = refine_control%max_vert
   endif

case (KEEP_NELEM)

   if (refine_control%max_elem /= huge(0)) then
      target = refine_control%max_elem
   else
      call get_grid_info(grid,procs,still_sequential,3301, &
                         total_nelem_leaf=itarget)
      target = itarget
      if (target > refine_control%max_elem) target = refine_control%max_elem
   endif

case (KEEP_NEQ)

   if (refine_control%max_dof /= huge(0)) then
      target = refine_control%max_dof
   else
      call get_grid_info(grid,procs,still_sequential,3301,total_dof=itarget)
      target = itarget
      if (target > refine_control%max_dof) target = refine_control%max_dof
   endif

case (KEEP_ERREST)

! TEMP not doing KEEP_ERREST
   call warning("have not yet decided how to handle refterm==KEEP_ERREST")
   target = 0

case (ONE_REF, ONE_REF_HALF_ERRIND)

   target = 0

case default

   call fatal("illegal value for refterm",procs=procs)
   stop

end select

end function set_termination_target

!        ----------
function adapt_done(grid,procs,refine_control,still_sequential,target, &
                    any_changes)
!        ----------

!----------------------------------------------------
! This routine checks to see if the target for terminating refinement is met.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
logical, intent(in) :: still_sequential
real(my_real), intent(in) :: target
logical, intent(in) :: any_changes
logical :: adapt_done
!----------------------------------------------------
! Local variables:

integer :: compare
!----------------------------------------------------
! Begin executable code

! Assume we're not done.

adapt_done = .false.

! Compare the refterm-dependent value to the target.

select case (refine_control%refterm)

case (DOUBLE_NVERT, KEEP_NVERT)

   call get_grid_info(grid,procs,still_sequential,3302,total_nvert=compare)
   if (compare >= target) adapt_done = .true.

case (DOUBLE_NELEM, KEEP_NELEM)

   call get_grid_info(grid,procs,still_sequential,3302,total_nelem_leaf=compare)
   if (compare >= target) adapt_done = .true.

case (DOUBLE_NEQ, KEEP_NEQ)

   call get_grid_info(grid,procs,still_sequential,3302,total_dof=compare)
   if (compare >= target) adapt_done = .true.

case (HALVE_ERREST, KEEP_ERREST)

   if (.not. grid%errind_up2date) then
      call all_error_indicators(grid,refine_control%error_estimator)
   endif
   if (compute_global_max_errind(grid,procs,still_sequential) <= target) then
      adapt_done = .true.
   endif

case (ONE_REF, ONE_REF_HALF_ERRIND)

   adapt_done = .true.

case default

   call fatal("illegal value for refterm")
   stop

end select

! check for reaching maximum levels or degree

if (.not. adapt_done) then
   if (refine_control%stop_on_maxlev) then
      call get_grid_info(grid,procs,still_sequential,3304,max_nlev=compare)
      if (compare >= refine_control%max_lev) adapt_done = .true.
   endif

   if (refine_control%stop_on_maxdeg) then
      call get_grid_info(grid,procs,still_sequential,3305,maxdeg=compare)
      if (compare >= refine_control%max_deg) adapt_done = .true.
   endif
endif

if (.not. adapt_done) then
   if (any_changes) then
      compare = 1
   else
      compare = 0
   endif
   compare = phaml_global_max(procs,compare,3306)
   if (compare == 0) adapt_done = .true.
endif

end function adapt_done

!          -------------------
subroutine mark_for_refinement(grid,procs,refine_control,still_sequential, &
                               leaf_elements,nleaf,desired_level, &
                               desired_degree,any_changes,coarsen_OK)
!          -------------------

!----------------------------------------------------
! This routine determines the new h-level and p for each element in a
! conforming adapted grid.
! TEMP this is only for refinement strategies that do a single refine/coarsen
!      by h or p
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
logical, intent(in) :: still_sequential
integer, intent(in) :: leaf_elements(:)
integer, intent(in) :: nleaf
integer, intent(out) :: desired_level(:), desired_degree(:)
logical, intent(out) :: any_changes
logical, intent(in) :: coarsen_OK
!----------------------------------------------------
! Local variables:

real(my_real) :: refine_cutoff, coarsen_cutoff
integer :: i, elem, new_p(2)
character(len=1) :: reftype(grid%biggest_elem)
logical :: h_coarsen, p_coarsen
!----------------------------------------------------
! Begin executable code

! Initialize so there is something in the unused entries.

desired_level = 0
desired_degree = 0
any_changes = .false.

! Update the error indicators.

call all_error_indicators(grid,refine_control%error_estimator,delayed=.true.)

! Determine the error indicator cutoffs for whether or not an element
! should be refined or coarsened.

call determine_cutoffs(grid,procs,refine_control,still_sequential, &
                       leaf_elements,nleaf,refine_cutoff,coarsen_cutoff)

! In 3D, if the refine cutoff is essentially 0, we have probably refined all
! elements, because the initial guess is setting the hierarchical coefficients
! to 0, so quit because we don't know what to refine.

if (global_element_kind == TETRAHEDRAL_ELEMENT .and. &
    refine_cutoff < 1.0e-100_my_real) return

! For each leaf element ...

!$omp parallel do &
!$omp  default(shared) &
!$omp  private(i,elem,new_p,h_coarsen,p_coarsen) &
!$omp  reduction(.or.: any_changes)

do i=1,nleaf
   elem = leaf_elements(i)

! If I don't own it, set the level to 1 to be set by the owner or
! compatibility, and leave the degree the same.

   if (.not. grid%element(elem)%iown) then

      desired_level(elem) = 1
      desired_degree(elem) = grid%element(elem)%degree

! If the error indicator for any eigenvector is sufficiently large, refine.

   elseif (any(grid%element_errind(elem,:) >= refine_cutoff)) then

! Determine if refinement is by h or p.

      call mark_reftype_one(elem,grid,refine_control,-1.0_my_real,reftype, &
                            .false.,new_p)

! Set the result accordingly.

      select case (reftype(elem))

      case ("h")
         desired_level(elem) = min(refine_control%max_lev, &
                                   grid%element(elem)%level+1)
         desired_degree(elem) = grid%element(elem)%degree

      case ("p")
         desired_level(elem) = grid%element(elem)%level
         desired_degree(elem) = min(refine_control%max_deg, &
                                    grid%element(elem)%degree+1)

      case ("n")
         desired_level(elem) = grid%element(elem)%level
         desired_degree(elem) = grid%element(elem)%degree

      case default
         ierr = PHAML_INTERNAL_ERROR
         call fatal("mark_reftype_one did not return h, p or n in mark_for_refinement")
         stop

      end select

! If not refining, determine if it should be h or p coarsened and set result.

   else

      if (coarsen_OK) then
         call determine_coarsening(grid,refine_control,elem,coarsen_cutoff, &
                                   h_coarsen,p_coarsen)
      else
         h_coarsen = .false.
         p_coarsen = .false.
      endif

      if (h_coarsen) then
         desired_level(elem) = max(1,grid%element(elem)%level-1)
         desired_degree(elem) = grid%element(elem)%degree
      elseif (p_coarsen) then
         desired_level(elem) = grid%element(elem)%level
         desired_degree(elem) = max(1,grid%element(elem)%degree-1)
      else
         desired_level(elem) = grid%element(elem)%level
         desired_degree(elem) = grid%element(elem)%degree
      endif
   endif

   if (desired_level(elem) /= grid%element(elem)%level .or. &
       desired_degree(elem) /= grid%element(elem)%degree) then
      any_changes = .true.
   endif

end do
!$omp end parallel do

! Reconcile the conflicting desired refinements/coarsenings to get a
! compatible grid.

call reconcile_requests(grid,procs,still_sequential,leaf_elements,nleaf, &
                        desired_level,desired_degree,any_changes)

end subroutine mark_for_refinement

!          -----------------
subroutine determine_cutoffs(grid,procs,refine_control,still_sequential, &
                             leaf_elements,nleaf,refine_cutoff,coarsen_cutoff)
!          -----------------

!----------------------------------------------------
! This routine determines the error indicator cutoff values for determining
! whether or not an element should be refined or coarsened.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
logical, intent(in) :: still_sequential
integer, intent(in) :: leaf_elements(:)
integer, intent(in) :: nleaf
real(my_real), intent(out) :: refine_cutoff, coarsen_cutoff
!----------------------------------------------------
! Local variables:

real(my_real), parameter :: fourth_root_2 = 1.1892071150027_my_real
real(my_real) :: global_max_errind, normsoln
integer :: nelem
!----------------------------------------------------
! Begin executable code

! The cutoffs depend on the refinement termination criterion.

select case(refine_control%refterm)

case(ONE_REF)

! For ONE_REF, use the given reftol divided by the square root of the number
! of elements, possibly scaled by the norm of the solution, for refinement.

   call get_grid_info(grid,procs,still_sequential,3303,total_nelem_leaf=nelem)
   refine_cutoff = refine_control%reftol/sqrt(real(nelem,my_real))

! TEMP only using first component of first eigenvalue for solution norm
   if (grid%errtype == RELATIVE_ERROR .and. &
       .not. (refine_control%reftype == HP_ADAPTIVE .and. &
       (refine_control%hp_strategy == HP_T3S .or. &
        refine_control%hp_strategy == HP_ALTERNATE))) then
      call norm_solution(grid,procs,still_sequential,1,1,energy=normsoln)
      if (normsoln /= 0.0_my_real) refine_cutoff = refine_cutoff*normsoln
   endif

! Set the coarsening cutoff at 1/100 the refinement cutoff.

   coarsen_cutoff = refine_cutoff/100

case(ONE_REF_HALF_ERRIND)

! For ONE_REF_HALF_ERRIND, the refine cutoff is the maximum error indicator
! divided by inc_factor.

! TEMP120127 if this works, add leaf_elements to other calls to compute_glob...
   global_max_errind = compute_global_max_errind(grid,procs,still_sequential, &
                                                 leaf_elements,nleaf)
   refine_cutoff = global_max_errind/refine_control%inc_factor

! Set the coarsening cutoff at 1/100 the maximum error indicator.

   coarsen_cutoff = global_max_errind/100

case default

! For all others, it's the maximum error indicator divided by the fourth root
! of 2, which is the bin width in the old code.

   global_max_errind = compute_global_max_errind(grid,procs,still_sequential)
   refine_cutoff = global_max_errind/fourth_root_2

! Set the coarsening cutoff at 1/100 the maximum error indicator.

   coarsen_cutoff = global_max_errind/100

end select

end subroutine determine_cutoffs

!          --------------------
subroutine determine_coarsening(grid,refine_control,elem,cutoff, &
                                h_coarsen,p_coarsen)
!          --------------------

!----------------------------------------------------
! This routine determines if an element should be h or p coarsened.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: elem
real(my_real), intent(in) :: cutoff
logical, intent(out) :: h_coarsen, p_coarsen
!----------------------------------------------------
! Local variables:

integer :: parent, siblings(4)
logical :: hcoarsen_possible, pcoarsen_possible
real(my_real) :: hcoarsen_errind, pcoarsen_errind
!----------------------------------------------------
! Begin executable code

! relatives

if (grid%element(elem)%level /= 1) then
   parent = hash_decode_key(grid%element(elem)%gid/MAX_CHILD,grid%elem_hash)
   siblings(1:2) = get_child_lid(grid%element(parent)%gid,ALL_CHILDREN, &
                                 grid%elem_hash)
   if (grid%element(parent)%mate == BOUNDARY) then
      siblings(3:4) = NO_CHILD
   else
      siblings(3:4) = get_child_lid( &
                     grid%element(parent)%mate,ALL_CHILDREN,grid%elem_hash)
   endif
endif

! Check conditions that allow h coarsening.
! derefine must be true.
! Cannot h coarsen level 1 elements.
! Only h coarsen elements that I own.
! Only do h coarsening with h or hp adaptive refinement.

hcoarsen_possible = refine_control%derefine .and. grid%element(elem)%iown &
                    .and. grid%element(elem)%level > 1
if (refine_control%reftype /= H_ADAPTIVE .and. &
    refine_control%reftype /= HP_ADAPTIVE) hcoarsen_possible = .false.

! Check conditions that allow p coarsening.
! derefine must be true.
! Cannot p coarsen an element that has p==1.
! Only p coarsen elements that I own.
! Only do p coarsening with p or hp adaptive.

pcoarsen_possible = refine_control%derefine .and. grid%element(elem)%iown &
                    .and. grid%element(elem)%degree > 1
if (refine_control%reftype /= P_ADAPTIVE .and. &
    refine_control%reftype /= HP_ADAPTIVE) pcoarsen_possible = .false.

! Compute coarsening error indicators

if (hcoarsen_possible) then
   hcoarsen_errind = hcoarsen_indicator(grid,parent,siblings)
else
   hcoarsen_errind = huge(0.0_my_real)
endif
if (pcoarsen_possible) then
   pcoarsen_errind = pcoarsen_indicator(grid,elem)
else
   pcoarsen_errind = huge(0.0_my_real)
endif

! If either coarsening indicator is small enough, perform the type of
! coarsening that has the least impact on the solution

h_coarsen = .false.
p_coarsen = .false.
if (hcoarsen_errind < pcoarsen_errind .and. &
    hcoarsen_errind < cutoff) then
   h_coarsen = .true.
elseif (pcoarsen_errind < cutoff) then
   p_coarsen = .true.
endif

end subroutine determine_coarsening

!          ------------------
subroutine reconcile_requests(grid,procs,still_sequential,leaf_elements,nleaf, &
                              desired_level,desired_degree,any_changes)
!          ------------------

!----------------------------------------------------
! This routine reconciles the desired new h-levels and p's so that the new
! desired levels and p's will result in a compatible grid across all processors
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: leaf_elements(:), nleaf
integer, intent(inout) :: desired_level(:),desired_degree(:)
logical, intent(inout) :: any_changes
!----------------------------------------------------
! Local variables:

integer :: part_bound_elements(grid%biggest_elem), npartbound, isend2(2), &
           new_desired(size(desired_level)), i, elem, j, neigh, p, &
           neighbors(NEIGHBORS_PER_ELEMENT,grid%biggest_elem), &
           neighneigh(NEIGHBORS_PER_ELEMENT), neigh_change, loop, &
           nirecv, nrrecv
logical :: changed, changed_this_pass
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
!----------------------------------------------------
! Begin executable code

! Make a list of owned leaf elements that have an unowned neighbor.

if (num_proc(procs) > 1 .and. .not. still_sequential) then
   call list_elements(grid,part_bound_elements,npartbound,own=.true., &
                      leaf=.true.,unowned_neigh=.true.)
else
   npartbound = 0
endif

! Exchange the desired p and h levels of partition-boundary elements with
! the neighbors.  The p's are needed so that the edges on the partition
! boundary get the right degree, but only have to be passed in this first
! message.

changed = .true.
loop = 0
if (num_proc(procs) > 1 .and. .not. still_sequential) then
   call exchange_desired(loop,grid,procs,changed,part_bound_elements, &
                         npartbound,any_changes,desired_level,desired_degree)
endif

! The 3D adaptive refinement algorithm does not maintain compatibility during
! individual refinements, so does not need to match up desired_level.

if (global_element_kind == TETRAHEDRAL_ELEMENT) return

! Identify the neighbors of each leaf element

!$omp parallel do &
!$omp  default(shared) &
!$omp  private(i,elem)

do i=1,nleaf
   elem = leaf_elements(i)
   neighbors(:,elem) = get_neighbors(elem,grid)
end do

!$omp end parallel do

! Repeat until no processor made a change in its desireds.

new_desired = desired_level

do while (changed)
   loop = loop + 1

! Repeat until there are no more changes to my desireds.

   changed = .false.

   do

      changed_this_pass = .false.

! For each leaf element ...

!$omp parallel do &
!$omp  default(shared) &
!$omp  private(i, elem, j, neigh, neighneigh, neigh_change) &
!$omp  reduction(.or.: changed_this_pass)

      do i=1,nleaf
         elem = leaf_elements(i)

! For each neighbor of elem ...

         do j=1,NEIGHBORS_PER_ELEMENT
            neigh = neighbors(j,elem)
            if (neigh == BOUNDARY) cycle

! Get the neighbors of the neighbor, to see if elem is across neigh's base.

            neighneigh = neighbors(:,neigh)

! See how much the desired level of the neighbor differs from its current level.

            neigh_change = desired_level(neigh)-grid%element(neigh)%level

! Increase the desired level of elem, if necessary to make a compatible grid.
! By always increasing, this gives a preference to refining over coarsening.

            if (neigh_change < 0) then ! coarsening of neighbor

               if (neighneigh(3) == elem) then ! shares base
                  if (desired_level(neigh) > new_desired(elem) .and. &
                      desired_level(neigh) > desired_level(elem)) then
                     new_desired(elem) = desired_level(neigh)
                     changed_this_pass = .true.
                  endif
               else ! does not share base
                  if (desired_level(neigh) > new_desired(elem) .and. &
                      desired_level(neigh) > desired_level(elem)) then
                     new_desired(elem) = desired_level(neigh)
                     changed_this_pass = .true.
                  endif
               endif

            elseif (neigh_change == 0) then ! neighbor stays the same

               if (neighneigh(3) == elem) then ! shares base
                  if (desired_level(neigh) > new_desired(elem)-1 .and. &
                      desired_level(neigh)-1 > desired_level(elem)) then
                     new_desired(elem) = desired_level(neigh)-1
                     changed_this_pass = .true.
                  endif
               else ! does not share base
                  if (desired_level(neigh) > new_desired(elem) .and. &
                      desired_level(neigh) > desired_level(elem)) then
                     new_desired(elem) = desired_level(neigh)
                     changed_this_pass = .true.
                  endif
               endif

            elseif (neigh_change == 1) then ! neighbor refined once

               if (neighneigh(3) == elem) then ! shares base
                  if (desired_level(neigh) > new_desired(elem)-1 .and. &
                      desired_level(neigh)-1 > desired_level(elem)) then
                     new_desired(elem) = desired_level(neigh)-1
                     changed_this_pass = .true.
                  endif
               else ! does not share base
                  if (desired_level(neigh) > new_desired(elem)-1 .and. &
                      desired_level(neigh)-1 > desired_level(elem)) then
                     new_desired(elem) = desired_level(neigh)-1
                     changed_this_pass = .true.
                  endif
               endif

            elseif (neigh_change == 2*(neigh_change/2)) then ! neighbor refined
                                                      ! an even number of times

               if (neighneigh(3) == elem) then ! shares base
                  if (desired_level(neigh) > new_desired(elem)-2 .and. &
                      desired_level(neigh)-2 > desired_level(elem)) then
                     new_desired(elem) = desired_level(neigh)-2
                     changed_this_pass = .true.
                  endif
               else ! does not share base
                  if (desired_level(neigh) > new_desired(elem)-1 .and. &
                      desired_level(neigh)-1 > desired_level(elem)) then
                     new_desired(elem) = desired_level(neigh)-1
                     changed_this_pass = .true.
                  endif
               endif

            else ! neighbor refined an odd number of times but more than once

               if (neighneigh(3) == elem) then ! shares base
                  if (desired_level(neigh) > new_desired(elem)-1 .and. &
                      desired_level(neigh)-1 > desired_level(elem)) then
                     new_desired(elem) = desired_level(neigh)-1
                     changed_this_pass = .true.
                  endif
               else ! does not share base
                  if (desired_level(neigh) > new_desired(elem)-2 .and. &
                      desired_level(neigh)-2 > desired_level(elem)) then
                     new_desired(elem) = desired_level(neigh)-2
                     changed_this_pass = .true.
                  endif
               endif

            endif

         end do ! neighbors
      end do ! leaf elements
!$omp end parallel do

! Check for any changes this time through the loop, and repeat if there are.

      if (.not. changed_this_pass) exit
      changed = .true.
      desired_level = new_desired
   end do

! Exchange the desired h levels of partition-boundary elements with
! the neighbors.

   if (num_proc(procs) > 1 .and. .not. still_sequential) then
      call exchange_desired(loop,grid,procs,changed,part_bound_elements, &
                            npartbound,any_changes,desired_level)
   endif

! inform the master about whether there are any changes

   if (my_proc(procs) == 1) then
      isend2 = 0
      if (changed) isend2(1) = 1
      if (any_changes) isend2(2) = 1
      call phaml_send(procs,MASTER,isend2,2,(/0.0_my_real/),0,3330+loop)
   elseif (my_proc(procs) == MASTER) then
      call phaml_recv(procs,p,irecv,nirecv,rrecv,nrrecv,3330+loop)
      changed = irecv(1) == 1
      any_changes = irecv(2) == 1
      deallocate(irecv)
   endif

end do ! until no processor made a change in its desireds

end subroutine reconcile_requests

!          ----------------
subroutine exchange_desired(loop,grid,procs,changed,part_bound_elements, &
                            npartbound,any_changes,desired_level,desired_degree)
!          ----------------

!----------------------------------------------------
! This routine exchanges the desired p and h levels of partition-boundary
! elements with the neighbors.  On input, if changed is false then there
! were no changes to desired_level on this processor so the data does not need
! to be passed.  On output, if changed is false then it is false on all
! processors, so we are done reconciling requests.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: loop
type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
logical, intent(inout) :: changed
integer, intent(in) :: part_bound_elements(:), npartbound
logical, intent(inout) :: any_changes
integer, intent(inout) :: desired_level(:)
integer, intent(inout), optional :: desired_degree(:)
!----------------------------------------------------
! Local variables:

integer :: my_processor, nproc, nisend, count, i, p, ip, nirecv, nrrecv, lid, &
           k, neigh(NEIGHBORS_PER_ELEMENT), astat, istart
logical :: on_boundary
integer, allocatable :: isend(:), nrecv(:)
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
!----------------------------------------------------
! Begin executable code

! Convenience variables.

my_processor = my_proc(procs)
nproc = num_proc(procs)

! Allocate message.

nisend = 3 + (1+KEY_SIZE)*npartbound
if (present(desired_degree)) then
   nisend = nisend + npartbound
endif
allocate(isend(nisend),nrecv(nproc),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in exchange_desired")
   stop
endif

! Pack the changed flag and number of elements.

if (changed) then
   isend(1) = 1
else
   isend(1) = 0
endif
if (any_changes) then
   isend(2) = 1
else
   isend(2) = 0
endif
isend(3) = npartbound

! For each element in the list of partition-boundary elements, pack the
! element ID and desired values.

if (changed) then
   count = 4
   do i=1,npartbound
      call hash_pack_key(grid%element(part_bound_elements(i))%gid,isend,count)
      count = count + KEY_SIZE
      isend(count) = desired_level(part_bound_elements(i))
      count = count + 1
      if (present(desired_degree)) then
         isend(count) = desired_degree(part_bound_elements(i))
         count = count + 1
      endif
   end do
else
   isend(4:nisend) = 0
endif

! Exchange the message with the other processors.

! TEMP eventually send it only to neighbor processors; need to be careful
!      about global reduction of changed in that case

if (my_processor /= MASTER) then
   call phaml_alltoall(procs,isend,nisend,irecv,nrecv,3310+loop)

! Receive and process the messages.

   istart = 0
   do ip=1,nproc
      if (ip == my_processor) then
         istart = istart + nrecv(ip)
         cycle
      endif

! Indicate if there was a change

      if (irecv(1+istart) == 1) changed = .true.
      if (irecv(2+istart) == 1) any_changes = .true.

! If there are changes in this message, desired should be the maximum of
! what the owner of each element thinks and what I have, if the element is
! on the partition boundary.
! TEMP I'm not sure if this disrupts coarsening.

      if (irecv(1+istart) == 1) then
         count = 4+istart
         do i=1,irecv(3+istart) ! irecv(3) is the number of elements sent
            lid = hash_decode_key(hash_unpack_key(irecv,count),grid%elem_hash)
            if (lid == HASH_NOT_FOUND) then
               count = count + KEY_SIZE+1
               if (present(desired_degree)) count = count + 1
            else
               neigh = get_neighbors(lid,grid)
               on_boundary = .false.
               do k=1,NEIGHBORS_PER_ELEMENT
                  if (neigh(k) /= boundary) then
                     if (grid%element(neigh(k))%iown) on_boundary = .true.
                  endif
               end do
               count = count + KEY_SIZE
               if (on_boundary) then
                  desired_level(lid) = max(desired_level(lid),irecv(count))
                  count = count + 1
                  if (present(desired_degree)) then
                     desired_degree(lid) = max(desired_degree(lid),irecv(count))
                     count = count + 1
                  endif
               else
                  count = count + 1
                  if (present(desired_degree)) count = count + 1
               endif
            endif
         end do
      end if
      istart = istart + nrecv(ip)
   end do
   deallocate(irecv)
endif

end subroutine exchange_desired

!          --------------
subroutine p_coarsen_grid(grid,refine_control,desired_degree,elem_list,nelem, &
                          any_changes)
!          --------------

!----------------------------------------------------
! This routine p-coarsens elements in elem_list whose current degree is larger
! than the desired degree.  NOTE: When this routine returns, the edge degree
! rule is not satisfied.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: desired_degree(:), elem_list(:), nelem
logical, intent(out) :: any_changes
!----------------------------------------------------
! Local variables:

integer :: i, errcode, delta_dof, delta_dof_own, total_delta_dof, &
           total_delta_dof_own
!----------------------------------------------------
! Begin executable code

any_changes = .false.

! p-coarsen the interior of each element that wants to be p coarsened.

total_delta_dof = 0
total_delta_dof_own = 0

!$omp parallel do schedule(dynamic) &
!$omp  default(shared) &
!$omp  private(i,errcode,delta_dof,delta_dof_own) &
!$omp  reduction(+ : total_delta_dof, total_delta_dof_own) &
!$omp  reduction(.or.: any_changes)

do i=1,nelem
   do while (grid%element(elem_list(i))%degree > desired_degree(elem_list(i)))
      call p_coarsen_elem_interior(grid,elem_list(i),errcode,refine_control, &
                                   delta_dof,delta_dof_own,delay_errind=.true.)
      if (errcode == 0) any_changes = .true.
      total_delta_dof = total_delta_dof + delta_dof
      total_delta_dof_own = total_delta_dof_own + delta_dof_own
   end do
end do
!$omp end parallel do

! do things that must be done OpenMP-sequentially after the parallel loop

grid%dof = grid%dof + total_delta_dof
grid%dof_own = grid%dof_own + total_delta_dof_own

end subroutine p_coarsen_grid

!          --------------
subroutine h_coarsen_grid(grid,leaf_list,nleaf,desired_level,desired_degree, &
                          refine_control,any_changes)
!          --------------

!----------------------------------------------------
! This routine h-coarsens elements in elem_list whose current level is larger
! than the desired level.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(inout) :: leaf_list(:), nleaf, desired_level(:), &
                          desired_degree(:)
type(refine_options) :: refine_control
logical, intent(out) :: any_changes
!----------------------------------------------------
! Local variables:

integer :: i, lev, coarsen_list(grid%biggest_elem), ncoarsen, errcode, &
           delta_nelem, delta_nelem_leaf, delta_nedge,  delta_nvert, &
           delta_nelem_leaf_own,  delta_nvert_own, delta_dof, delta_dof_own, &
           total_delta_nelem, total_delta_nelem_leaf, total_delta_nedge, &
           total_delta_nvert, total_delta_nelem_leaf_own, &
           total_delta_nvert_own, total_delta_dof, total_delta_dof_own
!----------------------------------------------------
! Begin executable code

any_changes = .false.

! For each level, finest to coarsest ...

do lev=grid%nlev,2,-1

! Make a list of the parents of elements of this level that want to be,
! coarsened, but only include one of an element and its mate.

   call make_h_coarsen_list(grid,leaf_list,nleaf,lev,desired_level, &
                            desired_degree,coarsen_list,ncoarsen)

! Coarsen each element on the list.

      total_delta_dof = 0
      total_delta_dof_own = 0
      total_delta_nelem = 0
      total_delta_nelem_leaf = 0
      total_delta_nelem_leaf_own = 0
      total_delta_nedge = 0
      total_delta_nvert = 0
      total_delta_nvert_own = 0

!$omp parallel do &
!$omp  default(shared) &
!$omp  private(i,errcode,delta_nelem, delta_nelem_leaf, delta_nedge, &
!$omp          delta_nvert, delta_nelem_leaf_own, delta_nvert_own, delta_dof, &
!$omp          delta_dof_own) &
!$omp reduction(+ : total_delta_nelem, total_delta_nelem_leaf, &
!$omp               total_delta_nedge, total_delta_nvert, &
!$omp               total_delta_nelem_leaf_own, total_delta_nvert_own, &
!$omp               total_delta_dof, total_delta_dof_own) &
!$omp  reduction(.or.: any_changes)

   do i=1,ncoarsen
      call h_coarsen_element(grid,coarsen_list(i),errcode,refine_control, &
                                  delta_nelem, delta_nelem_leaf, delta_nedge, &
                                  delta_nvert, delta_nelem_leaf_own, &
                                  delta_nvert_own, delta_dof, delta_dof_own, &
                                  delay_errind=.true.)

      if (errcode == 0) any_changes = .true.
      total_delta_dof = total_delta_dof + delta_dof
      total_delta_dof_own = total_delta_dof_own + delta_dof_own
      total_delta_nelem = total_delta_nelem + delta_nelem
      total_delta_nelem_leaf = total_delta_nelem_leaf + delta_nelem_leaf
      total_delta_nelem_leaf_own = total_delta_nelem_leaf_own + delta_nelem_leaf_own
      total_delta_nedge = total_delta_nedge + delta_nedge
      total_delta_nvert = total_delta_nvert + delta_nvert
      total_delta_nvert_own = total_delta_nvert_own + delta_nvert_own

   end do
!$omp end parallel do

! Do things that must be done OpenMP-sequentially after the parallel loop.

   grid%dof = grid%dof + total_delta_dof
   grid%dof_own = grid%dof_own + total_delta_dof_own
   grid%nelem = grid%nelem + total_delta_nelem
   grid%nelem_leaf = grid%nelem_leaf + total_delta_nelem_leaf
   grid%nelem_leaf_own = grid%nelem_leaf_own + total_delta_nelem_leaf_own
   grid%nedge = grid%nedge + total_delta_nedge
   grid%nvert = grid%nvert + total_delta_nvert
   grid%nvert_own = grid%nvert_own + total_delta_nvert_own

end do ! level

end subroutine h_coarsen_grid

!          -------------------
subroutine make_h_coarsen_list(grid,leaf_list,nleaf,lev,desired_level, &
                               desired_degree,coarsen_list,ncoarsen)
!          -------------------

!----------------------------------------------------
! This routine makes a list of the parents of elements of level lev that
! want to be coarsened, but only includes one of an element and its mate.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(inout) :: leaf_list(:), nleaf
integer, intent(in) :: lev
integer, intent(inout) :: desired_level(:), desired_degree(:)
integer, intent(out) :: coarsen_list(:), ncoarsen
!----------------------------------------------------
! Local variables:

integer :: i, elem, parent, mate, new_nleaf, sibling
! TEMP110504 this should be small_logical; search for others
logical :: onlist(grid%biggest_elem)
!----------------------------------------------------
! Begin executable code

ncoarsen = 0

! Mark all elements as not on the list

onlist = .false.

! For each element on the list that has level lev ...

new_nleaf = nleaf
do i=1,nleaf
   elem = leaf_list(i)
   if (grid%element(elem)%level /= lev) cycle

! See if it wants to be coarsened

   if (grid%element(elem)%level <= desired_level(elem)) cycle

! Identify the parent and mate

   parent = hash_decode_key(grid%element(elem)%gid/MAX_CHILD,grid%elem_hash)
   if (.not. grid%element(parent)%mate == BOUNDARY) then
      mate = hash_decode_key(grid%element(parent)%mate,grid%elem_hash)
   else
      mate = BOUNDARY
   endif

! If either the parent or the mate is already on the list, skip.

   if (onlist(parent)) cycle
   if (mate /= BOUNDARY) then
      if (onlist(mate)) cycle
   endif

! Add the parent to the list.

   ncoarsen = ncoarsen + 1
   coarsen_list(ncoarsen) = parent
   onlist(parent) = .true.

! Add the parent and mate to the leaf list, and assign desired level and degree

   leaf_list(new_nleaf+1) = parent
   desired_level(parent) = desired_level(elem)
   desired_degree(parent) = desired_degree(elem)
   if (mate == BOUNDARY) then
      new_nleaf = new_nleaf + 1
   else
      new_nleaf = new_nleaf + 2
      leaf_list(new_nleaf) = mate
      desired_level(mate) = desired_level(elem)
      sibling = get_child_lid(grid%element(parent)%mate,1,grid%elem_hash)
      desired_degree(mate) = desired_degree(sibling)
   endif

end do

nleaf = new_nleaf

end subroutine make_h_coarsen_list

!          -------------
subroutine p_refine_grid(grid,refine_control,desired_degree,elem_list,nelem, &
                         any_changes)
!          -------------

!----------------------------------------------------
! This routine p-refines elements in elem_list for which the current degree is
! less than the desired degree.  The edge rule does not have to be satisfied
! on entry, but it will be satisfied on return.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: desired_degree(:), elem_list(:), nelem
logical, intent(out) :: any_changes
!----------------------------------------------------
! Local variables:

! TEMP seems like edge_list need not be bigger than grid%nedge.  Verify and
!      look for similar places where biggest is used for dimension
integer :: total_delta_dof, total_delta_dof_own, i, errcode, delta_dof, &
           delta_dof_own, edge_list(grid%biggest_edge), nedge
!----------------------------------------------------
! Begin executable code

any_changes = .false.

! p-refine the interior of each element that wants to be p refined.

total_delta_dof = 0
total_delta_dof_own = 0

!$omp parallel &
!$omp  default(shared) &
!$omp  private(i,errcode,delta_dof,delta_dof_own) &
!$omp  reduction(+ : total_delta_dof, total_delta_dof_own) &
!$omp  reduction(.or.: any_changes)

!$omp do schedule(dynamic)
do i=1,nelem
   do while (grid%element(elem_list(i))%degree < desired_degree(elem_list(i)))
      call p_refine_element_interior(grid,refine_control,elem_list(i), &
                                     errcode,delta_dof,delta_dof_own)
      if (errcode == 0) any_changes = .true.
      total_delta_dof = total_delta_dof + delta_dof
      total_delta_dof_own = total_delta_dof_own + delta_dof_own
   end do
end do
!$omp end do

! do things that must be done OpenMP-sequentially after the parallel loop
!$omp single

grid%dof = grid%dof + total_delta_dof
grid%dof_own = grid%dof_own + total_delta_dof_own

! make a list of all the edges that no longer satisfy the edge rule

call list_edges_without_rule(grid,refine_control,elem_list,nelem,edge_list, &
                             nedge)

! enforce the edge rule on the listed edges

total_delta_dof = 0
total_delta_dof_own = 0

!$omp end single
!$omp do
do i=1,nedge
   call enforce_edge_rule(grid,refine_control,edge_list(i),delta_dof, &
                          delta_dof_own)
   total_delta_dof = total_delta_dof + delta_dof
   total_delta_dof_own = total_delta_dof_own + delta_dof_own
end do
!$omp end do

!$omp end parallel

! do things that must be done OpenMP-sequentially after the parallel loop

grid%dof = grid%dof + total_delta_dof
grid%dof_own = grid%dof_own + total_delta_dof_own

end subroutine p_refine_grid

!          -------------
subroutine h_refine_grid(grid,refine_control,solver_control,desired_level, &
                         desired_degree,elem_list,nelem,any_changes)
!          -------------

!----------------------------------------------------
! This routine h-refines elements in elem_list for which the current level is
! less than the desired level.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
integer, pointer :: desired_level(:), desired_degree(:), elem_list(:)
integer, intent(in) :: nelem
logical, intent(out) :: any_changes
!----------------------------------------------------
! Local variables:

integer :: lev, elem, errcode, mate, parent, i, nrefine, newsize, &
           one_vert_lid(2),one_edge_lid(6),one_elem_lid(4), &
           delta_dof, total_delta_dof, delta_dof_own, total_delta_dof_own, &
           delta_nelem, total_delta_nelem, delta_nelem_leaf, &
           total_delta_nelem_leaf, delta_nelem_leaf_own, &
           total_delta_nelem_leaf_own, delta_nedge, total_delta_nedge, &
           delta_nedge_own, total_delta_nedge_own, delta_nvert, &
           total_delta_nvert, delta_nvert_own, total_delta_nvert_own, &
           max_nlev, max_max_nlev, astat, refedge, nrefedge
integer, allocatable :: refine_list(:), vert_lid(:,:), edge_lid(:,:), &
                        elem_lid(:,:), refinement_edge_list(:)
logical, allocatable :: onlist(:)
logical :: reallocated

!----------------------------------------------------
! Begin executable code

any_changes = .false.

! Allocate refine_list and lid lists

allocate(refine_list(grid%biggest_elem), vert_lid(2,grid%biggest_elem), &
         edge_lid(6,grid%biggest_elem), elem_lid(4,grid%biggest_elem), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in h_refine_grid")
   stop
endif

! Allocate the list of refinement edges

if (global_element_kind == TETRAHEDRAL_ELEMENT) then
   allocate(refinement_edge_list(grid%biggest_edge),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in h_refine_grid")
      stop
   endif
   nrefedge = 0
endif

! Allocate onlist

allocate(onlist(size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in h_refine_grid")
   stop
endif

! Note: grid%dof etc must not be used in before_h_refine or after_h_refine
! because the reduction does not occur until the end of the parallel section

total_delta_dof = 0
total_delta_dof_own = 0
total_delta_nelem = 0
total_delta_nelem_leaf = 0
total_delta_nelem_leaf_own = 0
total_delta_nedge = 0
total_delta_nedge_own = 0
total_delta_nvert = 0
total_delta_nvert_own = 0
max_max_nlev = 0

!$omp parallel &
!$omp  default(shared) &
!$omp  private(lev,i,errcode,one_elem_lid,one_edge_lid,one_vert_lid, &
!$omp    delta_dof, delta_dof_own, delta_nelem, delta_nelem_leaf, &
!$omp    delta_nelem_leaf_own, delta_nedge, delta_nedge_own, &
!$omp    delta_nvert, delta_nvert_own, max_nlev, refedge) &
!$omp  reduction(+ : total_delta_dof, total_delta_dof_own, total_delta_nelem, &
!$omp    total_delta_nelem_leaf, total_delta_nelem_leaf_own, total_delta_nedge,&
!$omp    total_delta_nedge_own, total_delta_nvert, total_delta_nvert_own) &
!$omp  reduction(max : max_max_nlev) &
!$omp  reduction(.or.: any_changes)

! for each level, starting with the coarsest level ...

! TEMP3D I think they can be done all at once rather than by level

do lev=1,grid%nlev

!$omp barrier
!$omp single

! Make sure onlist is still big enough.

   if (size(onlist) /= size(grid%element)) then
      deallocate(onlist)
      allocate(onlist(size(grid%element)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in h_refine_grid")
         stop
      endif
   endif

! Create a list of elements on this level to be refined.
! There is a slight possibility that refine_list and the lid lists will
! come up short.  If that happens, double the allocation and start over.

   reallocated = .true.
   do while (reallocated)
      reallocated = .false.
      onlist(1:grid%biggest_elem) = .false.
      nrefine = 0
      do i=1,nelem
         elem = elem_list(i)
         if (grid%element(elem)%level == lev .and. &
             grid%element(elem)%level < desired_level(elem)) then
            if (onlist(elem)) cycle
            if (.not. grid%element(elem)%mate == BOUNDARY) then
               mate = hash_decode_key(grid%element(elem)%mate,grid%elem_hash)
            else
               mate = BOUNDARY
            endif
            if (mate /= BOUNDARY) then
               if (onlist(mate)) cycle
            endif
            nrefine = nrefine + 1
            if (nrefine > size(refine_list)) then
               newsize = 2*size(refine_list)
               deallocate(refine_list,vert_lid,edge_lid,elem_lid)
               allocate(refine_list(newsize), vert_lid(2,newsize), &
                        edge_lid(6,newsize), elem_lid(4,newsize), stat=astat)
               if (astat /= 0) then
                  ierr = ALLOC_FAILED
                  call fatal("memory allocation failed in h_refine_grid")
                  stop
               endif
               reallocated = .true.
               exit
            endif
            refine_list(nrefine) = elem
            onlist(elem) = .true.
         endif
      end do
   end do

! Get lids for the children, and other things that must be done
! OpenMP-sequentially before the parallel loop.

   if (nrefine /= 0) then
      call before_h_refine(grid,refine_list,nrefine,vert_lid=vert_lid, &
                           edge_lid=edge_lid,elem_lid=elem_lid, &
                           desired_level=desired_level, &
                           desired_degree=desired_degree, &
                           elem_list=elem_list)
   endif

!$omp end single

! If the list is empty, move on to the next level.

   if (nrefine == 0) cycle

!$omp do
   do i=1,nrefine
! TEMP I should be able to pass these array sections, but ifort gets SIGSEGV
      one_elem_lid = elem_lid(:,i)
      one_edge_lid = edge_lid(:,i)
      one_vert_lid = vert_lid(:,i)
      call h_refine_element(grid,refine_list(i),errcode,refedge,refine_control,&
                   solver_control,vert_child_lid=one_vert_lid, &
                   edge_child_lid=one_edge_lid,elem_child_lid=one_elem_lid, &
                   delta_dof=delta_dof,delta_dof_own=delta_dof_own, &
                   delta_nelem=delta_nelem,delta_nelem_leaf=delta_nelem_leaf, &
                   delta_nelem_leaf_own=delta_nelem_leaf_own, &
                   delta_nedge=delta_nedge,delta_nedge_own=delta_nedge_own, &
                   delta_nvert=delta_nvert,delta_nvert_own=delta_nvert_own, &
                   max_nlev=max_nlev,delay_errind=.true., &
                   desired_level=desired_level,desired_degree=desired_degree, &
                   elem_list=elem_list)

      if (errcode == 0) any_changes = .true.
      total_delta_dof = total_delta_dof + delta_dof
      total_delta_dof_own = total_delta_dof_own + delta_dof_own
      total_delta_nelem = total_delta_nelem + delta_nelem
      total_delta_nelem_leaf = total_delta_nelem_leaf + delta_nelem_leaf
      total_delta_nelem_leaf_own = total_delta_nelem_leaf_own + delta_nelem_leaf_own
      total_delta_nedge = total_delta_nedge + delta_nedge
      total_delta_nedge_own = total_delta_nedge_own + delta_nedge_own
      total_delta_nvert = total_delta_nvert + delta_nvert
      total_delta_nvert_own = total_delta_nvert_own + delta_nvert_own
      max_max_nlev = max(max_max_nlev,max_nlev)

! TEMP3D OPENMP this will not work in parallel
      if (global_element_kind == TETRAHEDRAL_ELEMENT) then
         if (refedge /= -1) then
            nrefedge = nrefedge + 1
            refinement_edge_list(nrefedge) = refedge
         endif
      endif

   end do
!$omp end do

! Do things that must be done OpenMP-sequentially after the parallel loop.

!$omp single

   call after_h_refine(grid,refine_list,nrefine,vert_lid,edge_lid,elem_lid)

!$omp end single

end do

!$omp end parallel

! total_* are reduced at end parallel, now add them to grid

grid%dof = grid%dof + total_delta_dof
grid%dof_own = grid%dof_own + total_delta_dof_own
grid%nelem = grid%nelem + total_delta_nelem
grid%nelem_leaf = grid%nelem_leaf + total_delta_nelem_leaf
grid%nelem_leaf_own = grid%nelem_leaf_own + total_delta_nelem_leaf_own
grid%nedge = grid%nedge + total_delta_nedge
grid%nedge_own = grid%nedge_own + total_delta_nedge_own
grid%nvert = grid%nvert + total_delta_nvert
grid%nvert_own = grid%nvert_own + total_delta_nvert_own
grid%nlev = max(grid%nlev,max_max_nlev)

! Free memory.

deallocate(onlist,refine_list,vert_lid,edge_lid,elem_lid)

! for tetrahedra, make sure there are no hanging nodes

if (global_element_kind == TETRAHEDRAL_ELEMENT) then
   call remove_hanging_nodes(grid,refinement_edge_list,nrefedge, &
                             refine_control,solver_control, &
                             desired_level,desired_degree,elem_list)
   deallocate(refinement_edge_list)
endif

end subroutine h_refine_grid

!          ---------------
subroutine enforce_overlap(grid,refine_control,solver_control,procs)
!          ---------------

!----------------------------------------------------
! This routine enforce the overlap conditions around the partition boundary.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
type(proc_info), intent(in) :: procs

!----------------------------------------------------
! Local variables:

integer :: part_bound_elements(grid%biggest_elem), npartbound, ni, nr, i, j, &
           k, count, p, nproc, my_processor, lid, err, lev, elem, astat
type(hash_key) :: gid, ancestor
integer, allocatable :: isend(:), nrecv(:)
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
logical :: on_boundary
logical, allocatable :: isboundaryvert(:)
!----------------------------------------------------
! Begin executable code

! Convenience variables.

my_processor = my_proc(procs)
nproc = num_proc(procs)

! Mark each vertex as being on the partition boundary or not by going through
! the elements and marking each vertex as being on the boundary if I own the
! element or vertex but not the other.

allocate(isboundaryvert(size(grid%vertex)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in enforce_overlap")
   stop
endif
isboundaryvert = .false.

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do i=1,VERTICES_PER_ELEMENT
         if (grid%element(elem)%iown .neqv. &
             grid%element(grid%vertex(grid%element(elem)%vertex(i))%assoc_elem)%iown) then
            isboundaryvert(grid%element(elem)%vertex(i)) = .true.
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

! Make a list of owned leaf elements that touch the partition boundary.

call list_elements(grid,part_bound_elements,npartbound,own=.true., &
                   leaf=.true.,bound_vert=.true.)

! Create a message containing the GIDs and degrees of those elements.

ni = npartbound*((VERTICES_PER_ELEMENT+1)*KEY_SIZE+1)
allocate(isend(ni),nrecv(nproc),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in enforce_overlap")
   stop
endif
count = 1
do i=1,npartbound
   do j=1,VERTICES_PER_ELEMENT
      call hash_pack_key(grid%vertex(grid%element(part_bound_elements(i))%vertex(j))%gid,isend,count)
      count = count + KEY_SIZE
   end do
   call hash_pack_key(grid%element(part_bound_elements(i))%gid,isend,count)
   count = count + KEY_SIZE
   isend(count) = grid%element(part_bound_elements(i))%degree
   count = count + 1
end do

! Exchange the message with the other processors.

if (my_processor /= MASTER) then
   call phaml_alltoall(procs,isend,ni,irecv,nrecv,3340)
endif

! Receive and process the messages.

if (my_processor /= MASTER) then

! For each element in the received message ...

   count = 1
   do j=1,sum(nrecv)/((VERTICES_PER_ELEMENT+1)*KEY_SIZE+1)

! See if the element is on my partition boundary (is one of the vertices on
! the boundary).

      on_boundary = .false.
      do k=1,VERTICES_PER_ELEMENT
         gid = hash_unpack_key(irecv,count)
         count = count + KEY_SIZE
         lid = hash_decode_key(gid,grid%vert_hash)
         if (lid /= HASH_NOT_FOUND) then
            if (isboundaryvert(lid)) then
               on_boundary = .true.
            endif
         endif
      end do

      gid = hash_unpack_key(irecv,count) ! sent element gid
      count = count + KEY_SIZE

! It is, enforce overlap.

      if (on_boundary) then

! Identify the first ancestor of the element that I have.

         ancestor = gid
         lid = hash_decode_key(gid,grid%elem_hash)
         do while (lid == HASH_NOT_FOUND)
            ancestor = ancestor/MAX_CHILD
            lid = hash_decode_key(ancestor,grid%elem_hash)
         end do

! If ancestor is not the element, then create the element.

         if (.not. ancestor == gid) then

            call create_element(grid,gid,refine_control,solver_control,err)
! TEMP120206 are there arrays I need to check for increasing size?
!            It appears isboundaryvert should be checked
            if (err /= 0) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("enforce_overlap: create_element returned error", &
                          intlist=(/err/))
               stop
            endif
         endif

! Set p for as sent.

         lid = hash_decode_key(gid,grid%elem_hash)
         if (lid == HASH_NOT_FOUND) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("enforce_overlap: lid is HASH_NOT_FOUND")
            stop
         endif
         call enforce_p(lid)

      endif

      count = count + 1
   end do ! elements in message

   deallocate(irecv)
endif ! not MASTER

deallocate(isboundaryvert,isend,nrecv)

contains
!-------

recursive subroutine enforce_p(pelem)
!                    ---------
integer :: pelem, child(MAX_CHILD)

! If elem is not a leaf, do the children.

child = get_child_lid(grid%element(pelem)%gid,ALL_CHILDREN,grid%elem_hash)
if (child(1) /= NO_CHILD) then
   call enforce_p(child(1))
   call enforce_p(child(2))

else

! Perform p-coarsening/p-refinement until the degree is right

   do while (grid%element(pelem)%degree > irecv(count))
      call p_coarsen_elem(grid,pelem,err,refine_control)
   end do
   do while (grid%element(pelem)%degree < irecv(count))
      call p_refine_elem(grid,pelem,refine_control,errcode=err)
      if (err /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call warning("enforce_overlap: p refinement failed", &
             intlist=(/err,pelem,grid%element(pelem)%degree,irecv(count)/))
         exit
      endif
   end do

endif

end subroutine enforce_p

end subroutine enforce_overlap

!          -----------------
subroutine pred_load_balance(grid,procs,refine_control,solver_control,lb, &
                             partition_method,still_sequential,balance_what, &
                             new_level,new_p,any_changes)
!          -----------------

!----------------------------------------------------
! This routine performs predictive load balancing based on new_level and new_p.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info) :: procs
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
type(Zoltan_Struct), pointer :: lb
integer, intent(in) :: partition_method
logical, intent(in) :: still_sequential
integer, intent(in) :: balance_what, new_level(:), new_p(:)
logical, intent(out) :: any_changes
!----------------------------------------------------
! Local variables:

type(hash_key), pointer :: export_gid(:)
integer, pointer :: export_proc(:)
logical :: repart, redist
integer :: nentity, minentity, maxentity, numexp, astat, p, nirecv, nrrecv
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)

!----------------------------------------------------
! Begin executable code

any_changes = .false.
nullify(export_gid,export_proc)

! see if repartitioning is really necessary

! some cases where partitioning is not needed:
!   master does not participate
!   number of processors is 1
!   program compiled for sequential execution
!   balancing nothing

if (my_proc(procs) == MASTER .or. num_proc(procs) == 1 .or. &
    PARALLEL==SEQUENTIAL .or. balance_what==BALANCE_NONE) then
   repart = .false.

! otherwise, partition if the load is sufficiently out of balance

else
   select case (balance_what)
   case (BALANCE_ELEMENTS)
      call get_grid_info(grid,procs,still_sequential,3360, &
                         nelem_leaf_own=nentity,no_master=.true.)
   case (BALANCE_VERTICES)
      call get_grid_info(grid,procs,still_sequential,3360, &
                         nvert_own=nentity,no_master=.true.)
   case (BALANCE_EQUATIONS)
      call get_grid_info(grid,procs,still_sequential,3360, &
                         dof_own=nentity,no_master=.true.)
   case default
      call fatal("bad value for balance what in subroutine pred_load_balance")
   end select
   minentity = phaml_global_min(procs,nentity,3370)
   maxentity = phaml_global_max(procs,nentity,3380)
   repart = ( (maxentity-minentity)/float(maxentity) > 0.05)
endif

! partition

if (repart) then
   call lightweight_set_weights(grid,balance_what,new_level,new_p)
   call partition(grid,procs,lb,refine_control,.true.,partition_method, &
                  balance_what,num_proc(procs),export_gid,export_proc, &
                  weights_already_set=.true.)

! redistribute

! see if enough elements are moved to be worthwhile
   if (associated(export_gid)) then
      numexp = size(export_gid)
   else
      numexp = 0
   endif
   numexp = phaml_global_sum(procs,numexp,3390)
   call get_grid_info(grid,procs,still_sequential,3395, &
                      total_nelem_leaf=nentity,no_master=.true.)
   redist = (numexp/float(nentity) > .05)
   if (redist) then
      call redistribute(grid,procs,refine_control,solver_control,export_gid, &
                        export_proc)
      call reconcile(grid,procs,refine_control,solver_control,still_sequential)
      any_changes = .true.
   endif
   if (associated(export_gid)) then
      deallocate(export_gid,export_proc)
   endif

endif

! inform the master about any changes

if (my_proc(procs) == 1) then
   if (any_changes) then
      call phaml_send(procs,MASTER,(/1/),1,(/0.0_my_real/),0,3391)
   else
      call phaml_send(procs,MASTER,(/0/),1,(/0.0_my_real/),0,3391)
   endif
elseif (my_proc(procs) == MASTER) then
   call phaml_recv(procs,p,irecv,nirecv,rrecv,nrrecv,3391)
   any_changes = irecv(1) == 1
   deallocate(irecv)
endif

end subroutine pred_load_balance

!========================================================================

!          ---------------
subroutine old_refine_adaptive(grid,procs,refine_control,solver_control,io_control,&
                           still_sequential,init_nvert,init_nelem,init_dof, &
                           loop,balance_what,predictive,no_time)
!          ---------------

!----------------------------------------------------
! This routine refines the grid in accordance with the parameters in
! refine_control.
! If derefine is true, the grid is first derefined to remove elements
! with very small error indicators.
! If refterm is DOUBLE_NVERT (NELEM, NEQ) then the number of vertices
! (elements, degrees of freedom) in the global grid is increased to
! init_nvert*inc_factor**loop. (init_nelem, init_dof)
! Each processor determines how many vertices (elements, dofs) by which to
! increase the grid based on the error indicators relative to those of other
! processors, unless load balancing is predictive in which case all processors
! get an equal share of the vertices (elements, dofs).
! If refterm is HALVE_ERREST then the maximum error indicator is reduced
! by the factor inc_factor.  However, if this leads to a very small change
! in the grid, additional refinement is performed.
! If refterm is KEEP_NVERT (NELEM, ERREST) then the number of vertices
! (elements) is kept constant, or the maximum error indicator is kept
! approximately constant.
! If refterm is ONE_REF then elements whose error indicator is larger than
! reftol are refined, and each element is refined at most once.
! If refterm is ONE_REF_HALF_ERRIND then elements whose error indicator is
! larger than maximum_error_indicator / inc_factor are refined once.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), target, intent(in) :: procs
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
type(io_options), intent(in) :: io_control
integer, intent(in) :: init_nvert, init_nelem, init_dof, loop, balance_what
logical, intent(in) :: still_sequential, predictive
logical, intent(in), optional :: no_time

!----------------------------------------------------
! Local variables:

integer :: elem, errcode, astat, nproc, my_processor, target, mate, &
           child(MAX_CHILD), refedge
real(my_real) :: global_max_errind
integer, pointer :: numhref(:), numpref(:), new_p(:,:)
! If len=1 is changed, also change it for temp_reftype in more_elements
character(len=1), pointer :: reftype(:)
type(errind_list) :: elist
logical :: return_to_elist, complete_elist, one_elist, target_met

!----------------------------------------------------
! Begin executable code

! Not modified for 3D

if (global_element_kind == TETRAHEDRAL_ELEMENT) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("old_refine_adaptive has not been written for tetrahedra")
   stop
endif

! MASTER doesn't participate

if (my_proc(procs) == MASTER) return

nproc = num_proc(procs)
my_processor = my_proc(procs)

! compute error indicators if needed

if (.not. grid%errind_up2date) then
   call all_error_indicators(grid,refine_control%error_estimator)
endif

! compute maximum error indicator

global_max_errind = compute_global_max_errind(grid,procs,still_sequential)

! set the target for terminating refinement

call set_target(target,grid,procs,refine_control,global_max_errind, &
                init_nvert,init_nelem,init_dof,loop,balance_what, &
                still_sequential,predictive)

! derefine the grid, if requested

if (refine_control%derefine) then
   call derefinement(grid,procs,refine_control,solver_control, &
                     global_max_errind,target)
endif

! set the number of times to refine each element

allocate(numhref(size(grid%element)),numpref(size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in refine",procs=procs)
   return
endif

call set_numref(grid,procs,still_sequential,refine_control,numhref,numpref)

! create the lists that group elements into bins based on the error indicator

allocate(elist%next_errind(size(grid%element)), &
         elist%prev_errind(size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in refine",procs=procs)
   return
endif

call make_elist(grid,procs,elist,refine_control,global_max_errind, &
                numhref,numpref,predictive,still_sequential)
elist%current_list = 1

! mark elements to be refined by h, by p or not, or leave undecided

allocate(reftype(size(grid%element)),new_p(2,size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in refine",procs=procs)
   return
endif
new_p = 0

call mark_reftype(grid,refine_control,global_max_errind,reftype,elist,new_p)

! set controls for the refine loop:
! return_to_elist, tells if children get put into an elist
! complete_elist, tells if you must finish the whole elist
! one_elist, tells if you quit after finishing the first elist

call set_loop_control(refine_control,return_to_elist,complete_elist,one_elist)

! main refinement loop

errcode = 0
do

! quit if the grid is full

   if (errcode /= 0) exit

! see if the target has been met

   target_met = check_target(grid,refine_control,elist,target)

! if the target is met and we don't complete lists, done

   if ((.not. complete_elist) .and. target_met) exit

! get the next element to refine

   elem = elist%head_errind(elist%current_list)

! if this is the end of the current list and either we do just the first list or
! the target has been met, done

   if (elem == END_OF_LIST .and. (one_elist .or. target_met)) exit

! if this is the end of the current list, find the beginning of the next
! nonempty list

   do while (elem==END_OF_LIST .and. elist%current_list<size(elist%head_errind))
      elist%current_list = elist%current_list + 1
      elem = elist%head_errind(elist%current_list)
   end do

! if all lists are empty, done

   if (elem == END_OF_LIST) exit

! if we return elements to the elists, complete the current elist, and have
! reached the last list, exit to avoid an infinite loop

   if (return_to_elist .and. complete_elist .and. &
       elist%current_list == size(elist%head_errind)) exit

! if we haven't determined the refinement type, do it now

   if (reftype(elem) == "u") then
      call mark_reftype_one(elem,grid,refine_control,global_max_errind, &
                            reftype,.false.,new_p(:,elem))
   endif

! refine the element and determine the type of refinement of the new element(s)

   if (reftype(elem) == "h") then
      call h_refine_element(grid,elem,errcode,refedge,refine_control, &
                            solver_control,elist,reftype,new_p,numhref, &
                            numpref,return_to_elist)
      reftype(elem) = "n"
      child = get_child_lid(grid%element(elem)%gid,ALL_CHILDREN,grid%elem_hash)
      call mark_reftype_one(child(1),grid,refine_control,1.0_my_real,reftype, &
                            .true.,new_p(:,child(1)))
      call mark_reftype_one(child(2),grid,refine_control,1.0_my_real,reftype, &
                            .true.,new_p(:,child(2)))
      if (grid%element(elem)%mate == BOUNDARY) then
         mate = BOUNDARY
      else
         mate = hash_decode_key(grid%element(elem)%mate,grid%elem_hash)
      endif
      if (mate /= BOUNDARY) then
         reftype(mate) = "n"
         child = get_child_lid(grid%element(mate)%gid,ALL_CHILDREN, &
                               grid%elem_hash)
         call mark_reftype_one(child(1),grid,refine_control,1.0_my_real, &
                               reftype,.true.,new_p(:,child(1)))
         call mark_reftype_one(child(2),grid,refine_control,1.0_my_real, &
                               reftype,.true.,new_p(:,child(2)))
      endif
   elseif (reftype(elem) == "p") then
      call p_refine_elem(grid,elem,refine_control,elist,numpref,return_to_elist)
      call mark_reftype_one(elem,grid,refine_control,1.0_my_real,reftype, &
                            .true.,new_p(:,elem))
   else
      call remove_from_errind_list(elem,elist)
   endif

end do ! main refine loop

! free memory

deallocate(elist%next_errind,elist%prev_errind,reftype,new_p,numhref,numpref, &
           stat=astat)

end subroutine old_refine_adaptive

!          ----------
subroutine set_target(target,grid,procs,refine_control,global_max_errind, &
                      init_nvert,init_nelem,init_dof,loop,balance_what, &
                      still_sequential,predictive)
!          ----------

!----------------------------------------------------
! This routine determines the value to use as the target for terminating
! refinement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(out) :: target
type(grid_type), intent(inout) :: grid
type (proc_info), target, intent(in) :: procs
type(refine_options), intent(in) :: refine_control
real(my_real), intent(in) :: global_max_errind
integer, intent(in) :: init_nvert, init_nelem, init_dof, loop, balance_what
logical, intent(in) :: still_sequential, predictive
!----------------------------------------------------
! Local variables:

real(my_real) :: my_fraction
integer :: my_num_big, lev, elem, total_num_big, total
!----------------------------------------------------
! Begin executable code

! cases where we don't need target

if (refine_control%refterm == ONE_REF .or. &
    refine_control%refterm == ONE_REF_HALF_ERRIND) then
   target = 0
   return
endif

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_T3S) then
   if (refine_control%t3s_reftype == H_UNIFORM) then
      target = 0
      return
   endif
endif

! cases that don't need my_fraction

if (refine_control%refterm == HALVE_ERREST) then
   target = max(1,nint(log(refine_control%inc_factor)/log(binw)))
   return
endif

if (refine_control%refterm == KEEP_ERREST) then
! TEMP not doing KEEP_ERREST
   call warning("have not yet decided how to handle refterm==KEEP_ERREST")
   target = 0
endif

! set the target for the termination criterion

! if load balancing is not predictive, determine my fraction by comparing the
! number of elements with large error indicators on this processor with those
! on other processors.  For now, large means it will be in the first two error
! indicator lists

if (still_sequential) then
   my_fraction = 1.0_my_real
elseif (predictive) then
   my_fraction = my_weight_fraction(grid,procs,predictive, &
                                    balance_what,refine_control)
else
   my_num_big = 0
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf) then
            if (maxval(grid%element_errind(elem,:))/grid%element(elem)%work > &
                global_max_errind/(binw**2)) then
               my_num_big = my_num_big + 1
            endif
         endif
         elem = grid%element(elem)%next
      end do
   end do
   total_num_big = phaml_global_sum(procs,my_num_big,3322)
   my_fraction = my_num_big/real(total_num_big,my_real)
endif

! set the target value depending on what type of termination is used

select case (refine_control%refterm)

case (DOUBLE_NVERT)

   total = init_nvert*refine_control%inc_factor**loop
   if (total > refine_control%max_vert) total = refine_control%max_vert
   target = ceiling(total*my_fraction)

case (DOUBLE_NELEM)

   total = init_nelem*refine_control%inc_factor**loop
   if (total > refine_control%max_elem) total = refine_control%max_elem
   target = ceiling(total*my_fraction)

case (DOUBLE_NEQ)

   total = init_dof*refine_control%inc_factor**loop
   if (total > refine_control%max_dof) total = refine_control%max_dof
   target = ceiling(total*my_fraction)

case (KEEP_NVERT)

   if (refine_control%max_vert /= huge(0)) then
      total = refine_control%max_vert
   else
      call get_grid_info(grid,procs,still_sequential,3345,total_nvert=total,&
                         no_master=.true.)
      if (total > refine_control%max_vert) total = refine_control%max_vert
   endif
   target = total*my_fraction

case (KEEP_NELEM)

   if (refine_control%max_elem /= huge(0)) then
      total = refine_control%max_elem
   else
      call get_grid_info(grid,procs,still_sequential,3345, &
                         total_nelem_leaf=total,no_master=.true.)
      if (total > refine_control%max_elem) total = refine_control%max_elem
   endif
   target = total*my_fraction

case (KEEP_NEQ)

   if (refine_control%max_dof /= huge(0)) then
      total = refine_control%max_dof
   else
      call get_grid_info(grid,procs,still_sequential,3345, &
                         total_dof=total,no_master=.true.)
      if (total > refine_control%max_dof) total = refine_control%max_dof
   endif
   target = total*my_fraction

case default

   call fatal("illegal value for refterm",procs=procs)
   stop

end select

end subroutine set_target

!        ------------------
function my_weight_fraction(grid,procs,predictive,balance_what,refine_control)
!        ------------------

!----------------------------------------------------
! This routine computes the fraction of the total weight of leaf elements that
! belongs to this processor.  For homogeneous load balancing, it should be
! 1/nproc.  This is really only useful if load balancing is done for a
! heterogeneous computing system.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: predictive
integer, intent(in) :: balance_what
type(refine_options), intent(in) :: refine_control
real(my_real) :: my_weight_fraction
!----------------------------------------------------
! Local variables:

real(my_real) :: my_total_weight, total_weight
integer :: lev, elem
!----------------------------------------------------
! Begin executable code

! set the current weights

call set_weights(grid,predictive,balance_what,refine_control,procs)

my_total_weight = 0.0_my_real
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf .and. grid%element(elem)%iown) then
         my_total_weight = my_total_weight + grid%element(elem)%weight
      endif
      elem = grid%element(elem)%next
   end do
end do

total_weight = phaml_global_sum(procs,my_total_weight,3321)

if (total_weight /= 0.0_my_real) then
   my_weight_fraction = my_total_weight/total_weight
else
   my_weight_fraction = 1.0_my_real
endif

end function my_weight_fraction

!          ------------
subroutine derefinement(grid,procs,refine_control,solver_control,max_errind, &
                        target)
!          ------------

!----------------------------------------------------
! This routine performs h and p derefinement of the grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
real(my_real), intent(in) :: max_errind
integer, intent(in) :: target
!----------------------------------------------------
! Local variables:

integer :: deref_list(2*size(grid%element)), degree_list(2*size(grid%element))
integer :: num_deref, lev, elem, parent, siblings(4), i, errcode, astat, &
           nproc, my_processor, p, next_elem, isub
! newcomm
integer :: nsend
integer, allocatable :: isend(:), nrecv(:)
integer, pointer :: irecv(:)
logical(small_logical) :: hcoarsen_checked(size(grid%element))
logical :: hcoarsen_possible, pcoarsen_possible, hcoarsen_occured, &
           pcoarsen_occured
real(my_real) :: hcoarsen_errind, pcoarsen_errind, small_errind
!----------------------------------------------------
! Begin executable code

! no derefinement for uniform refinement

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_T3S) then
   if (refine_control%t3s_reftype == H_UNIFORM) then
      return
   endif
endif

! unrefine all unrefineable elements whose h or p coarsening error indicator
! is small enough.  If that's not enough, increase the size allowed for the
! coarsening indicator and repeat.

num_deref = 0
small_errind = max_errind/100

do

   hcoarsen_checked = .false.

! go through the elements from finest level to coarsest level to allow
! parents of unrefined elements to also be unrefined

   do lev=grid%nlev,1,-1
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)

! relatives

         if (grid%element(elem)%level /= 1) then
            parent = hash_decode_key(grid%element(elem)%gid/MAX_CHILD, &
                                     grid%elem_hash)
            siblings(1:2) = get_child_lid(grid%element(parent)%gid, &
                                          ALL_CHILDREN,grid%elem_hash)
            if (grid%element(parent)%mate == BOUNDARY) then
               siblings(3:4) = NO_CHILD
            else
               siblings(3:4) = get_child_lid( &
                                 grid%element(parent)%mate,ALL_CHILDREN, &
                                 grid%elem_hash)
            endif
         endif

! Check conditions that allow h coarsening.
! If the element has already been checked for h coarsening (via a sibling)
! then don't need to check it again.
! Cannot h coarsen level 1 elements, elements that have children, elements
! whose siblings have children, or elements whose parent or parent's mate have
! children owned by different processors.  Only h coarsen elements that I own
! or I own a sibling.  Only do h coarsening with h or hp adaptive refinement.

         hcoarsen_possible = .not. hcoarsen_checked(elem)
         if (refine_control%reftype /= H_ADAPTIVE .and. &
             refine_control%reftype /= HP_ADAPTIVE) hcoarsen_possible = .false.
         if (grid%element(elem)%level == 1) hcoarsen_possible = .false.
         if (.not. grid%element(elem)%isleaf) hcoarsen_possible=.false.
         if (hcoarsen_possible) then
            do i=1,4
               if (i==3 .and. grid%element(parent)%mate == BOUNDARY) exit
               if (.not. grid%element(siblings(i))%isleaf) then
                  hcoarsen_possible = .false.
               endif
            end do
            if (grid%element(siblings(1))%iown .neqv. &
                grid%element(siblings(2))%iown) then
               hcoarsen_possible = .false.
            endif
            if (.not. grid%element(parent)%mate == BOUNDARY) then
               if (grid%element(siblings(3))%iown .neqv. &
                   grid%element(siblings(4))%iown) then
                  hcoarsen_possible = .false.
               endif
            endif
            if (grid%element(parent)%mate == BOUNDARY) then
               if (.not. grid%element(elem)%iown) hcoarsen_possible = .false.
            else
               if (.not. grid%element(elem)%iown .and. &
                .not. grid%element(siblings(3))%iown) hcoarsen_possible=.false.
            endif
         endif

! Check conditions that allow p coarsening.
! Cannot p coarsen an element that has p==1.  Only leaves need be p coarsened.
! Only p coarsen elements that I own.  Only do p coarsening with p or hp
! adaptive.

         pcoarsen_possible = grid%element(elem)%isleaf .and. &
                             grid%element(elem)%iown .and. &
                             grid%element(elem)%degree > 1
         if (refine_control%reftype /= P_ADAPTIVE .and. &
             refine_control%reftype /= HP_ADAPTIVE) pcoarsen_possible = .false.

! Compute coarsening error indicators

         if (hcoarsen_possible) then
            hcoarsen_errind = hcoarsen_indicator(grid,parent,siblings)
         else
            hcoarsen_errind = huge(0.0_my_real)
         endif
         if (pcoarsen_possible) then
            pcoarsen_errind = pcoarsen_indicator(grid,elem)
         else
            pcoarsen_errind = huge(0.0_my_real)
         endif

         hcoarsen_occured = .false.
         pcoarsen_occured = .false.

! If either coarsening indicator is small enough, perform the type of
! coarsening that has the least impact on the solution

         if (hcoarsen_errind < pcoarsen_errind .and. &
             hcoarsen_errind < small_errind) then

! find the next element that is not a sibling before doing h coarsening

            next_elem = grid%element(elem)%next
            do while (next_elem == siblings(1) .or. next_elem == siblings(2) &
                 .or. next_elem == siblings(3) .or. next_elem == siblings(4))
               next_elem = grid%element(next_elem)%next
            end do

! perform h coarsening

            call h_coarsen_element(grid,parent,errcode,refine_control)
            if (errcode == 0) then
               num_deref = num_deref + 1
               deref_list(num_deref) = parent
               degree_list(num_deref) = -1
               hcoarsen_occured = .true.
            endif

         elseif (pcoarsen_errind < small_errind) then

! perform p coarsening

            call p_coarsen_elem(grid,elem,errcode,refine_control)
            if (errcode == 0) then
               if (num_deref == 0) then
                  num_deref = num_deref + 1
                  deref_list(num_deref) = elem
                  degree_list(num_deref) = -2
               elseif (deref_list(num_deref) /= elem .or. &
                       degree_list(num_deref) /= -2) then
                  num_deref = num_deref + 1
                  deref_list(num_deref) = elem
                  degree_list(num_deref) = -2
               endif
               pcoarsen_occured = .true.
            endif
         endif


! if p coarsening occured, do the same element again.
! otherwise, next element

         if (.not. pcoarsen_occured) then
            if (grid%element(elem)%level /= 1) then
               hcoarsen_checked(siblings(1)) = .true.
               hcoarsen_checked(siblings(2)) = .true.
               if (.not. grid%element(parent)%mate == BOUNDARY) then
                  hcoarsen_checked(siblings(3)) = .true.
                  hcoarsen_checked(siblings(4)) = .true.
               endif
            endif
            if (hcoarsen_occured) then
               elem = next_elem
            else
               elem = grid%element(elem)%next
            endif
         endif
      end do ! next element
   end do ! next level

! see if enough derefinement has occured

   select case (refine_control%refterm)
   case (DOUBLE_NVERT, KEEP_NVERT)
      if (grid%nvert_own <= target) exit
   case (DOUBLE_NELEM, KEEP_NELEM)
      if (grid%nelem_leaf_own <= target) exit
   case (DOUBLE_NEQ, KEEP_NEQ)
      if (grid%dof_own <= target) exit
   case default
      exit
   end select
   small_errind = 2*small_errind
   if (small_errind > max_errind) exit

end do

! Send the lists of derefined elements to other processors and send the
! current degree of elements that were p-derefined and -1 to indicate
! h-derefinement

nproc = num_proc(procs)
my_processor = my_proc(procs)

nsend = (KEY_SIZE+1)*num_deref
allocate(isend(nsend),nrecv(nproc),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in derefinement",procs=procs)
   stop
endif

do i=1,num_deref
   call hash_pack_key(grid%element(deref_list(i))%gid,isend, &
                  1+(i-1)*(KEY_SIZE+1))
   if (degree_list(i) == -1) then
      isend(i*(KEY_SIZE+1)) = -1
   else
      isend(i*(KEY_SIZE+1)) = grid%element(deref_list(i))%degree
   endif
enddo

call phaml_alltoall(procs,isend,nsend,irecv,nrecv,3350)

! Derefine elements that I have but were derefined by the owner.
! Also note which ones cannot be derefined.

num_deref = 0
isub = 1
do p=1,nproc
   if (p == my_processor) then
      isub = isub + nrecv(p)
      cycle
   endif
   do i=1,nrecv(p)/(KEY_SIZE+1)
      elem = hash_decode_key(hash_unpack_key(irecv,isub),grid%elem_hash)
      if (elem /= HASH_NOT_FOUND) then
         if (irecv(isub+KEY_SIZE) == -1) then
            call h_coarsen_element(grid,elem,errcode,refine_control)
            if (errcode /= 0 .and. errcode /= -2) then
               num_deref = num_deref + 1
               deref_list(num_deref) = elem
               degree_list(num_deref) = -p
            endif
         else
            call p_coarsen_elem(grid,elem,errcode,refine_control)
            if (errcode /= 0) then
               num_deref = num_deref + 1
               deref_list(num_deref) = elem
               degree_list(num_deref) = 1
            endif
         endif
      endif
      isub = isub + KEY_SIZE+1
   end do
end do
if (associated(irecv)) deallocate(irecv)

! send the list of elements which could not be derefined.  For p refinement
! send the degree; for h refinement send -owner so only the owner rerefines

deallocate(isend)
nsend = (KEY_SIZE+1)*num_deref
allocate(isend(nsend),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in derefinement",procs=procs)
   stop
endif

do i=1,num_deref
   call hash_pack_key(grid%element(deref_list(i))%gid,isend, &
                      (1+(i-1)*(KEY_SIZE+1)))
   if (degree_list(i) < 0) then
      isend(i*(KEY_SIZE+1)) = degree_list(i)
   else
      isend(i*(KEY_SIZE+1)) = grid%element(deref_list(i))%degree
   endif
end do

call phaml_alltoall(procs,isend,nsend,irecv,nrecv,3351)

! rerefine elements that were returned as underefineable by some processor

isub = 1
do p=1,nproc
   do i=1,nrecv(p)/(KEY_SIZE+1)
      if (irecv(isub+KEY_SIZE) == -my_processor) then
         call create_element(grid,2*hash_unpack_key(irecv,isub), &
                             refine_control,solver_control,errcode)
      elseif (irecv(isub+KEY_SIZE) > 0) then
         elem = hash_decode_key(hash_unpack_key(irecv,isub),grid%elem_hash)
         if (elem /= HASH_NOT_FOUND) then
            if (grid%element(elem)%isleaf) then
               do while (grid%element(elem)%degree < irecv(isub+KEY_SIZE))
                  call p_refine_elem(grid,elem,refine_control)
               end do
            endif
         endif
      endif
      isub = isub + KEY_SIZE+1
   end do
end do
if (associated(irecv)) deallocate(irecv)

deallocate(isend,nrecv)

end subroutine derefinement

!        ------------------
function pcoarsen_indicator(grid,elem)
!        ------------------

!----------------------------------------------------
! This routine computes the p coarsening error indicator for element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real) :: pcoarsen_indicator
!----------------------------------------------------
! Local variables:

integer :: i, j, k, deg
real(my_real) :: temp
!----------------------------------------------------
! Begin executable code

! Compute the l2 norm of the p-hierarchic coefficients of the current
! degree of this element.  Don't worry about whether the edges are actually
! of the same degree.
! For systems or multiple eigenvectors, use the max of the individual solutions.

deg = grid%element(elem)%degree

if (deg <= 1) then
   pcoarsen_indicator = huge(0.0_my_real)
   return
endif

pcoarsen_indicator = 0

do j=1,grid%system_size
 do k=1,max(1,grid%num_eval)
   temp = 0
   do i = element_dof(deg-1)+1, element_dof(deg)
      temp = temp + grid%element(elem)%solution(i,j,k)**2
   end do
   do i=1,EDGES_PER_ELEMENT
      if (grid%edge(grid%element(elem)%edge(i))%degree >= deg) then
         temp = temp + grid%edge(grid%element(elem)%edge(i))%solution(deg-1,j,k)**2
      endif
   end do
   pcoarsen_indicator = max(pcoarsen_indicator,sqrt(temp))
 end do
end do

end function pcoarsen_indicator

!        ------------------
function hcoarsen_indicator(grid,parent,children)
!        ------------------

!----------------------------------------------------
! This routine computes the h coarsening error indicator for the elements in
! children, which should be siblings that are children of parent and its mate.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: parent,children(:)
real(my_real) :: hcoarsen_indicator
!----------------------------------------------------
! Local variables:

real(my_real), allocatable :: solution(:,:)
real(my_real) :: xvert(4),yvert(4),temp
integer :: i,j,k,degree,astat,isub
!----------------------------------------------------
! Begin executable code

! set degree to be the maximum degree, and allocate space for the local
! copy of the solution coefficients

degree = 0
do i=1,4
  if (i==3 .and. children(3)==NO_CHILD) exit
  do j=1,EDGES_PER_ELEMENT
     degree = max(degree,grid%edge(grid%element(children(i))%edge(j))%degree)
  end do
end do

allocate(solution((degree+1)**2+degree**2,grid%nsoln),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in hcoarsen_indicator")
   return
endif

! copy solution coefficients to the local variable.  see subroutine
! phier2nodal for the order of the solution components

solution = 0
solution(1,:) = reshape( &
                grid%vertex_solution(grid%element(children(1))%vertex(1),:,:), &
                (/grid%nsoln/))
solution(2,:) = reshape( &
                grid%vertex_solution(grid%element(children(2))%vertex(1),:,:), &
                (/grid%nsoln/))
solution(3,:) = reshape( &
                grid%vertex_solution(grid%element(children(1))%vertex(2),:,:), &
                (/grid%nsoln/))
if (children(3) /= NO_CHILD) then
   solution(4,:) = reshape( &
                grid%vertex_solution(grid%element(children(3))%vertex(2),:,:), &
                (/grid%nsoln/))
endif
isub = 5
do i=1,degree-1
   if (grid%edge(grid%element(children(1))%edge(1))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(1))%edge(1))%solution(i,:,:),(/grid%nsoln/))
   endif
   isub = isub+1
end do
solution(isub,:) = reshape( &
                grid%vertex_solution(grid%element(children(1))%vertex(3),:,:), &
                (/grid%nsoln/))
isub = isub+1
do i=1,degree-1
   if (grid%edge(grid%element(children(1))%edge(2))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(1))%edge(2))%solution(i,:,:),(/grid%nsoln/))
   endif
   isub = isub + 1
end do
do i=1,degree-1
   if (grid%edge(grid%element(children(1))%edge(3))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(1))%edge(3))%solution(i,:,:),(/grid%nsoln/))
   endif
   isub = isub + 1
end do
k = 1
do j=1,degree-2
   do i=1,j
      if (grid%element(children(1))%degree >= j+2) then
         solution(isub,:) = reshape(grid%element(children(1))%solution(k,:,:),(/grid%nsoln/))
      endif
      isub = isub + 1
      k = k+1
   end do
end do
do i=1,degree-1
   if (grid%edge(grid%element(children(2))%edge(2))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(2))%edge(2))%solution(i,:,:),(/grid%nsoln/))
   endif
   isub = isub + 1
end do
do i=1,degree-1
   if (grid%edge(grid%element(children(2))%edge(3))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(2))%edge(3))%solution(i,:,:),(/grid%nsoln/))
   endif
   isub = isub + 1
end do
k = 1
do j=1,degree-2
   do i=1,j
      if (grid%element(children(2))%degree >= j+2) then
         solution(isub,:) = reshape(grid%element(children(2))%solution(k,:,:),(/grid%nsoln/))
      endif
      isub = isub + 1
      k = k+1
   end do
end do
if (children(3) /= NO_CHILD) then
   do i=1,degree-1
      if (grid%edge(grid%element(children(3))%edge(1))%degree >= i+1) then
         solution(isub,:) = reshape(grid%edge(grid%element(children(3))%edge(1))%solution(i,:,:),(/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   do i=1,degree-1
      if (grid%edge(grid%element(children(3))%edge(3))%degree >= i+1) then
         solution(isub,:) = reshape(grid%edge(grid%element(children(3))%edge(3))%solution(i,:,:),(/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   k = 1
   do j=1,degree-2
      do i=1,j
         if (grid%element(children(3))%degree >= j+2) then
            solution(isub,:) = reshape(grid%element(children(3))%solution(k,:,:),(/grid%nsoln/))
         endif
         isub = isub + 1
         k = k+1
      end do
   end do
   do i=1,degree-1
      if (grid%edge(grid%element(children(4))%edge(3))%degree >= i+1) then
         solution(isub,:) = reshape(grid%edge(grid%element(children(4))%edge(3))%solution(i,:,:),(/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   k = 1
   do j=1,degree-2
      do i=1,j
         if (grid%element(children(4))%degree >= j+2) then
            solution(isub,:) = reshape(grid%element(children(4))%solution(k,:,:),(/grid%nsoln/))
         endif
         isub = isub + 1
         k = k+1
      end do
   end do
endif

! set the outer vertices of the parents

xvert(1) = grid%vertex(grid%element(parent)%vertex(1))%coord%x
yvert(1) = grid%vertex(grid%element(parent)%vertex(1))%coord%y
xvert(2) = grid%vertex(grid%element(parent)%vertex(2))%coord%x
yvert(2) = grid%vertex(grid%element(parent)%vertex(2))%coord%y
xvert(3) = grid%vertex(grid%element(parent)%vertex(3))%coord%x
yvert(3) = grid%vertex(grid%element(parent)%vertex(3))%coord%y
if (children(3) == NO_CHILD) then
   xvert(4) = huge(0.0_my_real)
   yvert(4) = huge(0.0_my_real)
else
   xvert(4) = grid%vertex(grid%element(children(3))%vertex(2))%coord%x
   yvert(4) = grid%vertex(grid%element(children(3))%vertex(2))%coord%y
endif

! convert solution to nodal basis coefficients, and then to h-hierarchical
! basis coefficients

call phier2nodal(solution,xvert,yvert,degree)
call nodal2hhier(solution,xvert,yvert,degree)

! compute the l2 norm of the red node h-hierarchic coefficients; for systems
! and multiple eigenvectors use the max over all solutions
! See subroutine phier2nodal for the order of the nodes

hcoarsen_indicator = 0

do j=1,grid%nsoln
   temp = 0
! along edge between children(1) and children(2), including central node
   isub = 5
   do i=isub,isub+2*((degree-1)/2),2
      temp = temp + solution(i,j)**2
   end do
   if (degree > 1) then
! along edge between children(1) and children(3)
      isub = degree+5
      do i=isub,isub+2*((degree-2)/2),2
         temp = temp + solution(i,j)**2
      end do
! along edge between children(2) and children(4)
      isub = 3+3*degree + element_dof(degree)
      do i=isub,isub+2*((degree-2)/2),2
         temp = temp + solution(i,j)**2
      end do
! along edge between children(3) and children(4)
      if (children(3) /= NO_CHILD) then
         isub = 1 + 5*degree + 2*element_dof(degree)
         do i=isub,isub+2*((degree-2)/2),2
            temp = temp + solution(i,j)**2
         end do
      endif
      if (degree > 2) then
! interior of children(1)
         isub = 6 + 3*(degree-1)
         do k=1,(degree-1)/2
            do i=1,degree-2*k
               temp = temp + solution(isub,j)**2
               isub = isub + 1
            end do
            isub = isub + degree - 2*k - 1
         end do
! interior of children(2)
         isub = 1 + 5*degree + element_dof(degree)
         do k=1,(degree-1)/2
            do i=1,degree-2*k
               temp = temp + solution(isub,j)**2
               isub = isub + 1
            end do
            isub = isub + degree - 2*k - 1
         end do
         if (children(3) /= NO_CHILD) then
! interior of children(3)
            isub = -1 + 7*degree + 2*element_dof(degree)
            do k=1,(degree-1)/2
               do i=1,degree-2*k
                  temp = temp + solution(isub,j)**2
                  isub = isub + 1
               end do
               isub = isub + degree - 2*k - 1
            end do
! interior of children(4)
            isub = -2 + 8*degree + 3*element_dof(degree)
            do k=1,(degree-1)/2
               do i=1,degree-2*k
                  temp = temp + solution(isub,j)**2
                  isub = isub + 1
               end do
               isub = isub + degree - 2*k - 1
            end do
         endif ! children(3) /= NO_CHILD
      endif ! degree > 2
   endif ! degree > 1

   hcoarsen_indicator = max(hcoarsen_indicator,sqrt(temp))
end do

deallocate(solution)

end function hcoarsen_indicator

!          ----------
subroutine make_elist(grid,procs,elist,refcont,global_max_errind, &
                      numhref,numpref,predictive_load_balance,still_sequential)
!          ----------

!----------------------------------------------------
! This routine groups the leaf elements into size(head_errind) lists.
! List 1 contains those with error indicator between maxerrind and
! maxerrind/binw, list 2 between maxerrind/binw and maxerrind/binw**2, etc.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type(errind_list), intent(inout) :: elist
type(refine_options), intent(in) :: refcont
real(my_real), intent(in) :: global_max_errind
integer, intent(in) :: numhref(:), numpref(:)
logical, intent(in) :: predictive_load_balance, still_sequential

!----------------------------------------------------
! Local variables:

real(my_real), allocatable :: errind(:)
integer :: lev, elem, i, astat, total
real(my_real) :: max_errind, reftol, normsoln
logical :: is_uniform

!----------------------------------------------------
! Begin executable code

allocate(errind(grid%biggest_elem),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in make_elist",procs=procs)
   return
endif
errind = 0.0_my_real

! set error indicators and find maximum

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%iown) then
         if (grid%element(elem)%isleaf) then
            errind(elem) = &
                  maxval(grid%element_errind(elem,:))/grid%element(elem)%work
         endif
      endif
      elem = grid%element(elem)%next
   end do
end do

! use global max errind from before derefinement for computing bins

elist%max_errind = global_max_errind
max_errind = global_max_errind

if (refcont%refterm == ONE_REF) then
   call get_grid_info(grid,procs,still_sequential,3344, &
                      total_nelem_leaf=total,no_master=.true.)
   reftol = refcont%reftol/sqrt(real(total,my_real))
! TEMP only using first component of first eigenvalue
   if (grid%errtype == RELATIVE_ERROR .and. .not. (refcont%reftype == HP_ADAPTIVE .and. &
    (refcont%hp_strategy == HP_T3S .or. refcont%hp_strategy == HP_ALTERNATE))) then
      call norm_solution(grid,procs,still_sequential,1,1,energy=normsoln)
      if (normsoln /= 0.0_my_real) reftol = reftol*normsoln
   endif
endif
if (refcont%refterm == ONE_REF_HALF_ERRIND) then
   reftol = global_max_errind/refcont%inc_factor
endif

! create lists

elist%head_errind = END_OF_LIST
elist%tail_errind = END_OF_LIST
elist%next_errind = NOT_ON_LIST
elist%prev_errind = NOT_ON_LIST

! determine if this is uniform refinement

is_uniform = .false.
if (global_max_errind == 0.0_my_real) is_uniform = .true.
if (refcont%reftype == HP_ADAPTIVE .and. &
    refcont%hp_strategy == HP_T3S) then
   if (refcont%t3s_reftype == H_UNIFORM) is_uniform = .true.
endif

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
! if numhref or numpref is positive, it goes in the first bin
         if (numhref(elem) > 0 .or. numpref(elem) > 0) then
            i = 1
! if both numhref and numpref are 0, it goes in the last bin
         elseif (numhref(elem) == 0 .and. numpref(elem) == 0) then
            i = size(elist%head_errind)
! uniform refinement puts everything in first bin
         elseif (is_uniform) then
            i = 1
         else
            select case (refcont%refterm)
! ONE_REF put large error indicators in first bin and everything else in second
            case (ONE_REF, ONE_REF_HALF_ERRIND)
               if (grid%element(elem)%iown .and. &
                   maxval(grid%element_errind(elem,:)) > reftol) then
                  i = 1
               else
                  i = 2
               endif
            case default
! otherwise find the right bin
               do i=1,size(elist%head_errind)-1
                  if (errind(elem) > max_errind/(binw**i)) exit
               end do
            end select
         endif
         if (elist%head_errind(i) == END_OF_LIST) then
            elist%head_errind(i) = elem
            elist%tail_errind(i) = elem
            elist%next_errind(elem) = END_OF_LIST
            elist%prev_errind(elem) = END_OF_LIST
         else
            elist%next_errind(elem) = elist%head_errind(i)
            elist%prev_errind(elist%head_errind(i)) = elem
            elist%head_errind(i) = elem
            elist%prev_errind(elem) = END_OF_LIST
         endif
      endif
      elem = grid%element(elem)%next
   end do
end do
deallocate(errind,stat=astat)

end subroutine make_elist

!          ----------------
subroutine set_loop_control(refine_control,return_to_elist,complete_elist, &
                            one_elist)
!          ----------------

!----------------------------------------------------
! This routine sets the variables that control the refine loop:
! return_to_elist, tells if children get put into an elist
! complete_elist, tells if you must finish the whole elist
! one_elist, tells if you quit after finishing the first elist
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(refine_options), intent(in) :: refine_control
logical, intent(out) :: return_to_elist, complete_elist, one_elist
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

select case (refine_control%reftype)
case (H_ADAPTIVE, P_ADAPTIVE)
   select case (refine_control%refterm)
   case (DOUBLE_NVERT, DOUBLE_NELEM, DOUBLE_NEQ, HALVE_ERREST, KEEP_NVERT, &
         KEEP_NELEM, KEEP_NEQ, KEEP_ERREST)
      return_to_elist = .true.
      complete_elist = .false.
      one_elist = .false.
   case (ONE_REF, ONE_REF_HALF_ERRIND)
      return_to_elist = .false.
      complete_elist = .true.
      one_elist = .true.
   end select
case (HP_ADAPTIVE)
   select case (refine_control%hp_strategy)
   case (HP_BIGGER_ERRIND, HP_APRIORI, HP_PRIOR2P_E, HP_PRIOR2P_H1, &
         HP_TYPEPARAM, HP_COEF_DECAY, HP_COEF_ROOT, HP_SMOOTH_PRED, &
         HP_NEXT3P, HP_STEEPEST_SLOPE)
      select case (refine_control%refterm)
      case (DOUBLE_NVERT, DOUBLE_NELEM, DOUBLE_NEQ, HALVE_ERREST, KEEP_NVERT, &
            KEEP_NELEM, KEEP_NEQ, KEEP_ERREST)
         return_to_elist = .true.
         complete_elist = .false.
         one_elist = .false.
      case (ONE_REF, ONE_REF_HALF_ERRIND)
         return_to_elist = .false.
         complete_elist = .true.
         one_elist = .true.
      end select
   case (HP_T3S)
      return_to_elist = .true.
      complete_elist = .true.
      one_elist = .true.
   case (HP_ALTERNATE)
      return_to_elist = .false.
      complete_elist = .true.
      one_elist = .true.
   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("unexpected value for hp_strategy in set_loop_control")
      stop
   end select
case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unexpected value for reftype in set_loop_control")
   stop
end select

end subroutine set_loop_control

!        ------------
function check_target(grid,refine_control,elist,target)
!        ------------

!----------------------------------------------------
! This routine checks to see if the target for terminating refinement is met
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(refine_options), intent(in) :: refine_control
type(errind_list), intent(in) :: elist
integer, intent(in) :: target
logical :: check_target
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

check_target = .false.

select case (refine_control%refterm)
case (DOUBLE_NVERT, KEEP_NVERT)
   if (grid%nvert_own >= target) check_target = .true.
case (DOUBLE_NELEM, KEEP_NELEM)
   if (grid%nelem_leaf_own >= target) check_target = .true.
case (DOUBLE_NEQ, KEEP_NEQ)
   if (grid%dof_own >= target) check_target = .true.
case (HALVE_ERREST, KEEP_ERREST)
   if (elist%current_list > target) check_target = .true.
case (ONE_REF, ONE_REF_HALF_ERRIND)
   check_target = .false.
case default
   call fatal("illegal value for refterm")
   stop
end select

end function check_target

!          ---------
subroutine reconcile(grid,procs,refine_control,solver_control,sequent)
!          ---------

!----------------------------------------------------
! This routine reconciles the refinements among the partitions
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
logical, intent(in) :: sequent

!----------------------------------------------------
! Local variables:

integer :: lev, elem, count, part, i, j, elemlid, errcode, vertlid, &
           astat, nproc, vert, my_master, refedge
! newcomm
integer :: nsend
integer, allocatable :: isend(:), nrecv(:)
integer, pointer :: irecv(:)
type(hash_key) :: gid
logical, pointer :: isboundary(:), boundary_elem(:)
logical :: doit, isboundary_periodic

!----------------------------------------------------
! Begin executable code

if (sequent .or. num_proc(procs) == 1) return

grid_changed = .true.

if (my_proc(procs) == MASTER) return

! start timing the reconciliation process

call reset_watch((/precon,cprecon/))
call start_watch((/precon,trecon/))

nproc = num_proc(procs)

allocate(nrecv(nproc),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reconcile",procs=procs)
   return
endif

! make a list of the unowned elements I have refined along with the degree
! if it was p refined or -1 if it was h refined

count = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%hrefined_unowned) count = count + 1
      if (grid%element(elem)%prefined_unowned) count = count + 1
      elem = grid%element(elem)%next
   end do
end do

nsend = (1+KEY_SIZE)*count
allocate(isend(nsend),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reconcile",procs=procs)
   return
endif

count = 1
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%hrefined_unowned) then
         call hash_pack_key(grid%element(elem)%gid,isend,count)
         isend(count+KEY_SIZE) = -1
         grid%element(elem)%hrefined_unowned = .false.
         count = count + KEY_SIZE+1
      endif
      if (grid%element(elem)%prefined_unowned) then
         call hash_pack_key(grid%element(elem)%gid,isend,count)
         isend(count+KEY_SIZE) = grid%element(elem)%degree
         grid%element(elem)%prefined_unowned = .false.
         count = count + KEY_SIZE+1
      endif
      elem = grid%element(elem)%next
   end do
end do
 
! exchange the list with other partitions

call start_watch((/cprecon,ctrecon/))
call phaml_alltoall(procs,isend,nsend,irecv,nrecv,3470)
call stop_watch((/cprecon,ctrecon/))

errcode = 0
deallocate(isend,stat=astat)

! scan the list looking for elements that I own and have not refined,
! and refine those elements

count = 1
outer: do part=1,nproc
   do i=1,nrecv(part)/(KEY_SIZE+1)
      elemlid = hash_decode_key(hash_unpack_key(irecv,count),grid%elem_hash)
      if (elemlid == HASH_NOT_FOUND) then
         count = count + KEY_SIZE+1
         cycle
      endif
      if (.not. grid%element(elemlid)%iown) then
         count = count + KEY_SIZE+1
         cycle
      endif
      if (irecv(count+KEY_SIZE) == -1) then
         if (grid%element(elemlid)%isleaf) then
            call h_refine_element(grid,elemlid,errcode,refedge,refine_control, &
                                  solver_control)
            if (errcode /= 0) then
               call warning("refinement failed during reconciliation")
               exit outer
            endif
         endif
      else
         if (grid%element(elemlid)%isleaf) then
            do while (grid%element(elemlid)%degree < irecv(count+KEY_SIZE))
               call p_refine_elem(grid,elemlid,refine_control)
            end do
         endif
      endif
      count = count + KEY_SIZE+1
   end do
end do outer

if (associated(irecv)) deallocate(irecv,stat=astat)

! enforce overlap so that any element that touches my partition from the
! outside exists and has high enough degrees to match the one that owns it

! mark each vertex as being on the partition boundary or not by going through
! the elements and marking each vertex as being on the boundary if I own the
! element or vertex but not the other

allocate(isboundary(size(grid%vertex)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reconcile",procs=procs)
   return
endif
isboundary = .false.

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do i=1,VERTICES_PER_ELEMENT
         if (grid%element(elem)%iown .neqv. &
             grid%element(grid%vertex(grid%element(elem)%vertex(i))%assoc_elem)%iown) then
            isboundary(grid%element(elem)%vertex(i)) = .true.
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

! periodic boundary vertices should be isboundary if either the master or
! any slave is

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do i=1,VERTICES_PER_ELEMENT
         vert = grid%element(elem)%vertex(i)
         if (.not. is_periodic_vert(vert,grid)) cycle
         isboundary_periodic = isboundary(vert)
         my_master = grid%vertex(vert)%next
         do while (my_master /= vert)
            isboundary_periodic = isboundary_periodic .or. isboundary(my_master)
            my_master = grid%vertex(my_master)%next
         end do
         if (isboundary_periodic) then
            isboundary(vert) = .true.
            my_master = grid%vertex(vert)%next
            do while (my_master /= vert)
               isboundary(my_master) = .true.
               my_master = grid%vertex(my_master)%next
            end do
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

! make a list of owned leaves that touch the partition boundary by going
! through the elements and finding those that have a marked vertex.  Send
! the element along with the vertices and the element and edge degrees

count = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
         do i=1,VERTICES_PER_ELEMENT
            if (isboundary(grid%element(elem)%vertex(i))) then
               count = count + 1
               exit
            endif
         end do
      endif
      elem = grid%element(elem)%next
   end do
end do

nsend=4*count*(KEY_SIZE+1)
allocate(isend(nsend),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reconcile",procs=procs)
   return
endif

count = 1
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
         do i=1,VERTICES_PER_ELEMENT
            if (isboundary(grid%element(elem)%vertex(i))) then
               call hash_pack_key(grid%element(elem)%gid,isend,count)
               call hash_pack_key(grid%vertex(grid%element(elem)%vertex(1))%gid, &
                                  isend,count+KEY_SIZE)
               call hash_pack_key(grid%vertex(grid%element(elem)%vertex(2))%gid, &
                                  isend,count+2*KEY_SIZE)
               call hash_pack_key(grid%vertex(grid%element(elem)%vertex(3))%gid, &
                                  isend,count+3*KEY_SIZE)
               isend(count+4*KEY_SIZE) = grid%element(elem)%degree
               isend(count+4*KEY_SIZE+1:count+4*KEY_SIZE+3) = &
                                  grid%edge(grid%element(elem)%edge)%degree
               count = count + 4*(KEY_SIZE+1)
               exit
            endif
         end do
      endif
      elem = grid%element(elem)%next
   end do
end do

! exchange lists with other processors

call start_watch((/cprecon,ctrecon/))
call phaml_alltoall(procs,isend,nsend,irecv,nrecv,3480)
call stop_watch((/cprecon,ctrecon/))

deallocate(isend,stat=astat)

! go through the received lists of elements.  If any vertex of the element is
! marked as being on the partition boundary, h refine to create the element if
! it doesn't exist, and p refine if the element or edge degrees are too small.
! Also, keep track of which ones are on the partition boundary for later.

allocate(boundary_elem(size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reconcile",procs=procs)
   return
endif
boundary_elem = .false.
count = 1
do part=1,nproc
   do i=1,nrecv(part)/(4*(KEY_SIZE+1))
      doit = .false.
      do j=1,VERTICES_PER_ELEMENT
         vertlid = hash_decode_key(hash_unpack_key(irecv,count+j*KEY_SIZE), &
                                   grid%vert_hash)
         if (vertlid == HASH_NOT_FOUND) cycle
         if (isboundary(vertlid)) then
            doit = .true.
            exit
         endif
      end do
      if (doit) then
         gid = hash_unpack_key(irecv,count)
         elemlid = hash_decode_key(gid,grid%elem_hash)
         if (elemlid == HASH_NOT_FOUND) then
            call create_element(grid,gid,refine_control,solver_control,errcode)
            if (errcode /= 0) then
               call warning("refinement for overlap failed")
               count = count + 4*(KEY_SIZE+1)
               cycle
            endif
            if (size(isboundary) /= size(grid%vertex)) then
               call realloc_isboundary(isboundary,size(grid%vertex))
            endif
            if (size(boundary_elem) /= size(grid%element)) then
               call realloc_boundary_elem(boundary_elem,size(grid%element))
            endif
            elemlid = hash_decode_key(gid,grid%elem_hash)
         endif
         if (grid%element(elemlid)%isleaf) then
            boundary_elem(elemlid) = .true.
            do while (grid%element(elemlid)%degree < irecv(count+4*KEY_SIZE))
               call p_refine_elem(grid,elemlid,refine_control)
            end do
         endif
         do j=1,EDGES_PER_ELEMENT
            if (grid%edge(grid%element(elemlid)%edge(j))%degree < &
                irecv(count+4*KEY_SIZE+j)) then
               call adjust_edge_degree(grid,grid%element(elemlid)%edge(j), &
                                       irecv(count+4*KEY_SIZE+j))
            endif
         end do
      endif
      count = count + 4*(KEY_SIZE+1)
   end do
end do

if (associated(irecv)) deallocate(irecv,stat=astat)
deallocate(isboundary,stat=astat)

deallocate(nrecv,stat=astat)

! if doing adaptive p, p coarsen each unowned leaf as much as possible
! TEMP I may be p coarsening elements that need to be p refined for overlap,
!      which adds extra work and looses the solution

if (refine_control%reftype == P_ADAPTIVE .or. &
    refine_control%reftype == HP_ADAPTIVE) then
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf .and. &
             (.not. grid%element(elem)%iown) .and. &
             (.not. boundary_elem(elem))) then
            errcode = 0
            do while (errcode == 0)
               call p_coarsen_elem(grid,elem,errcode,refine_control)
            end do
         endif
         elem = grid%element(elem)%next
      end do
   end do
endif

deallocate(boundary_elem)
grid%errind_up2date = .false.

! stop timing the reconciliation process

call stop_watch((/precon,trecon/))

end subroutine reconcile

!          ------------------
subroutine realloc_isboundary(isboundary,newsize)
!          ------------------

!----------------------------------------------------
! This routine reallocates isboundary to newsize
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical, pointer :: isboundary(:)
integer, intent(in) :: newsize
!----------------------------------------------------
! Local variables:

logical, pointer :: newvar(:)
integer :: astat
!----------------------------------------------------
! Begin executable code

allocate(newvar(newsize),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in realloc_isboundary")
   stop
endif
newvar(1:size(isboundary)) = isboundary
newvar(size(isboundary)+1:newsize) = .false.
deallocate(isboundary)
isboundary => newvar

end subroutine realloc_isboundary

!          ------------------
subroutine realloc_boundary_elem(boundary_elem,newsize)
!          ------------------

!----------------------------------------------------
! This routine reallocates boundary_elem to newsize
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical, pointer :: boundary_elem(:)
integer, intent(in) :: newsize
!----------------------------------------------------
! Local variables:

logical, pointer :: newvar(:)
integer :: astat
!----------------------------------------------------
! Begin executable code

allocate(newvar(newsize),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in realloc_boundary_elem")
   stop
endif
newvar(1:size(boundary_elem)) = boundary_elem
newvar(size(boundary_elem)+1:newsize) = .false.
deallocate(boundary_elem)
boundary_elem => newvar

end subroutine realloc_boundary_elem

!          ------------------
subroutine adjust_edge_degree(grid,edge,newdeg)
!          ------------------

!----------------------------------------------------
! This routine sets the edge degree to the given value.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: edge, newdeg
!----------------------------------------------------
! Local variables:

integer :: olddeg, edge_size, oldsize, astat, d1, d2, d3, i
real(my_real), pointer :: temp1(:,:,:)
!----------------------------------------------------
! Begin executable code

olddeg = grid%edge(edge)%degree

! make sure allocated memory is large enough.

edge_size = newdeg-1
if (associated(grid%edge(edge)%solution)) then
   oldsize = size(grid%edge(edge)%solution,dim=1)
else
   oldsize = 0
endif
if (oldsize < edge_size) then
   allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)), &
            stat=astat)
   if (astat /= 0) then
      call fatal("allocation failed in adjust_edge_degree")
      stop
   endif
   temp1 = 0.0_my_real
   if (oldsize > 0) temp1(1:oldsize,:,:) = grid%edge(edge)%solution
   deallocate(grid%edge(edge)%solution, stat=astat)
   grid%edge(edge)%solution => temp1
   if (grid%have_true) then
      nullify(temp1)
      allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)), &
               stat=astat)
      if (astat /= 0) then
         call fatal("allocation failed in adjust_edge_degree")
         stop
      endif
      temp1 = 0.0_my_real
      if (oldsize > 0) temp1(1:oldsize,:,:) = grid%edge(edge)%exact
      deallocate(grid%edge(edge)%exact, stat=astat)
      grid%edge(edge)%exact => temp1
   endif
endif
if (grid%oldsoln_exists) then
   if (associated(grid%edge(edge)%oldsoln)) then
      nullify(temp1)
      allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in adjust_edge_degree")
         stop
      endif
      temp1 = 0.0_my_real
      d1 = min(size(grid%edge(edge)%oldsoln,dim=1),size(temp1,dim=1))
      d2 = min(size(grid%edge(edge)%oldsoln,dim=2),size(temp1,dim=2))
      d3 = min(size(grid%edge(edge)%oldsoln,dim=3),size(temp1,dim=3))
      temp1(1:d1,1:d2,1:d3) = grid%edge(edge)%oldsoln(1:d1,1:d2,1:d3)
      deallocate(grid%edge(edge)%oldsoln)
      grid%edge(edge)%oldsoln => temp1
   endif
endif

! Set the edge degree.  Also set Dirichlet boundary conditions on Dirichlet
! edges that change degree and set solution to 0 for other new components

grid%edge(edge)%degree = newdeg
grid%dof = grid%dof + grid%system_size*(newdeg - olddeg)
if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
   grid%dof_own = grid%dof_own + grid%system_size*(newdeg - olddeg)
endif
if (grid%have_true) then
   grid%edge(edge)%exact = 0.0_my_real
endif
do i=1,grid%system_size
   if (grid%edge_type(edge,i) == DIRICHLET) then
      call edge_exact(grid,edge,i,"d")
   else
      grid%edge(edge)%solution(newdeg-1,i,:) = 0.0_my_real
   endif
   if (grid%have_true) call edge_exact(grid,edge,i,"t")
end do

end subroutine adjust_edge_degree

end module refine_adapt_mod
