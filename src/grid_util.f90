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

module grid_util

!----------------------------------------------------
! This module contains grid utility routines.
! communication tags in this module are of the form 35xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use hash_mod
use gridtype_mod
use quadrature_rules
use basis_functions
use message_passing
!----------------------------------------------------

implicit none
private
public allocate_grid, &        ! not thread safe
       deallocate_grid, &      ! not thread safe
       realloc_solution, &     ! not thread safe; TEMP OPENMP can be
       reset_dirich_exact, &   ! not thread safe; TEMP OPENMP can be
       copy_grid, &            ! not thread safe
       copy_old, &             ! not thread safe
       edge_exact, &           ! thread safe if vertex not changing
       face_exact, &           ! thread safe
       elem_exact, &           ! thread safe if vertex and edge not changing
       point_on_edge, &        ! thread safe
       level2_mate, &          ! thread safe
       matching_periodic_vert,&! thread safe
       matching_periodic_edge,&! thread safe
       is_periodic_vert, &     ! thread safe
       is_periodic_edge, &     ! thread safe
       is_periodic_face, &     ! thread safe
       is_periodic_vert_master,& ! thread safe
       is_periodic_edge_master,& ! thread safe
       is_periodic_face_master,& ! thread safe
       is_periodic_vert_slave, & ! thread safe
       is_periodic_edge_slave, & ! thread safe
       is_periodic_face_slave, & ! thread safe
       find_containing_element,&!thread safe
       get_child_lid, &        ! thread safe
       get_child_gid, &        ! thread safe
       get_neighbors, &        ! thread safe
       get_vertex_elements, &  ! thread safe
       get_edge_elements, &    ! thread safe
       get_face_elements, &    ! thread safe
       get_grid_info, &        ! not thread safe
       get_vertices_and_solution, & ! thread safe
       check_pde_terms, &      ! thread safe
       check_triangle_shapes, &! thread safe
       get_next_free_vert, &   ! not thread save
       get_next_free_edge, &   ! not thread save
       get_next_free_elem, &   ! not thread save
       put_next_free_vert, &   ! not thread save
       put_next_free_edge, &   ! not thread save
       put_next_free_elem, &   ! not thread save
       extend_nlev, &          ! not thread safe
       more_verts, &           ! not thread safe
       more_edges, &           ! not thread safe
       more_elements, &        ! not thread safe
       phier2nodal, &          ! thread safe
       nodal2phier, &          ! thread safe
       nodal2hhier, &          ! thread safe
       element_diameter, &     ! thread safe
       element_volume, &       ! thread safe
       list_elements, &        ! not thread safe
       list_edges_without_rule, & ! not thread safe
       element_dof, &          ! thread safe
       face_dof, &             ! thread safe
       count_dof, &            ! not thread safe
       compute_global_max_errind, & ! not thread safe; is OpenMP parallel
       zcoord, &               ! thread safe
       set_zcoord, &           ! thread safe
       my_pdecoefs, &          ! thread safe
       my_bconds, &            ! thread safe
       my_iconds, &            ! thread safe
       my_trues, &             ! thread safe
       my_truexs, &            ! thread safe
       my_trueys, &            ! thread safe
       my_truezs, &            ! thread safe
       my_boundary_point, &    ! thread safe
       my_boundary_npiece, &   ! thread safe
       my_phaml_integral_kernel, & ! thread safe
       my_regularity           ! thread safe

!----------------------------------------------------
! Non-module procedures used are:

interface

   subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
   use global
   real(my_real), intent(in) :: x,y
   real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                 c(:,:),rs(:)
   end subroutine pdecoefs

   subroutine bconds(x,y,bmark,itype,c,rs)
   use global
   real(my_real), intent(in) :: x,y
   integer, intent(in) :: bmark
   integer, intent(out) :: itype(:)
   real(my_real), intent(out) :: c(:,:),rs(:)
   end subroutine bconds

   function iconds(x,y,comp,eigen)
   use global
   real(my_real), intent(in) :: x,y
   integer, intent(in) :: comp, eigen
   real(my_real) :: iconds
   end function iconds

   function trues(x,y,comp,eigen) ! real (my_real)
   use global
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp, eigen
   real (my_real) :: trues
   end function trues

   function truexs(x,y,comp,eigen) ! real (my_real)
   use global
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp, eigen
   real (my_real) :: truexs
   end function truexs

   function trueys(x,y,comp,eigen) ! real (my_real)
   use global
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp, eigen
   real (my_real) :: trueys
   end function trueys

   subroutine boundary_point(ipiece,s,x,y)
   use global
   integer, intent(in) :: ipiece
   real(my_real), intent(in) :: s
   real(my_real), intent(out) :: x,y
   end subroutine boundary_point

   function boundary_npiece(hole)
   integer, intent(in) :: hole
   integer :: boundary_npiece
   end function boundary_npiece

   subroutine boundary_param(start,finish)
   use global
   real(my_real), intent(out) :: start(:), finish(:)
   end subroutine boundary_param

   function phaml_integral_kernel(kernel,x,y)
   use global
   integer, intent(in) :: kernel
   real(my_real), intent(in) :: x,y
   real(my_real) :: phaml_integral_kernel
   end function phaml_integral_kernel

   function regularity(x,y)
   use global
   real(my_real), intent(in) :: x(*),y(*)
   real(my_real) :: regularity
   end function regularity

end interface

!----------------------------------------------------
! Generic procedures

private get_child_lid_scalar, get_child_lid_array, &
        get_child_gid_scalar, get_child_gid_array, &
        zcoord_scalar, zcoord_array

interface get_child_lid
   module procedure get_child_lid_scalar, get_child_lid_array
end interface

interface get_child_gid
   module procedure get_child_gid_scalar, get_child_gid_array
end interface

interface zcoord
   module procedure zcoord_scalar, zcoord_array
end interface

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

! Undocumented feature to make the grid and current element lid available to
! pde_coefs to, for example, access tags.  But don't tell anyone!
! pde_coefs need to use both global and grid_util.

type(grid_type), save, pointer, public :: secret_grid

!----------------------------------------------------

contains

!          -------------
subroutine allocate_grid(grid,nvert,nev,type,degree)
!          -------------

!----------------------------------------------------
! This routine performs initial allocation of memory in the grid data structure
! and initializes it with an empty grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout), target :: grid
integer, intent(in) :: nvert, nev, type, degree

!----------------------------------------------------
! Local variables:

integer :: i, astat1, astat2, astat3, dstat, npiece, nevdim

!----------------------------------------------------
! Begin executable code

nevdim = max(1,nev)

! grid data structure

allocate(grid%element(4*nvert),grid%edge(4*nvert),grid%vertex(nvert), &
         grid%face(0),stat=astat1)
allocate(grid%errest_energy(nevdim), &
         grid%errest_Linf(grid%system_size*nevdim), &
         grid%errest_L2(grid%system_size*nevdim), &
         grid%errest_eigenvalue(nevdim),stat=astat2)
if (astat1 /= 0 .or. astat2 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf, grid%errest_L2,grid%errest_eigenvalue, &
              stat=dstat)
   return
endif
nullify(grid%initial_neighbor)
call hash_table_init(grid%elem_hash,size(grid%element))
if (ierr == ALLOC_FAILED) then
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2, &
              grid%errest_eigenvalue,stat=dstat)
   return
endif
call hash_table_init(grid%face_hash,1)
if (ierr == ALLOC_FAILED) then
   call hash_table_destroy(grid%elem_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2, &
              grid%errest_eigenvalue,stat=dstat)
   return
endif
call hash_table_init(grid%edge_hash,size(grid%edge))
if (ierr == ALLOC_FAILED) then
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%face_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2, &
              grid%errest_eigenvalue,stat=dstat)
   return
endif
call hash_table_init(grid%vert_hash,size(grid%vertex))
if (ierr == ALLOC_FAILED) then
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%face_hash)
   call hash_table_destroy(grid%elem_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2, &
              grid%errest_eigenvalue,stat=dstat)
   return
endif
grid%next_free_elem = 1
grid%next_free_face = 1
grid%next_free_edge = 1
grid%next_free_vert = 1
allocate(grid%head_level_elem(8),stat=astat1)
if (astat1 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%head_level_elem,stat=dstat)
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%face_hash)
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%vert_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2, &
              grid%errest_eigenvalue,stat=dstat)
   return
endif
grid%head_level_elem = END_OF_LIST
grid%nelem = 0
grid%nelem_leaf = 0
grid%nelem_leaf_own = 0
grid%nface = 0
grid%nface_own = 0
grid%nedge = 0
grid%nedge_own = 0
grid%nvert = 0
grid%nvert_own = 0
grid%dof = 0
grid%dof_own = 0
grid%nlev = 0

! elements

allocate(grid%element_errind(size(grid%element),max(1,nev)),stat=astat1)
if (astat1 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%element_errind,stat=dstat)
   deallocate(grid%head_level_elem,stat=dstat)
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%face_hash)
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%vert_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2, &
              grid%errest_eigenvalue,stat=dstat)
   return
endif

! faces

nullify(grid%face_type)

! edges

allocate(grid%edge_type(size(grid%edge),grid%system_size),stat=astat1)
if (astat1 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%edge_type,stat=dstat)
   deallocate(grid%element_errind,stat=dstat)
   deallocate(grid%head_level_elem,stat=dstat)
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%face_hash)
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%vert_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2, &
              grid%errest_eigenvalue,stat=dstat)
   return
endif

! vertices

allocate(grid%vertex_type(size(grid%vertex),grid%system_size),stat=astat1)
if (astat1 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%edge_type,stat=dstat)
   deallocate(grid%vertex_type,stat=dstat)
   deallocate(grid%element_errind,stat=dstat)
   deallocate(grid%head_level_elem,stat=dstat)
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%face_hash)
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%vert_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2, &
              grid%errest_eigenvalue,stat=dstat)
   return
endif
allocate(grid%vertex_solution(size(grid%vertex),grid%system_size,nevdim), &
         stat=astat1)
if (astat1 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%edge_type,grid%vertex_type,stat=dstat)
   deallocate(grid%element_errind,stat=dstat)
   deallocate(grid%head_level_elem,stat=dstat)
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%face_hash)
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%vert_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2, &
              grid%errest_eigenvalue,stat=dstat)
   return
endif
nullify(grid%vertex_exact,grid%vertex_oldsoln)
grid%nsoln = grid%system_size*nevdim
grid%num_eval = nev
grid%have_true = .false.

! boundary parameter ranges, if using boundary subroutines

if (boundary_npiece(0) > 0) then
   npiece = boundary_npiece(0)
   i = 1
   do while (boundary_npiece(i) > 0)
      npiece = npiece + boundary_npiece(i)
      i = i+1
   end do
   allocate(grid%bp_start(npiece),grid%bp_finish(npiece),stat=astat1)
   if (astat1 /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in allocate_grid")
      stop
   endif
   call boundary_param(grid%bp_start,grid%bp_finish)
else
   nullify(grid%bp_start,grid%bp_finish)
endif

! eigenvalue monitor

if (type == EIGENVALUE) then
   allocate(grid%eigen_results%eigenvalue(nev), &
            grid%eigen_results%eigensolver_l2_resid(nev), &
            grid%eigen_results%eigensolver_errbound(nev), &
            stat=astat1)
   if (astat1 /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in allocate_grid")
      deallocate(grid%vertex_solution)
      deallocate(grid%edge_type,grid%vertex_type,stat=dstat)
      deallocate(grid%element_errind,stat=dstat)
      deallocate(grid%head_level_elem,stat=dstat)
      call hash_table_destroy(grid%elem_hash)
      call hash_table_destroy(grid%face_hash)
      call hash_table_destroy(grid%edge_hash)
      call hash_table_destroy(grid%vert_hash)
      deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
                 grid%errest_Linf,grid%errest_L2, &
                 grid%errest_eigenvalue,stat=dstat)
      if (associated(grid%eigen_results%eigenvalue)) &
         deallocate(grid%eigen_results%eigenvalue)
      if (associated(grid%eigen_results%eigensolver_l2_resid)) &
         deallocate(grid%eigen_results%eigensolver_l2_resid)
      if (associated(grid%eigen_results%eigensolver_errbound)) &
         deallocate(grid%eigen_results%eigensolver_errbound)
      return
   endif
   grid%eigen_results%eigenvalue = 0.0_my_real
   grid%eigen_results%eigensolver_l2_resid = 0.0_my_real
   grid%eigen_results%eigensolver_errbound = 0.0_my_real
else
   nullify(grid%eigen_results%eigenvalue, &
           grid%eigen_results%eigensolver_l2_resid, &
           grid%eigen_results%eigensolver_errbound)
endif

grid%errind_up2date = .false.
grid%oldsoln_exists = .false.

!call count_grid_memory(grid)

secret_grid => grid

end subroutine allocate_grid

!          ---------------
subroutine deallocate_grid(grid)
!          ---------------

!----------------------------------------------------
! This routine deallocates memory for a grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid

!----------------------------------------------------
! Local variables:

integer :: i, dstat

!----------------------------------------------------
! Begin executable code

if (associated(grid%vertex_type)) deallocate(grid%vertex_type,stat=dstat)
if (associated(grid%vertex_solution)) deallocate(grid%vertex_solution,stat=dstat)
if (associated(grid%vertex_exact)) deallocate(grid%vertex_exact,stat=dstat)
if (associated(grid%vertex_oldsoln)) deallocate(grid%vertex_oldsoln,stat=dstat)
if (associated(grid%edge_type)) deallocate(grid%edge_type,stat=dstat)
do i=1,grid%biggest_edge
   if (associated(grid%edge(i)%solution)) deallocate(grid%edge(i)%solution,stat=dstat)
   if (associated(grid%edge(i)%exact)) deallocate(grid%edge(i)%exact,stat=dstat)
   if (associated(grid%edge(i)%oldsoln)) deallocate(grid%edge(i)%oldsoln,stat=dstat)
end do
if (associated(grid%element_errind)) deallocate(grid%element_errind,stat=dstat)
do i=1,grid%biggest_elem
   if (associated(grid%element(i)%solution)) deallocate(grid%element(i)%solution,stat=dstat)
   if (associated(grid%element(i)%exact)) deallocate(grid%element(i)%exact,stat=dstat)
   if (associated(grid%element(i)%oldsoln)) deallocate(grid%element(i)%oldsoln,stat=dstat)
end do
deallocate(grid%element, grid%edge, grid%vertex, grid%head_level_elem, &
            grid%errest_energy, grid%errest_Linf, &
           grid%errest_L2, grid%errest_eigenvalue,stat=dstat)
if (associated(grid%eigen_results%eigenvalue)) then
   deallocate(grid%eigen_results%eigenvalue, &
              grid%eigen_results%eigensolver_l2_resid, &
              grid%eigen_results%eigensolver_errbound)
endif
if (associated(grid%initial_neighbor)) &
   deallocate(grid%initial_neighbor, stat=dstat)
call hash_table_destroy(grid%elem_hash)
call hash_table_destroy(grid%face_hash)
call hash_table_destroy(grid%edge_hash)
call hash_table_destroy(grid%vert_hash)
if (associated(grid%bp_start)) deallocate(grid%bp_start,grid%bp_finish,stat=dstat)

nullify(secret_grid)

end subroutine deallocate_grid

!          ----------------
subroutine realloc_solution(grid,degree,system_size,neval)
!          ----------------

!----------------------------------------------------
! This routine reallocates the solution components of the grid for
! the given degree, system size and number of eigenvalues, keeping the first
! nsolut old solutions and filling in 0.0 if nsolut is bigger than
! the old number of solutions
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: degree, system_size, neval
!----------------------------------------------------
! Local variables:

integer :: i, ncopy, ncopyss, ncopyev, astat, nsolut, neq_edge, neq_bubble, &
           neq_copy, deg, lev, elem, edge
real(my_real), pointer :: temp3(:,:,:)
!----------------------------------------------------
! Begin executable code

nsolut = system_size*neval
ncopy = min(nsolut,grid%nsoln)
ncopyss = min(size(grid%vertex_solution,2),system_size)
ncopyev = min(size(grid%vertex_solution,3),neval)

deallocate(grid%errest_energy, grid%errest_Linf, grid%errest_L2, &
           grid%errest_eigenvalue,stat=astat)
allocate(grid%errest_energy(neval), grid%errest_Linf(nsolut), &
         grid%errest_L2(nsolut), grid%errest_eigenvalue(neval),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in realloc_solution")
   return
endif

deallocate(grid%element_errind,stat=astat)
allocate(grid%element_errind(size(grid%element),neval),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in realloc_solution")
   return
endif

! vertices

if (size(grid%vertex_solution,2) /= system_size .or. &
    size(grid%vertex_solution,3) /= neval) then
   allocate(temp3(size(grid%vertex),system_size,neval),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in realloc_solution")
      return
   endif
   temp3(1:grid%biggest_vert,1:ncopyss,1:ncopyev) = grid%vertex_solution(1:grid%biggest_vert,1:ncopyss,1:ncopyev)
   temp3(1:grid%biggest_vert,ncopyss+1:system_size,1:ncopyev) = 0.0_my_real
   temp3(1:grid%biggest_vert,1:ncopyss,ncopyev+1:neval) = 0.0_my_real
   temp3(1:grid%biggest_vert,ncopyss+1:system_size,ncopyev+1:neval)=0.0_my_real
   deallocate(grid%vertex_solution,stat=astat)
   grid%vertex_solution => temp3
   nullify(temp3)
endif

if (grid%have_true) then
 if (size(grid%vertex_exact,2) /= system_size .or. &
     size(grid%vertex_exact,3) /= neval) then
   allocate(temp3(size(grid%vertex),system_size,neval),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in realloc_solution")
      return
   endif
   temp3(1:grid%biggest_vert,1:ncopyss,1:ncopyev) = &
      grid%vertex_exact(1:grid%biggest_vert,1:ncopyss,1:ncopyev)
   temp3(1:grid%biggest_vert,ncopyss+1:system_size,1:ncopyev) = 0.0_my_real
   temp3(1:grid%biggest_vert,1:ncopyss,ncopyev+1:neval) = 0.0_my_real
   temp3(1:grid%biggest_vert,ncopyss+1:system_size,ncopyev+1:neval)=0.0_my_real
   deallocate(grid%vertex_exact,stat=astat)
   grid%vertex_exact => temp3
   nullify(temp3)
 endif
endif

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)

! edges

      do i=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(i)
         if (degree > 0) then
            deg = degree
            if (grid%edge(edge)%degree > 0) then
               grid%dof = grid%dof - (grid%edge(edge)%degree-1) * &
                                   count(grid%edge_type(edge,:)/=PERIODIC_SLAVE)
               if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
                  grid%dof_own = grid%dof_own - &
                                 (grid%edge(edge)%degree-1) * &
                                 count(grid%edge_type(edge,:)/=PERIODIC_SLAVE)
               endif
            endif
            grid%edge(edge)%degree = degree
            grid%dof = grid%dof + (grid%edge(edge)%degree-1) * &
                                   count(grid%edge_type(edge,:)/=PERIODIC_SLAVE)
            if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
               grid%dof_own = grid%dof_own + &
                              (grid%edge(edge)%degree-1) * &
                              count(grid%edge_type(edge,:)/=PERIODIC_SLAVE)
            endif
         else
            deg = grid%edge(edge)%degree
         endif
         neq_edge = deg-1
         if (deg <= 1) then
            if (associated(grid%edge(edge)%solution)) then
               deallocate(grid%edge(edge)%solution)
            endif
            if (associated(grid%edge(edge)%exact)) then
               deallocate(grid%edge(edge)%exact)
            endif
            nullify(grid%edge(edge)%solution,grid%edge(edge)%exact)
         else
            if (associated(grid%edge(edge)%solution)) then
               if (size(grid%edge(edge)%solution,dim=1) == neq_edge .and. &
                   size(grid%edge(edge)%solution,dim=2) == system_size .and. &
                   size(grid%edge(edge)%solution,dim=3) == neval) cycle
            endif
            if (associated(grid%edge(edge)%solution)) then
               neq_copy = min(neq_edge,size(grid%edge(edge)%solution,dim=1))
            else
               neq_copy = 0
            endif
            allocate(temp3(neq_edge,system_size,neval),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("memory allocation failed in realloc_solution")
               return
            endif
            temp3 = 0.0_my_real
            if (neq_copy > 0) then
               temp3(1:neq_copy,1:ncopyss,1:ncopyev) = &
                  grid%edge(edge)%solution(1:neq_copy,1:ncopyss,1:ncopyev)
            endif
            if (associated(grid%edge(edge)%solution)) then
               deallocate(grid%edge(edge)%solution,stat=astat)
            endif
            grid%edge(edge)%solution => temp3
            nullify(temp3)
            if (grid%have_true) then
               if (associated(grid%edge(edge)%exact)) then
                  neq_copy = min(neq_edge,size(grid%edge(edge)%exact,dim=1))
               else
                  neq_copy = 0
               endif
               allocate(temp3(neq_edge,system_size,neval),stat=astat)
               if (astat /= 0) then
                  ierr = ALLOC_FAILED
                  call fatal("memory allocation failed in realloc_solution")
                  return
               endif
               temp3 = 0.0_my_real
               if (neq_copy > 0) then
                  temp3(1:neq_copy,1:ncopyss,1:ncopyev) = &
                     grid%edge(edge)%exact(1:neq_copy,1:ncopyss,1:ncopyev)
               endif
               if (associated(grid%edge(edge)%exact)) then
                  deallocate(grid%edge(edge)%exact,stat=astat)
               endif
               grid%edge(edge)%exact => temp3
               nullify(temp3)
            endif
         endif
      end do

! element

      if (degree > 0) then
         deg = degree
         if (grid%element(elem)%degree > 0) then
            grid%dof = grid%dof - system_size * &
                 element_dof(grid%element(elem)%degree)
            if (grid%element(elem)%iown) then
               grid%dof_own = grid%dof_own - system_size * &
                 element_dof(grid%element(elem)%degree)
            endif
         endif
         grid%element(elem)%degree = degree
         grid%dof = grid%dof + system_size * &
                 element_dof(grid%element(elem)%degree)
         if (grid%element(elem)%iown) then
            grid%dof_own = grid%dof_own + system_size * &
                 element_dof(grid%element(elem)%degree)
         endif
      else
         deg = grid%element(elem)%degree
      endif
      neq_bubble = element_dof(deg)
      if (deg <= 2) then
         if (associated(grid%element(elem)%solution)) then
            deallocate(grid%element(elem)%solution)
         endif
         if (associated(grid%element(elem)%exact)) then
            deallocate(grid%element(elem)%exact)
         endif
         nullify(grid%element(elem)%solution,grid%element(elem)%exact)
      else
         if (associated(grid%element(elem)%solution)) then
            if (size(grid%element(elem)%solution,dim=1) == neq_bubble .and. &
                size(grid%element(elem)%solution,dim=2) == system_size .and. &
                size(grid%element(elem)%solution,dim=3) == neval) then
               elem = grid%element(elem)%next
               cycle
            endif
         endif
         if (associated(grid%element(elem)%solution)) then
            neq_copy = min(neq_bubble,size(grid%element(elem)%solution,dim=1))
         else
            neq_copy = 0
         endif
         allocate(temp3(neq_bubble,system_size,neval),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in realloc_solution")
            return
         endif
         temp3 = 0.0_my_real
         if (neq_copy > 0) then
            temp3(1:neq_copy,1:ncopyss,1:ncopyev) = &
               grid%element(elem)%solution(1:neq_copy,1:ncopyss,1:ncopyev)
         endif
         if (associated(grid%element(elem)%solution)) then
            deallocate(grid%element(elem)%solution,stat=astat)
         endif
         grid%element(elem)%solution => temp3
         nullify(temp3)
         if (grid%have_true) then
            if (associated(grid%element(elem)%exact)) then
               neq_copy = min(neq_bubble,size(grid%element(elem)%exact,dim=1))
            else
               neq_copy = 0
            endif
            allocate(temp3(neq_bubble,system_size,neval),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("memory allocation failed in realloc_solution")
               return
            endif
            temp3 = 0.0_my_real
            if (neq_copy > 0) then
               temp3(1:neq_copy,1:ncopyss,1:ncopyev) = &
                  grid%element(elem)%exact(1:neq_copy,1:ncopyss,1:ncopyev)
            endif
            if (associated(grid%element(elem)%exact)) then
               deallocate(grid%element(elem)%exact,stat=astat)
            endif
            grid%element(elem)%exact => temp3
            nullify(temp3)
         endif
      endif
      elem = grid%element(elem)%next
   end do
end do

grid%nsoln = nsolut
grid%errind_up2date = .false.

end subroutine realloc_solution

!          ------------------
subroutine reset_dirich_exact(grid)
!          ------------------

!----------------------------------------------------
! This routine resets Dirichlet boundary conditions and the exact solution
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: lev, vert, edge, elem, i, j, ivert
integer :: bctype(grid%system_size)
real(my_real) :: bcrhs(grid%system_size), &
                 bccoef(grid%system_size,grid%system_size)
logical(small_logical) :: visited_vert(grid%biggest_vert), &
                          visited_edge(grid%biggest_edge)
!----------------------------------------------------
! Begin executable code

! for each vertex, set Dirichlet b.c. and exact

visited_vert = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do ivert=1,VERTICES_PER_ELEMENT
         vert = grid%element(elem)%vertex(ivert)
         if (visited_vert(vert)) cycle
         visited_vert(vert) = .true.
         if (any(grid%vertex_type(vert,:) == DIRICHLET) .or. &
             any(grid%vertex_type(vert,:) == PERIODIC_SLAVE_DIR) .or. &
             any(grid%vertex_type(vert,:) == PERIODIC_MASTER_DIR)) then
            call my_bconds(grid%vertex(vert)%coord%x,grid%vertex(vert)%coord%y,&
                           zcoord(grid%vertex(vert)%coord), &
                           grid%vertex(vert)%bmark,bctype,bccoef,bcrhs)
            do i=1,grid%system_size
              if (grid%vertex_type(vert,i) == DIRICHLET .or. &
                  grid%vertex_type(vert,i) == PERIODIC_SLAVE_DIR .or. &
                  grid%vertex_type(vert,i) == PERIODIC_MASTER_DIR) then
               grid%vertex_solution(vert,i,:) = bcrhs(i)
              endif
            end do
         endif
         if (grid%have_true) then
            do j=1,max(1,grid%num_eval)
               do i=1,grid%system_size
                  grid%vertex_exact(vert,i,j) = &
                     my_trues(grid%vertex(vert)%coord%x, &
                              grid%vertex(vert)%coord%y, &
                              zcoord(grid%vertex(vert)%coord),i,j)
               end do
            end do
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

! for each leaf element
!   for each edge of element that has not already been done
!      set Dirichlet b.c. and exact
!   set exact

visited_edge = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
         do j=1,EDGES_PER_ELEMENT
            edge = grid%element(elem)%edge(j)
            if (visited_edge(edge)) cycle
            visited_edge(edge) = .true.
            do i=1,grid%system_size
               if (grid%edge_type(edge,i) == DIRICHLET) then
                  call edge_exact(grid,edge,i,"d")
               endif
               if (grid%have_true) call edge_exact(grid,edge,i,"t")
            end do
         end do
         if (grid%have_true) then
            do i=1,grid%system_size
               call elem_exact(grid,elem,i,"t")
            end do
         endif
      endif
      elem = grid%element(elem)%next
   end do
end do

end subroutine reset_dirich_exact

!          ---------
subroutine copy_grid(old_grid,new_grid)
!          ---------

!----------------------------------------------------
! This routine makes an exact copy of old_grid in new_grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: old_grid
type(grid_type), intent(out) :: new_grid
!----------------------------------------------------
! Local variables:

integer :: i, astat
!----------------------------------------------------
! Begin executable code

! elements

if (associated(old_grid%element)) then
   allocate(new_grid%element(size(old_grid%element)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   do i=1,old_grid%biggest_elem
      new_grid%element(i)%gid = old_grid%element(i)%gid
      new_grid%element(i)%mate = old_grid%element(i)%mate
      new_grid%element(i)%weight = old_grid%element(i)%weight
      if (associated(old_grid%element(i)%solution)) then
         allocate(new_grid%element(i)%solution( &
                  size(old_grid%element(i)%solution,dim=1), &
                  size(old_grid%element(i)%solution,dim=2), &
                  size(old_grid%element(i)%solution,dim=3)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in copy_grid")
            return
         endif
         new_grid%element(i)%solution = old_grid%element(i)%solution
      else
         nullify(new_grid%element(i)%solution)
      endif
      if (associated(old_grid%element(i)%exact)) then
         allocate(new_grid%element(i)%exact( &
                  size(old_grid%element(i)%exact,dim=1), &
                  size(old_grid%element(i)%exact,dim=2), &
                  size(old_grid%element(i)%exact,dim=3)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in copy_grid")
            return
         endif
         new_grid%element(i)%exact = old_grid%element(i)%exact
      else
         nullify(new_grid%element(i)%exact)
      endif
      if (associated(old_grid%element(i)%oldsoln)) then
         allocate(new_grid%element(i)%oldsoln( &
                  size(old_grid%element(i)%oldsoln,dim=1), &
                  size(old_grid%element(i)%oldsoln,dim=2), &
                  size(old_grid%element(i)%oldsoln,dim=3)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in copy_grid")
            return
         endif
         new_grid%element(i)%oldsoln = old_grid%element(i)%oldsoln
      else
         nullify(new_grid%element(i)%oldsoln)
      endif
      new_grid%element(i)%work = old_grid%element(i)%work
      new_grid%element(i)%sp_eta_pred = old_grid%element(i)%sp_eta_pred
      new_grid%element(i)%vertex = old_grid%element(i)%vertex
      new_grid%element(i)%edge = old_grid%element(i)%edge
      new_grid%element(i)%degree = old_grid%element(i)%degree
      new_grid%element(i)%level = old_grid%element(i)%level
      new_grid%element(i)%in = old_grid%element(i)%in
      new_grid%element(i)%out = old_grid%element(i)%out
      new_grid%element(i)%order = old_grid%element(i)%order
      new_grid%element(i)%tags = old_grid%element(i)%tags
      new_grid%element(i)%next = old_grid%element(i)%next
      new_grid%element(i)%previous = old_grid%element(i)%previous
      new_grid%element(i)%isleaf = old_grid%element(i)%isleaf
      new_grid%element(i)%oldleaf = old_grid%element(i)%oldleaf
      new_grid%element(i)%iown = old_grid%element(i)%iown
      new_grid%element(i)%hrefined_unowned = old_grid%element(i)%hrefined_unowned
      new_grid%element(i)%prefined_unowned = old_grid%element(i)%prefined_unowned
   end do
else
   nullify(new_grid%element)
endif

! faces

allocate(new_grid%face(0),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in copy_grid")
   return
endif
nullify(new_grid%face_type)

! edges

if (associated(old_grid%edge)) then
   allocate(new_grid%edge(size(old_grid%edge)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   do i=1,old_grid%biggest_edge
      new_grid%edge(i)%gid = old_grid%edge(i)%gid
      new_grid%edge(i)%vertex = old_grid%edge(i)%vertex
      new_grid%edge(i)%bmark = old_grid%edge(i)%bmark
      new_grid%edge(i)%degree = old_grid%edge(i)%degree
      new_grid%edge(i)%assoc_elem = old_grid%edge(i)%assoc_elem
      new_grid%edge(i)%tags = old_grid%edge(i)%tags
      if (associated(old_grid%edge(i)%solution)) then
         allocate(new_grid%edge(i)%solution( &
                  size(old_grid%edge(i)%solution,dim=1), &
                  size(old_grid%edge(i)%solution,dim=2), &
                  size(old_grid%edge(i)%solution,dim=3)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in copy_grid")
            return
         endif
         new_grid%edge(i)%solution = old_grid%edge(i)%solution
      else
         nullify(new_grid%edge(i)%solution)
      endif
      if (associated(old_grid%edge(i)%exact)) then
         allocate(new_grid%edge(i)%exact( &
                  size(old_grid%edge(i)%exact,dim=1), &
                  size(old_grid%edge(i)%exact,dim=2), &
                  size(old_grid%edge(i)%exact,dim=3)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in copy_grid")
            return
         endif
         new_grid%edge(i)%exact = old_grid%edge(i)%exact
      else
         nullify(new_grid%edge(i)%exact)
      endif
      if (associated(old_grid%edge(i)%oldsoln)) then
         allocate(new_grid%edge(i)%oldsoln( &
                  size(old_grid%edge(i)%oldsoln,dim=1), &
                  size(old_grid%edge(i)%oldsoln,dim=2), &
                  size(old_grid%edge(i)%oldsoln,dim=3)),stat=astat)
         new_grid%edge(i)%oldsoln = old_grid%edge(i)%oldsoln
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in copy_grid")
            return
         endif
      else
         nullify(new_grid%edge(i)%oldsoln)
      endif
      new_grid%edge(i)%next = old_grid%edge(i)%next
   end do
else
   nullify(new_grid%edge)
endif

! vertices

if (associated(old_grid%vertex)) then
   allocate(new_grid%vertex(size(old_grid%vertex)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   do i=1,old_grid%biggest_vert
      new_grid%vertex(i)%gid = old_grid%vertex(i)%gid
      new_grid%vertex(i)%coord = old_grid%vertex(i)%coord
      new_grid%vertex(i)%bparam = old_grid%vertex(i)%bparam
      new_grid%vertex(i)%bmark = old_grid%vertex(i)%bmark
      new_grid%vertex(i)%assoc_elem = old_grid%vertex(i)%assoc_elem
      new_grid%vertex(i)%tags = old_grid%vertex(i)%tags
      new_grid%vertex(i)%next = old_grid%vertex(i)%next
   end do
else
   nullify(new_grid%vertex)
endif

! remaining fields

call hash_table_copy(old_grid%elem_hash,new_grid%elem_hash)
call hash_table_copy(old_grid%face_hash,new_grid%face_hash)
call hash_table_copy(old_grid%edge_hash,new_grid%edge_hash)
call hash_table_copy(old_grid%vert_hash,new_grid%vert_hash)
new_grid%boundbox_min = old_grid%boundbox_min
new_grid%boundbox_max = old_grid%boundbox_max
if (associated(old_grid%vertex_solution)) then
   allocate(new_grid%vertex_solution( &
            size(old_grid%vertex_solution,dim=1), &
            size(old_grid%vertex_solution,dim=2), &
            size(old_grid%vertex_solution,dim=3)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%vertex_solution = old_grid%vertex_solution
else
   nullify(new_grid%vertex_solution)
endif
if (associated(old_grid%vertex_exact)) then
   allocate(new_grid%vertex_exact( &
            size(old_grid%vertex_exact,dim=1), &
            size(old_grid%vertex_exact,dim=2), &
            size(old_grid%vertex_exact,dim=3)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%vertex_exact = old_grid%vertex_exact
else
   nullify(new_grid%vertex_exact)
endif
if (associated(old_grid%vertex_oldsoln)) then
   allocate(new_grid%vertex_oldsoln( &
            size(old_grid%vertex_oldsoln,dim=1), &
            size(old_grid%vertex_oldsoln,dim=2), &
            size(old_grid%vertex_oldsoln,dim=3)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%vertex_oldsoln = old_grid%vertex_oldsoln
else
   nullify(new_grid%vertex_oldsoln)
endif
if (associated(old_grid%element_errind)) then
   allocate(new_grid%element_errind( &
            size(old_grid%element_errind,dim=1), &
            size(old_grid%element_errind,dim=2)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%element_errind = old_grid%element_errind
else
   nullify(new_grid%element_errind)
endif
if (associated(old_grid%errest_energy)) then
   allocate(new_grid%errest_energy(size(old_grid%errest_energy)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%errest_energy = old_grid%errest_energy
else
   nullify(new_grid%errest_energy)
endif
if (associated(old_grid%errest_Linf)) then
   allocate(new_grid%errest_Linf(size(old_grid%errest_Linf)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%errest_Linf = old_grid%errest_Linf
else
   nullify(new_grid%errest_Linf)
endif
if (associated(old_grid%errest_L2)) then
   allocate(new_grid%errest_L2(size(old_grid%errest_L2)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%errest_L2 = old_grid%errest_L2
else
   nullify(new_grid%errest_L2)
endif
if (associated(old_grid%errest_eigenvalue)) then
   allocate(new_grid%errest_eigenvalue(size(old_grid%errest_eigenvalue)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%errest_eigenvalue = old_grid%errest_eigenvalue
else
   nullify(new_grid%errest_eigenvalue)
endif
new_grid%max_blen = old_grid%max_blen
if (associated(old_grid%bp_start)) then
   allocate(new_grid%bp_start(size(old_grid%bp_start)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%bp_start = old_grid%bp_start
else
   nullify(new_grid%bp_start)
endif
if (associated(old_grid%bp_finish)) then
   allocate(new_grid%bp_finish(size(old_grid%bp_finish)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%bp_finish = old_grid%bp_finish
else
   nullify(new_grid%bp_finish)
endif
if (associated(old_grid%edge_type)) then
   allocate(new_grid%edge_type( &
            size(old_grid%edge_type,dim=1), &
            size(old_grid%edge_type,dim=2)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%edge_type = old_grid%edge_type
else
   nullify(new_grid%edge_type)
endif
if (associated(old_grid%vertex_type)) then
   allocate(new_grid%vertex_type( &
            size(old_grid%vertex_type,dim=1), &
            size(old_grid%vertex_type,dim=2)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%vertex_type = old_grid%vertex_type
else
   nullify(new_grid%vertex_type)
endif
if (associated(old_grid%initial_neighbor)) then
   allocate(new_grid%initial_neighbor( &
            size(old_grid%initial_neighbor,dim=1), &
            size(old_grid%initial_neighbor,dim=2)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%initial_neighbor = old_grid%initial_neighbor
else
   nullify(new_grid%initial_neighbor)
endif
if (associated(old_grid%head_level_elem)) then
   allocate(new_grid%head_level_elem(size(old_grid%head_level_elem)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%head_level_elem = old_grid%head_level_elem
else
   nullify(new_grid%head_level_elem)
endif
new_grid%next_free_elem = old_grid%next_free_elem
new_grid%next_free_face = old_grid%next_free_face
new_grid%next_free_edge = old_grid%next_free_edge
new_grid%next_free_vert = old_grid%next_free_vert
new_grid%partition = old_grid%partition
new_grid%system_size = old_grid%system_size
new_grid%num_eval = old_grid%num_eval
new_grid%nsoln = old_grid%nsoln
new_grid%nelem = old_grid%nelem
new_grid%nelem_leaf = old_grid%nelem_leaf
new_grid%nelem_leaf_own = old_grid%nelem_leaf_own
new_grid%nface = old_grid%nface
new_grid%nface_own = old_grid%nface_own
new_grid%nedge = old_grid%nedge
new_grid%nedge_own = old_grid%nedge_own
new_grid%nvert = old_grid%nvert
new_grid%nvert_own = old_grid%nvert_own
new_grid%nlev = old_grid%nlev
new_grid%dof = old_grid%dof
new_grid%dof_own = old_grid%dof_own
new_grid%biggest_elem = old_grid%biggest_elem
new_grid%biggest_face = old_grid%biggest_face
new_grid%biggest_edge = old_grid%biggest_edge
new_grid%biggest_vert = old_grid%biggest_vert
new_grid%errtype = old_grid%errtype
new_grid%errind_up2date = old_grid%errind_up2date
new_grid%oldsoln_exists = old_grid%oldsoln_exists
new_grid%have_true = old_grid%have_true
new_grid%triangle_files = old_grid%triangle_files
if (associated(old_grid%eigen_results%eigenvalue)) then
   allocate(new_grid%eigen_results%eigenvalue(size(old_grid%eigen_results%eigenvalue)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%eigen_results%eigenvalue = old_grid%eigen_results%eigenvalue
else
   nullify(new_grid%eigen_results%eigenvalue)
endif
if (associated(old_grid%eigen_results%eigensolver_l2_resid)) then
   allocate(new_grid%eigen_results%eigensolver_l2_resid(size(old_grid%eigen_results%eigensolver_l2_resid)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%eigen_results%eigensolver_l2_resid = old_grid%eigen_results%eigensolver_l2_resid
else
   nullify(new_grid%eigen_results%eigensolver_l2_resid)
endif
if (associated(old_grid%eigen_results%eigensolver_errbound)) then
   allocate(new_grid%eigen_results%eigensolver_errbound(size(old_grid%eigen_results%eigensolver_errbound)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_grid")
      return
   endif
   new_grid%eigen_results%eigensolver_errbound = old_grid%eigen_results%eigensolver_errbound
else
   nullify(new_grid%eigen_results%eigensolver_errbound)
endif
new_grid%eigen_results%ncv = old_grid%eigen_results%ncv
new_grid%eigen_results%maxit = old_grid%eigen_results%maxit
new_grid%eigen_results%niter = old_grid%eigen_results%niter
new_grid%eigen_results%nconv = old_grid%eigen_results%nconv

end subroutine copy_grid

!          --------
subroutine copy_old(grid)
!          --------

!----------------------------------------------------
! This routine copies solution into oldsoln
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: lev, elem, astat, i, edge
!----------------------------------------------------
! Begin executable code

do lev=1,grid%nlev

! for each element

   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
         if (associated(grid%element(elem)%solution)) then

! if the element is a leaf and has a solution, then allocate oldsoln to the
! right size and copy solution to it

            if (associated(grid%element(elem)%oldsoln)) then
               if (size(grid%element(elem)%oldsoln,dim=1) /= &
                   size(grid%element(elem)%solution,dim=1) .or. &
                   size(grid%element(elem)%oldsoln,dim=2) /= &
                   size(grid%element(elem)%solution,dim=2) .or. &
                   size(grid%element(elem)%oldsoln,dim=3) /= &
                   size(grid%element(elem)%solution,dim=3)) then
                  deallocate(grid%element(elem)%oldsoln,stat=astat)
               endif
            endif
            if (.not. associated(grid%element(elem)%oldsoln)) then
               allocate(grid%element(elem)%oldsoln( &
                          size(grid%element(elem)%solution,dim=1), &
                          size(grid%element(elem)%solution,dim=2), &
                          size(grid%element(elem)%solution,dim=3)),stat=astat)
               if (astat /= 0) then
                  ierr = ALLOC_FAILED
                  call fatal("memory allocation failed in copy_old")
                  stop
               endif
            endif
            grid%element(elem)%oldsoln = grid%element(elem)%solution
            grid%element(elem)%oldleaf = .true.

         else

! if element is a leaf and does not have a solution, make sure it does not
! have oldsoln either

            if (associated(grid%element(elem)%oldsoln)) then
               deallocate(grid%element(elem)%oldsoln,stat=astat)
            endif
            grid%element(elem)%oldleaf = .true.

         endif
      else

! if element is not a leaf, make sure it does not have an oldsoln or say
! that it is an oldleaf

         if (associated(grid%element(elem)%oldsoln)) then
            deallocate(grid%element(elem)%oldsoln,stat=astat)
         endif
         grid%element(elem)%oldleaf = .false.

      endif

! for each edge

      do i=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(i)

         if (grid%element(elem)%isleaf) then
            if (associated(grid%edge(edge)%solution)) then

! if the element is a leaf and the edge solution is allocated, make sure the
! edge oldsoln is the same size and copy solution to it

               if (associated(grid%edge(edge)%oldsoln)) then
                  if (size(grid%edge(edge)%oldsoln,dim=1) /= &
                      size(grid%edge(edge)%solution,dim=1) .or. &
                      size(grid%edge(edge)%oldsoln,dim=2) /= &
                      size(grid%edge(edge)%solution,dim=2) .or. &
                      size(grid%edge(edge)%oldsoln,dim=3) /= &
                      size(grid%edge(edge)%solution,dim=3)) then
                     deallocate(grid%edge(edge)%oldsoln,stat=astat)
                  endif
               endif
               if (.not. associated(grid%edge(edge)%oldsoln)) then
                  allocate(grid%edge(edge)%oldsoln( &
                             size(grid%edge(edge)%solution,dim=1), &
                             size(grid%edge(edge)%solution,dim=2), &
                             size(grid%edge(edge)%solution,dim=3)),stat=astat)
                  if (astat /= 0) then
                     ierr = ALLOC_FAILED
                     call fatal("memory allocation failed in copy_old")
                     stop
                  endif
               endif
               grid%edge(edge)%oldsoln = grid%edge(edge)%solution

            else

! if element is a leaf and edge does not have a solution, make sure it does not
! have oldsoln either

               if (associated(grid%edge(edge)%oldsoln)) then
                  deallocate(grid%edge(edge)%oldsoln,stat=astat)
               endif
            endif

! if element is not a leaf, don't change the edge because the neighbor might
! be a leaf and want the edge set, and it won't cause any damage to leave it

         endif
      end do ! next edge

      elem = grid%element(elem)%next
   end do ! next element

end do ! next level

! make sure the oldsoln is the same size as solution and copy

if (associated(grid%vertex_oldsoln)) then
   if (size(grid%vertex_oldsoln,dim=1) /= &
       size(grid%vertex_solution,dim=1) .or. &
       size(grid%vertex_oldsoln,dim=2) /= &
       size(grid%vertex_solution,dim=2) .or. &
       size(grid%vertex_oldsoln,dim=3) /= &
       size(grid%vertex_solution,dim=3)) then
      deallocate(grid%vertex_oldsoln,stat=astat)
   endif
endif
if (.not. associated(grid%vertex_oldsoln)) then
   allocate(grid%vertex_oldsoln(size(grid%vertex_solution,1), &
            size(grid%vertex_solution,2),size(grid%vertex_solution,3)), &
            stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_old")
      stop
   endif
endif
grid%vertex_oldsoln = grid%vertex_solution

grid%oldsoln_exists = .true.

end subroutine copy_old

!          ----------
subroutine edge_exact(grid,edge,sysrank,what)
!          ----------

!----------------------------------------------------
! This routine sets "exact" solutions on an edge.  what indicates what is
! being set, and can be "d" for Dirichlet boundary conditions, "i" for
! initial conditions, or "t" for the true solution.
! It assumes the coefficients for the linear bases at the endpoints of
! edge are already set, and computes the coefficients of the high order
! bases along that edge as a least squares fit to the function.
! For multiple eigenvalue problems, all eigenfunctions are set.  It is
! assumed that they all satisfy the same boundary conditions (bconds is
! defined for component=1..system_size) but have different initial conditions
! and true solutions (iconds and trues are defined for component=1..system_size
! and eigen=1..num_eval)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: edge, sysrank
character(len=1), intent(in) :: what
!----------------------------------------------------
! Local variables:

integer :: i, j, nqp, nev, jerr, deg(4), astat, n, ss
real(my_real) :: xv(VERTICES_PER_ELEMENT),yv(VERTICES_PER_ELEMENT)
real(my_real), pointer :: weight(:), xq(:), yq(:), zq(:)
real(my_real), allocatable :: basis(:,:), a(:), b(:,:), true(:,:)
integer, allocatable :: ipiv(:)
real(my_real) :: c(grid%system_size,grid%system_size),rs(grid%system_size)
integer :: itype(grid%system_size)

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! Begin executable code

! nothing to do if only linear bases on this edge

if (grid%edge(edge)%degree < 2) return

! Dirichlet boundary conditions override initial conditions

if (what == "i" .and. grid%edge_type(edge,sysrank) == DIRICHLET) return

! useful constants

nev = max(1,grid%num_eval)
ss = grid%system_size

! local copy of the vertex coordinates, plus a fake third vertex to make
! a triangle

xv(1) = grid%vertex(grid%edge(edge)%vertex(1))%coord%x
xv(2) = grid%vertex(grid%edge(edge)%vertex(2))%coord%x
xv(3) = (xv(1)+xv(2))/2
yv(1) = grid%vertex(grid%edge(edge)%vertex(1))%coord%y
yv(2) = grid%vertex(grid%edge(edge)%vertex(2))%coord%y
yv(3) = (yv(1)+yv(2))/2
if (abs(xv(2)-xv(1)) > abs(yv(2)-yv(1))) then
   yv(3) = yv(3) + 1
else
   xv(3) = xv(3) + 1
endif

! get the quadrature rule for the edge

call quadrature_rule_line(min(MAX_QUAD_ORDER_LINE,grid%edge(edge)%degree+1), &
                          xv(1:2), yv(1:2), nqp, weight, xq, yq, jerr)
allocate(zq(size(xq))); zq=0 ! TEMP3D
if (jerr /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Error getting quadrature rule in edge_exact",intlist=(/jerr/))
   stop
endif

! evaluate basis functions at the quadrature points

allocate(basis(3+grid%edge(edge)%degree-1,nqp),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in edge_exact")
   stop
endif

deg = 1
deg(3) = grid%edge(edge)%degree

call p_hier_basis_func(xq,yq,xv,yv,deg,"a",basis)

! evaluate the function at the quadrature points

allocate(true(nqp,nev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in edge_exact")
   stop
endif

select case (what)
case ("d")
   do i=1,nqp
      call my_bconds(xq(i),yq(i),zq(i),grid%edge(edge)%bmark,itype,c,rs)
      true(i,:) = rs(sysrank)
   end do
case ("i")
   do j=1,nev
      do i=1,nqp
         true(i,j) = my_iconds(xq(i),yq(i),zq(i),sysrank,j)
      end do
   end do
case ("t")
   do j=1,nev
      do i=1,nqp
         true(i,j) = my_trues(xq(i),yq(i),zq(i),sysrank,j)
         if (true(i,j) == huge(0.0_my_real)) true(i,j) = 0.0_my_real
      end do
   end do
end select

! take off the linear bases part of the function

do j=1,nev
 if (what == "t") then
  if (associated(grid%vertex_exact)) then
   true(:,j) = true(:,j) - &
    grid%vertex_exact(grid%edge(edge)%vertex(1),sysrank,j)*basis(1,:)-&
    grid%vertex_exact(grid%edge(edge)%vertex(2),sysrank,j)*basis(2,:)
  endif
 else
  true(:,j) = true(:,j) - &
   grid%vertex_solution(grid%edge(edge)%vertex(1),sysrank,j)*basis(1,:) - &
   grid%vertex_solution(grid%edge(edge)%vertex(2),sysrank,j)*basis(2,:)
 endif
end do

! set up least squares linear system, storing the matrix in LAPACK symmetric
! packed form

n = grid%edge(edge)%degree - 1

allocate(a((n*(n+1))/2),b(n,nev),ipiv(n),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in edge_exact")
   stop
endif

do i=1,n
   do j=i,n
      a(i+(j*(j-1))/2) = sum(weight*basis(3+i,:)*basis(3+j,:))
   end do
   do j=1,nev
      b(i,j) = sum(weight*basis(3+i,:)*true(:,j))
   end do
end do

! solve the least squares system

if (my_real == kind(1.0e0)) then
   call sspsv("U",n,nev,a,ipiv,b,n,jerr)
elseif (my_real == kind(1.0d0)) then
   call dspsv("U",n,nev,a,ipiv,b,n,jerr)
else
   ierr = PHAML_INTERNAL_ERROR
      call fatal("my_real is neither single nor double precision. Can't call LAPACK")
   stop 
endif

! copy the solution to the edge data structure

do j=1,nev
   if (what == "t") then
      if (grid%have_true) then
         do i=1,n
            grid%edge(edge)%exact(i,sysrank,j) = b(i,j)
         end do
      endif
   else
      do i=1,n
         grid%edge(edge)%solution(i,sysrank,j) = b(i,j)
      end do
   endif
end do

deallocate(weight,xq,yq,zq,basis,a,b,true,ipiv)

end subroutine edge_exact

!          ----------
subroutine face_exact(grid,edge,sysrank,what)
!          ----------

!----------------------------------------------------
! In 3D this routine sets "exact" solutions for the face basis functions.
! In 2D it shouldn't be called.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: edge, sysrank
character(len=1), intent(in) :: what
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call warning("face_exact was called in 2D")

end subroutine face_exact

!          ----------
subroutine elem_exact(grid,elem,sysrank,what)
!          ----------

!----------------------------------------------------
! This routine sets "exact" solutions for bubble basis functions.
! what indicates what is being set, and can be "i" for initial conditions
! or "t" for the true solution.
! It assumes the coefficients for the bases at the vertices and edges
! of elem are already set, and computes the
! coefficients of the high order bubble bases in elem as a least squares
! fit to the function.
! For multiple eigenvalue problems, all eigenfunctions are set.  It is
! assumed that they have different initial conditions
! and true solutions (iconds and trues are defined for component=1..system_size
! and eigen=1..num_eval)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem, sysrank
character(len=1), intent(in) :: what
!----------------------------------------------------
! Local variables:

integer :: i, j, k, nqp, nev, ss, jerr, astat, n, degree, nbasis, isub, &
           deg(EDGES_PER_ELEMENT+1)
real(my_real) :: xv(VERTICES_PER_ELEMENT),yv(VERTICES_PER_ELEMENT)
real(my_real), pointer :: weight(:), xq(:), yq(:), zq(:)
real(my_real), allocatable :: basis(:,:), wbasis(:,:), a(:), b(:,:), true(:,:),&
                              afull(:,:)
integer, allocatable :: ipiv(:)

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! Begin executable code

! must be at least cubics to have bubble bases

if (grid%element(elem)%degree < 3) return

! useful constants

nev = max(1,grid%num_eval)
ss = grid%system_size

! local copy of the vertex coordinates

xv = grid%vertex(grid%element(elem)%vertex)%coord%x
yv = grid%vertex(grid%element(elem)%vertex)%coord%y

! get the quadrature rule for the triangle

degree = grid%element(elem)%degree
do i=1,EDGES_PER_ELEMENT
   degree = max(degree,grid%edge(grid%element(elem)%edge(i))%degree)
end do
if (degree > 1) then
   degree = min(MAX_QUAD_ORDER_TRI,2*degree-2)
endif
call quadrature_rule_tri(degree,xv,yv,nqp,weight,xq,yq,jerr,stay_in=.true.)
allocate(zq(size(xq))); zq=0 ! TEMP3D
if (jerr /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Error getting quadrature rule in elem_exact",intlist=(/jerr/))
   stop
endif

! evaluate basis functions at the quadrature points

nbasis = VERTICES_PER_ELEMENT
do i=1,EDGES_PER_ELEMENT 
   deg(i) = grid%edge(grid%element(elem)%edge(i))%degree
   nbasis = nbasis + max(0,deg(i)-1)
end do
deg(EDGES_PER_ELEMENT+1) = grid%element(elem)%degree
nbasis = nbasis + element_dof(deg(EDGES_PER_ELEMENT+1))
allocate(basis(nbasis,nqp),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in elem_exact")
   stop
endif

call p_hier_basis_func(xq,yq,xv,yv,deg,"a",basis)

! evaluate the function at the quadrature points

allocate(true(nqp,nev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in elem_exact")
   stop
endif

do j=1,nev
   do i=1,nqp
      if (what == "t") then
         true(i,j) = my_trues(xq(i),yq(i),zq(i),sysrank,j)
         if (true(i,j) == huge(0.0_my_real)) true(i,j) = 0.0_my_real
      else
         true(i,j) = my_iconds(xq(i),yq(i),zq(i),sysrank,j)
      endif
   end do
end do

! take off the vertex- and edge-basis parts of the function

do i=1,nev
 if (what == "t") then
  if (associated(grid%vertex_exact)) then
   true(:,i) = true(:,i) - &
        grid%vertex_exact(grid%element(elem)%vertex(1),sysrank,i)*basis(1,:) - &
        grid%vertex_exact(grid%element(elem)%vertex(2),sysrank,i)*basis(2,:) - &
        grid%vertex_exact(grid%element(elem)%vertex(3),sysrank,i)*basis(3,:)
   isub = 3
   do j=1,EDGES_PER_ELEMENT
      do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
         isub = isub + 1
         if (grid%have_true) then
            true(:,i) = true(:,i) - &
               grid%edge(grid%element(elem)%edge(j))%exact(k,sysrank,i)*basis(isub,:)
         endif
      end do
   end do
  endif
 else
   true(:,i) = true(:,i) - &
      grid%vertex_solution(grid%element(elem)%vertex(1),sysrank,i)*basis(1,:) -&
      grid%vertex_solution(grid%element(elem)%vertex(2),sysrank,i)*basis(2,:) -&
      grid%vertex_solution(grid%element(elem)%vertex(3),sysrank,i)*basis(3,:)
   isub = 3
   do j=1,EDGES_PER_ELEMENT
      do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
        isub = isub + 1
        true(:,i) = true(:,i) - &
         grid%edge(grid%element(elem)%edge(j))%solution(k,sysrank,i)*basis(isub,:)
      end do
   end do
 endif
end do

! set up least squares linear system, storing the matrix in LAPACK symmetric
! packed form

n = element_dof(grid%element(elem)%degree)

allocate(a((n*(n+1))/2),b(n,nev),ipiv(n),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in elem_exact")
   stop
endif

isub = 3
do i=1,EDGES_PER_ELEMENT 
   isub = isub + max(0,grid%edge(grid%element(elem)%edge(i))%degree-1)
end do

allocate(afull(n,n),wbasis(nqp,n),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in elem_exact")
   stop
endif
afull = 0.0_my_real
basis(1:n,:) = basis(isub+1:nbasis,:)
do i=1,nqp
   wbasis(i,1:n) = weight(i)*basis(1:n,i)
end do

if (my_real == kind(1.0)) then
   call sgemm("N","N",n,n,nqp,1.0,basis,nbasis,wbasis,nqp, &
              0.0,afull,n)
elseif (my_real == kind(1.0d0)) then
   call dgemm("N","N",n,n,nqp,1.0d0,basis,nbasis,wbasis,nqp, &
              0.0d0,afull,n)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither single nor double precision. Can't call GEMM")
   stop 
endif

do i=1,n
   do j=i,n
      a(i+(j*(j-1))/2) = afull(i,j)
   end do
end do

if (my_real == kind(1.0)) then
   call sgemm("T","N",n,nev,nqp,1.0,wbasis,nqp,true,nqp,0.0,b,n)
elseif (my_real == kind(1.0d0)) then
   call dgemm("T","N",n,nev,nqp,1.0d0,wbasis,nqp,true,nqp,0.0d0,b,n)
endif

! solve the least squares system

if (my_real == kind(1.0e0)) then
   call sspsv("U",n,nev,a,ipiv,b,n,jerr)
elseif (my_real == kind(1.0d0)) then
   call dspsv("U",n,nev,a,ipiv,b,n,jerr)
else
   ierr = PHAML_INTERNAL_ERROR
      call fatal("my_real is neither single nor double precision. Can't call LAPACK")
   stop
endif

! copy the solution to the element data structure

do j=1,nev
   if (what == "t") then
      if (grid%have_true) then
         do i=1,n
            grid%element(elem)%exact(i,sysrank,j) = b(i,j)
         end do
      endif
   else
      do i=1,n
         grid%element(elem)%solution(i,sysrank,j) = b(i,j)
      end do
   endif
end do

deallocate(weight,xq,yq,zq,basis,wbasis,a,b,true,ipiv,afull)

end subroutine elem_exact

!          -------------
subroutine point_on_edge(grid,x1,y1,bmark1,bparam1,x2,y2,bmark2,bparam2,x3,y3, &
                         emark,f,x,y,bparam)
!          -------------

!----------------------------------------------------
! This routine determines the (x,y) coordinate and bparam for a point on the
! (possibly curved) edge [(x1,y1),(x2,y2)] of the triangle with that edge and
! vertex (x3,y3).  The bmarks and bparams of the endpoints, and edge bparam
! are also input.  The point is 0 < f < 1 of the way from point 1 to 2.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: x1,y1,bparam1,x2,y2,bparam2,x3,y3
integer, intent(in) :: bmark1,bmark2,emark
real(my_real), intent(in) :: f
real(my_real), intent(out) :: x,y,bparam
!----------------------------------------------------
! Local variables:

integer :: piece
real(my_real) :: p1, p2, pm, temp, a1, a2, a3, xm, ym
!----------------------------------------------------
! Begin executable code

! if the boundary subroutines are not used or it is an interior edge,
! then return the midpoint and bparam doesn't matter

if (boundary_npiece(0) <= 0 .or. emark == 0) then
   x = (1-f)*x1 + f*x2
   y = (1-f)*y1 + f*y2
   bparam = 0
   return
endif

! determine parameters for the endpoints on piece emark

! first point is on same piece as edge
if (bmark1 == emark) then
   p1 = bparam1
! first point is a singleton.  A singleton cannot be the last piece of a
! hole, so we can look for the next larger index that is not a singleton.
! This will be emark iff bmark1 comes before emark.
elseif (bparam1 == grid%bp_start(bmark1) .and. &
        bparam1 == grid%bp_finish(bmark1)) then
   piece = bmark1+1
   do while (grid%bp_start(piece) == grid%bp_finish(piece))
      piece = piece + 1
   end do
   if (piece == emark) then
      p1 = grid%bp_start(emark)
   else
      p1 = grid%bp_finish(emark)
   endif
! first point is the beginning of it's piece, so also end of edge piece
elseif (bparam1 == grid%bp_start(bmark1)) then
   p1 = grid%bp_finish(emark)
! first point is the end of it's piece, so also the begining of edge piece
elseif (bparam1 == grid%bp_finish(bmark1)) then
   p1 = grid%bp_start(emark)
! failed to find what it is
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Failed to find edge endpoint in edge piece.")
endif

if (bmark2 == emark) then
   p2 = bparam2
elseif (bparam2 == grid%bp_start(bmark2) .and. &
        bparam2 == grid%bp_finish(bmark2)) then
   piece = bmark2+1
   do while (grid%bp_start(piece) == grid%bp_finish(piece))
      piece = piece + 1
   end do
   if (piece == emark) then
      p2 = grid%bp_start(emark)
   else
      p2 = grid%bp_finish(emark)
   endif
elseif (bparam2 == grid%bp_start(bmark2)) then
   p2 = grid%bp_finish(emark)
elseif (bparam2 == grid%bp_finish(bmark2)) then
   p2 = grid%bp_start(emark)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Failed to find edge endpoint in edge piece.")
endif

! find the equation of the line that goes through the desired point of the edge
! and the opposite vertex, in the form a1*X + a2*Y + a3 = 0

xm = (1-f)*x1 + f*x2
ym = (1-f)*y1 + f*y2
a1 = y3-ym
a2 = xm-x3
a3 = x3*ym - y3*xm

! if vertex 1 in the line equation is positive, swap p1 and p2 so that
! p1 goes with the negative one

if (a1*x1+a2*y1+a3 > 0.0_my_real) then
   temp = p1
   p1 = p2
   p2 = temp
endif

! use a bisection root finder to find the point where the line intersects
! the boundary

do
   pm = (p1+p2)/2
   call boundary_point(emark,pm,xm,ym)
   if (abs(a1*xm+a2*ym+a3) < 100*epsilon(0.0_my_real)) then
      x = xm
      y = ym
      bparam = pm
      return
   elseif (a1*xm+a2*ym+a3 < 0.0_my_real) then
      p1 = pm
   else
      p2 = pm
   endif
end do

end subroutine point_on_edge

!        -----------
function level2_mate(grid,parent,k)
!        -----------

!----------------------------------------------------
! This routine determines the mate of a triangle on level 2, which
! should be a child of the k'th neighbor of the level 1 parent
! RESTRICTION triangles
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: parent,k
type(hash_key) :: level2_mate
!----------------------------------------------------
! Local variables:

integer :: neigh, i
!----------------------------------------------------
! Begin executable code

neigh = grid%initial_neighbor(k,parent)

! easy case -- boundary

if (neigh == BOUNDARY) then
   level2_mate = BOUNDARY
   return
endif

! determine which neighbor parent is of the neighbor

do i=1,NEIGHBORS_PER_ELEMENT
   if (grid%initial_neighbor(i,neigh) == parent) exit
end do
if (i > NEIGHBORS_PER_ELEMENT) then
   call fatal("initial_neighbor is not reflexive")
   stop
endif
if (i == NEIGHBORS_PER_ELEMENT) then
   call fatal("initial_neighbor mates are not reflexive")
   stop
endif

! compute the mate

level2_mate = 2*grid%element(neigh)%gid + (2-i)

end function level2_mate

!        ----------------------
function matching_periodic_vert(grid,vert)
!        ----------------------

!----------------------------------------------------
! This routine returns the vertex that matches up with periodic boundary
! vertex vert.
! In the case of doubly periodic boundary vertices, if vert is slave then the
! master is returned, and if vert is the master, the first slave is returned.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: vert
integer :: matching_periodic_vert
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

matching_periodic_vert = -1

! If vert is a periodic slave, then one or more nexts give the master.

if (is_periodic_vert_slave(vert,grid)) then
   matching_periodic_vert = grid%vertex(vert)%next
   do while (is_periodic_vert_slave(matching_periodic_vert,grid))
      matching_periodic_vert = grid%vertex(matching_periodic_vert)%next
   end do

! Otherwise, next gives a matching vertex

elseif (is_periodic_vert_master(vert,grid)) then
   matching_periodic_vert = grid%vertex(vert)%next

else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Non-periodic vertex passed into matching_periodic_vert")
   stop

endif

end function matching_periodic_vert

!        ----------------------
function matching_periodic_edge(grid,edge)
!        ----------------------

!----------------------------------------------------
! This routine returns the edge that matches up with periodic boundary
! edge edge.  Element elem contains either edge or the matching edge.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: edge
integer :: matching_periodic_edge
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

matching_periodic_edge = -1

! If edge is periodic then next gives the matching edge, for both a master
! and a slave.

if (is_periodic_edge(edge,grid)) then
   matching_periodic_edge = grid%edge(edge)%next

else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Non-periodic edge passed into matching_periodic_edge")
   stop

endif

end function matching_periodic_edge

!        ----------------
function is_periodic_vert(vert,grid,comp)
!        ----------------

!----------------------------------------------------
! This routine returns true if the given vertex is periodic.
! If comp is present, it only checks component comp; otherwise it is periodic
! if any component is periodic.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
integer, intent(in) :: vert
type(grid_type), intent(in) :: grid
integer, intent(in), optional :: comp
logical :: is_periodic_vert
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (.not. grid%any_periodic) then
   is_periodic_vert = .false.
   return
endif

if (present(comp)) then
   is_periodic_vert=grid%vertex_type(vert,comp)==PERIODIC_MASTER     .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_MASTER_DIR .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_MASTER_NAT .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_MASTER_MIX .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_SLAVE      .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_SLAVE_DIR  .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_SLAVE_NAT  .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_SLAVE_MIX
else
   is_periodic_vert=any(grid%vertex_type(vert,:)==PERIODIC_MASTER)     .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_MASTER_DIR) .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_MASTER_NAT) .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_MASTER_MIX) .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_SLAVE)      .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_SLAVE_DIR)  .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_SLAVE_NAT)  .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_SLAVE_MIX)
endif

end function is_periodic_vert

!        ----------------
function is_periodic_edge(edge,grid,comp)
!        ----------------

!----------------------------------------------------
! This routine returns true if the given edge is periodic.
! If comp is present, it only checks component comp; otherwise it is periodic
! if any component is periodic.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
integer, intent(in) :: edge
type(grid_type), intent(in) :: grid
integer, intent(in), optional :: comp
logical :: is_periodic_edge
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (.not. grid%any_periodic) then
   is_periodic_edge = .false.
   return
endif

if (present(comp)) then
   is_periodic_edge=grid%edge_type(edge,comp)==PERIODIC_MASTER .or. &
                    grid%edge_type(edge,comp)==PERIODIC_SLAVE
else
   is_periodic_edge=any(grid%edge_type(edge,:)==PERIODIC_MASTER) .or. &
                    any(grid%edge_type(edge,:)==PERIODIC_SLAVE)
endif

end function is_periodic_edge

!        ----------------
function is_periodic_face(face,grid,comp)
!        ----------------

!----------------------------------------------------
! In 3D this routine returns true if the given face is periodic.
! In 2D it returns false.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
integer, intent(in) :: face
type(grid_type), intent(in) :: grid
integer, intent(in), optional :: comp
logical :: is_periodic_face
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

is_periodic_face = .false.

end function is_periodic_face

!        -----------------------
function is_periodic_vert_master(vert,grid,comp)
!        -----------------------

!----------------------------------------------------
! This routine returns true if the given vertex is periodic and the master.
! If comp is present, it only checks component comp; otherwise it is periodic
! if any component is periodic.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
integer, intent(in) :: vert
type(grid_type), intent(in) :: grid
integer, intent(in), optional :: comp
logical :: is_periodic_vert_master
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (.not. grid%any_periodic) then
   is_periodic_vert_master = .false.
   return
endif

if (present(comp)) then
   is_periodic_vert_master = &
                    grid%vertex_type(vert,comp)==PERIODIC_MASTER     .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_MASTER_DIR .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_MASTER_NAT .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_MASTER_MIX
else
   is_periodic_vert_master = &
                    any(grid%vertex_type(vert,:)==PERIODIC_MASTER)     .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_MASTER_DIR) .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_MASTER_NAT) .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_MASTER_MIX)
endif

end function is_periodic_vert_master

!        -----------------------
function is_periodic_edge_master(edge,grid,comp)
!        -----------------------

!----------------------------------------------------
! This routine returns true if the given edge is periodic and the master.
! If comp is present, it only checks component comp; otherwise it is periodic
! if any component is periodic.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
integer, intent(in) :: edge
type(grid_type), intent(in) :: grid
integer, intent(in), optional :: comp
logical :: is_periodic_edge_master
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (.not. grid%any_periodic) then
   is_periodic_edge_master = .false.
   return
endif

if (present(comp)) then
   is_periodic_edge_master = grid%edge_type(edge,comp)==PERIODIC_MASTER
else
   is_periodic_edge_master = any(grid%edge_type(edge,:)==PERIODIC_MASTER)
endif

end function is_periodic_edge_master

!        -----------------------
function is_periodic_face_master(face,grid,comp)
!        -----------------------

!----------------------------------------------------
! In 3D this routine returns true if the given face is periodic.
! In 2D it returns false.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
integer, intent(in) :: face
type(grid_type), intent(in) :: grid
integer, intent(in), optional :: comp
logical :: is_periodic_face_master
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

is_periodic_face_master = .false.

end function is_periodic_face_master

!        ----------------------
function is_periodic_vert_slave(vert,grid,comp)
!        ----------------------

!----------------------------------------------------
! This routine returns true if the given vertex is periodic and a slave.
! If comp is present, it only checks component comp; otherwise it is periodic
! if any component is periodic.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
integer, intent(in) :: vert
type(grid_type), intent(in) :: grid
integer, intent(in), optional :: comp
logical :: is_periodic_vert_slave
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (.not. grid%any_periodic) then
   is_periodic_vert_slave = .false.
   return
endif

if (present(comp)) then
   is_periodic_vert_slave = &
                    grid%vertex_type(vert,comp)==PERIODIC_SLAVE     .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_SLAVE_DIR .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_SLAVE_NAT .or. &
                    grid%vertex_type(vert,comp)==PERIODIC_SLAVE_MIX
else
   is_periodic_vert_slave = &
                    any(grid%vertex_type(vert,:)==PERIODIC_SLAVE)     .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_SLAVE_DIR) .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_SLAVE_NAT) .or. &
                    any(grid%vertex_type(vert,:)==PERIODIC_SLAVE_MIX)
endif

end function is_periodic_vert_slave

!        ----------------------
function is_periodic_edge_slave(edge,grid,comp)
!        ----------------------

!----------------------------------------------------
! This routine returns true if the given edge is periodic and a slave.
! If comp is present, it only checks component comp; otherwise it is periodic
! if any component is periodic.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
integer, intent(in) :: edge
type(grid_type), intent(in) :: grid
integer, intent(in), optional :: comp
logical :: is_periodic_edge_slave
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (.not. grid%any_periodic) then
   is_periodic_edge_slave = .false.
   return
endif

if (present(comp)) then
   is_periodic_edge_slave = grid%edge_type(edge,comp)==PERIODIC_SLAVE
else
   is_periodic_edge_slave = any(grid%edge_type(edge,:)==PERIODIC_SLAVE)
endif

end function is_periodic_edge_slave

!        ----------------------
function is_periodic_face_slave(face,grid,comp)
!        ----------------------

!----------------------------------------------------
! In 3D this routine returns true if the given face is periodic.
! In 2D it returns false.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
integer, intent(in) :: face
type(grid_type), intent(in) :: grid
integer, intent(in), optional :: comp
logical :: is_periodic_face_slave
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

is_periodic_face_slave = .false.

end function is_periodic_face_slave

!          -----------------------
subroutine find_containing_element(x,y,have_it,elem,grid,starting_elem,z)
!          -----------------------

!----------------------------------------------------
! This routine looks for an element containing the point (x,y).  If it
! determines that (x,y) is not in an element owned by this processor, it
! quits looking and returns (have_it=0,elem=0).  Otherwise it returns
! have_it=1 and the element index in elem.  If the point falls on an
! element boundary, it is indeterminant as to which of the containing
! elements it returns.
! If starting_elem is present, the point must be in that (possibly non-leaf)
! element.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(out) :: have_it,elem
type(grid_type), intent(in) :: grid
integer, optional :: starting_elem
real(my_real), optional, intent(in) :: z

!----------------------------------------------------
! Local variables:

integer :: e,root,i,neigh(NEIGHBORS_PER_ELEMENT)
real(my_real) :: bc(VERTICES_PER_ELEMENT,1)
integer, parameter :: roundoff_fudge = 1000
!----------------------------------------------------
! Begin executable code

! find an element in the initial grid that contains (x,y)

if (present(starting_elem)) then
   root = starting_elem
else

! first look for it by moving from a triangle to a neighbor in the direction
! of the point, which is characterized by a negative barycentric coordinate

   root = -10
   e = grid%head_level_elem(1)
   do
! if all barycentric coordinates are positive, it's in there
      call barycentric((/x/),(/y/), &
                       grid%vertex(grid%element(e)%vertex)%coord%x, &
                       grid%vertex(grid%element(e)%vertex)%coord%y, &
                       bc,no_det=.true.)
      if (all(bc >= -roundoff_fudge*epsilon(0.0_my_real))) then
         root = e
         exit
      endif
      neigh = grid%initial_neighbor(:,e)
      do i=1,VERTICES_PER_ELEMENT
         if (bc(i,1) < 0) then
            e = neigh(i)
            if (e /= BOUNDARY) exit
         endif
      end do
      if (e == BOUNDARY) exit
   end do
endif ! present starting_elem

! if that failed, go through all the initial elements

if (root == -10) then
   e = grid%head_level_elem(1)
   do while (e /= END_OF_LIST)
      call barycentric((/x/),(/y/), &
                       grid%vertex(grid%element(e)%vertex)%coord%x, &
                       grid%vertex(grid%element(e)%vertex)%coord%y, &
                       bc,no_det=.true.)
      if (all(bc >= -roundoff_fudge*epsilon(0.0_my_real))) then
         root = e
         exit
      endif
      e = grid%element(e)%next
   end do
endif

! go down the refinement tree to find the leaf element that contains (x,y)

if (root /= -10) then
   call find_containing_leaf(x,y,root,have_it,elem,grid,z)
else
   have_it = 0
   elem = 0
endif

end subroutine find_containing_element

!                    --------------------
recursive subroutine find_containing_leaf(x,y,root,have_it,elem,grid,z)
!                    --------------------

!----------------------------------------------------
! This routine recursively goes down the refinement tree, starting at root,
! to find a leaf element containing (x,y) as long as (x,y) may be in an
! element I own.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: root
integer, intent(out) :: have_it,elem
type(grid_type), intent(in) :: grid
real(my_real), optional, intent(in) :: z

!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
real(my_real) :: bc(VERTICES_PER_ELEMENT,1)
integer, parameter :: roundoff_fudge = 1000
!----------------------------------------------------
! Begin executable code

! if root is a leaf, we've found it; but verify ownership

allc = ALL_CHILDREN
children = get_child_lid(grid%element(root)%gid,allc,grid%elem_hash)
if (children(1) == NO_CHILD) then
   if (grid%element(root)%iown) then
      have_it = 1
      elem = root
   else
      have_it = 0
      elem = 0
   endif
   return
endif

! otherwise, look at the barycentric coordinates of (x,y) in each child,
! until one is found where they are all positive

have_it = 0
elem = 0
do i=1,MAX_CHILD
   call barycentric((/x/),(/y/), &
                    grid%vertex(grid%element(children(i))%vertex)%coord%x, &
                    grid%vertex(grid%element(children(i))%vertex)%coord%y, &
                    bc,no_det=.true.)
   if (all(bc >= -roundoff_fudge*epsilon(0.0_my_real))) then
      call find_containing_leaf(x,y,children(i),have_it,elem,grid,z)
      exit
   endif
end do

end subroutine find_containing_leaf

!        --------------------
function get_child_gid_scalar(gid,child)
!        --------------------

!----------------------------------------------------
! This function returns the global ID of child number child of
! element number gid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: gid
integer, intent(in) :: child
type(hash_key) :: get_child_gid_scalar

!----------------------------------------------------
! Local variables: 

!----------------------------------------------------
! Begin executable code

get_child_gid_scalar = MAX_CHILD*gid+(child-1)

end function get_child_gid_scalar

!        -------------------
function get_child_gid_array(gid,child)
!        -------------------

!----------------------------------------------------
! This function returns the global ID of the children of element
! number gid listed in child
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: gid
integer, intent(in) :: child(:)
type(hash_key) :: get_child_gid_array(size(child))

!----------------------------------------------------
! Local variables: 

integer :: i

!----------------------------------------------------
! Begin executable code

do i=1,size(child)
   get_child_gid_array(i) = MAX_CHILD*gid+(child(i)-1)
end do

end function get_child_gid_array

!        --------------------
function get_child_lid_scalar(gid,child,table)
!        --------------------

!----------------------------------------------------
! This function returns the local ID of child number child of
! element number gid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: gid
integer, intent(in) :: child
type(hash_table), intent(in) :: table
integer :: get_child_lid_scalar

!----------------------------------------------------
! Local variables: 

!----------------------------------------------------
! Begin executable code

if (hash_overflow(gid,MAX_CHILD,child-1)) then
   get_child_lid_scalar = NO_CHILD
else
   get_child_lid_scalar = hash_decode_key(MAX_CHILD*gid+(child-1),table)
   if (get_child_lid_scalar == HASH_NOT_FOUND) then
      get_child_lid_scalar = NO_CHILD
   endif
endif

end function get_child_lid_scalar

!        -------------------
function get_child_lid_array(gid,child,table)
!        -------------------

!----------------------------------------------------
! This function returns the local ID of the children of element
! number gid listed in child
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: gid
integer, intent(in) :: child(:)
type(hash_table), intent(in) :: table
integer :: get_child_lid_array(size(child))

!----------------------------------------------------
! Local variables: 

integer :: i

!----------------------------------------------------
! Begin executable code

do i=1,size(child)
   if (hash_overflow(gid,MAX_CHILD,child(i)-1)) then
      get_child_lid_array(i) = NO_CHILD
   else
      get_child_lid_array(i)=hash_decode_key(MAX_CHILD*gid+(child(i)-1), &
                                             table)
      if (get_child_lid_array(i) == HASH_NOT_FOUND) then
         get_child_lid_array(i) = NO_CHILD
      endif
   endif
end do

end function get_child_lid_array

!        -------------
function get_neighbors(lid,grid)
!        -------------

!----------------------------------------------------
! This function returns the local IDs of the neighbors of element lid.
! Order is (sibling,child of parent's mate,mate) i.e. opposite vert 1,2,3.
! On level 1, the order is the same as in initial_neighbor.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: lid
type(grid_type), intent(in) :: grid
integer :: get_neighbors(NEIGHBORS_PER_ELEMENT)

!----------------------------------------------------
! Local variables: 

integer :: parent_lid, mate_lid, child_lid(MAX_CHILD)
type(hash_key) :: gid, parent, child1, parent_mate_gid, gid_neighbor
type(hash_key) :: children(MAX_CHILD)
integer :: allc(MAX_CHILD)
integer :: i, j, k, vert, other_vert(1), edge, other_edge(EDGES_PER_ELEMENT)

!----------------------------------------------------
! Begin executable code

! In the case of periodic boundary conditions, identify the matching edges
! and vertex

do i=1,EDGES_PER_ELEMENT
   edge = grid%element(lid)%edge(i)
   if (is_periodic_edge(edge,grid)) then
      other_edge(i) = matching_periodic_edge(grid,edge)
   else
      other_edge(i) = edge
   endif
end do
vert = grid%element(lid)%vertex(1)
if (is_periodic_vert(vert,grid)) then
   other_vert(1) = matching_periodic_vert(grid,vert)
else
   other_vert(1) = vert
endif

! Neighbors of a level 1 element.

if (grid%element(lid)%level == 1) then

! Start with the neighbors on level 1.

   get_neighbors = grid%initial_neighbor(:,lid)

! For each neighbor ...

   do i=1,NEIGHBORS_PER_ELEMENT

! If the initial neighbor was refined, get the child that shares an edge.

      if (get_neighbors(i) /= BOUNDARY) then
         gid_neighbor = grid%element(get_neighbors(i))%gid
         if (hash_overflow(gid_neighbor,MAX_CHILD,1)) then
            child_lid = HASH_NOT_FOUND
         else
            allc = ALL_CHILDREN
            children = get_child_gid(gid_neighbor,allc)
            child_lid(1) = hash_decode_key(children(1),grid%elem_hash)
            child_lid(2) = hash_decode_key(children(2),grid%elem_hash)
         endif
         if (child_lid(1) /= HASH_NOT_FOUND) then
            do j=1,EDGES_PER_ELEMENT
               do k=1,MAX_CHILD
                  if (any(grid%element(child_lid(k))%edge == &
                          grid%element(lid)%edge(j)) .or. &
                      any(grid%element(child_lid(k))%edge == other_edge(j))) then
                     get_neighbors(i) = child_lid(k)
                     exit
                  endif
               end do
            end do
         endif
      endif
   end do

! Elements that are not level 1

else

   gid = grid%element(lid)%gid

! Identify the parent and the first child of the parent

   parent = gid/MAX_CHILD
   parent_lid = hash_decode_key(parent,grid%elem_hash)
   child1 = MAX_CHILD*parent

! The sibling is the child of the parent that is not this element

   if (child1 == gid) then
      gid_neighbor = child1+1
   else
      gid_neighbor = child1
   endif

! If the sibling was refined, get the child whose vertex 1 is my vertex 2

   if (hash_overflow(gid_neighbor,MAX_CHILD,1)) then
      child_lid = HASH_NOT_FOUND
   else
      allc = ALL_CHILDREN
      children = get_child_gid(gid_neighbor,allc)
      child_lid(1) = hash_decode_key(children(1),grid%elem_hash)
   endif
   if (child_lid(1) /= HASH_NOT_FOUND) then
      if (grid%element(child_lid(1))%vertex(1) == &
          grid%element(lid)%vertex(2)) then
         get_neighbors(1) = child_lid(1)
         gid_neighbor = BOUNDARY ! so I don't hash decode it later
      else
         gid_neighbor = children(2)
      endif
   endif
   if (.not. gid_neighbor == BOUNDARY) then
      get_neighbors(1) = hash_decode_key(gid_neighbor,grid%elem_hash)
   endif

! Get the neighboring child of the parent's mate

   parent_mate_gid = grid%element(parent_lid)%mate
   if (parent_mate_gid == BOUNDARY) then
      get_neighbors(2) = BOUNDARY
      gid_neighbor = BOUNDARY
   else
      if (child1 == gid) then
         gid_neighbor = MAX_CHILD*parent_mate_gid
      else
         gid_neighbor = MAX_CHILD*parent_mate_gid+1
      endif

! If the neighbor was refined, get the child whose vertex 1 is my vertex 1.

      if (hash_overflow(gid_neighbor,MAX_CHILD,1)) then
         child_lid(1) = HASH_NOT_FOUND
      else
         allc = ALL_CHILDREN
         children = get_child_gid(gid_neighbor,allc)
         child_lid(1) = hash_decode_key(children(1),grid%elem_hash)
      endif
      if (child_lid(1) /= HASH_NOT_FOUND) then
         if (grid%element(child_lid(1))%vertex(1)==grid%element(lid)%vertex(1).or.&
             grid%element(child_lid(1))%vertex(1)==other_vert(1)) then
            get_neighbors(2) = child_lid(1)
            gid_neighbor = BOUNDARY
         else
            gid_neighbor = children(2)
         endif
      endif
   endif
   if (.not. gid_neighbor == BOUNDARY) then
      get_neighbors(2) = hash_decode_key(gid_neighbor,grid%elem_hash)
   endif

! Get the mate

   if (grid%element(lid)%mate == BOUNDARY) then
      get_neighbors(3) = BOUNDARY
   else
      mate_lid = hash_decode_key(grid%element(lid)%mate,grid%elem_hash)
      if (mate_lid == HASH_NOT_FOUND) then
         get_neighbors(3)=hash_decode_key(grid%element(lid)%mate/MAX_CHILD, &
                                          grid%elem_hash)
      else
         get_neighbors(3) = mate_lid
      endif
   endif
endif

end function get_neighbors


!                    -------------------
recursive subroutine get_vertex_elements(grid,vert,elements,include_periodic)
!                    -------------------

!----------------------------------------------------
! This routine returns a list of all the elements around vertex vert.
! If include_periodic is present and true, it also includes elements
! around periodic mates.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: vert
integer, pointer :: elements(:)
logical, intent(in), optional :: include_periodic
!----------------------------------------------------
! Local variables:

! assume there aren't more than 100 neighbors
integer, parameter :: MAX_NEIGH = 100
integer :: nelem, loc_elements(MAX_NEIGH), neighbors(NEIGHBORS_PER_ELEMENT), &
           neigh, ineigh, i, next, astat, periodic_vert
logical :: loc_periodic
integer, pointer :: periodic_elements(:)
!----------------------------------------------------
! Begin executable code

if (present(include_periodic)) then
   loc_periodic = include_periodic
else
   loc_periodic = .false.
endif

! begin with associated element

loc_elements(1) = grid%vertex(vert)%assoc_elem
nelem = 1
next = 1

do while (next <= nelem)

! get the neighbors of the next element in the list

   neighbors = get_neighbors(loc_elements(next),grid)
   next = next + 1

! for each neighbor ...

   outer: do ineigh=1,NEIGHBORS_PER_ELEMENT
      neigh = neighbors(ineigh)
      if (neigh == BOUNDARY) cycle

! if it is already on the list, move on

      do i=1,nelem
         if (loc_elements(i) == neigh) cycle outer
      end do

! if any vertex is vert, add it to the list

      if (any(grid%element(neigh)%vertex == vert)) then
         nelem = nelem + 1
         if (nelem > MAX_NEIGH) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("too many elements in get_vertex_elements")
            stop
         endif
         loc_elements(nelem) = neigh
      endif

   end do outer

end do

! if the vertex is periodic and periodic elements are to be included,
! include the elements around each periodic mate

if (loc_periodic) then
   if (is_periodic_vert(vert,grid)) then

      periodic_vert = grid%vertex(vert)%next
      do while (periodic_vert /= vert)
         call get_vertex_elements(grid,periodic_vert,periodic_elements)
         call merge_periodic_vert
         deallocate(periodic_elements)
         periodic_vert = grid%vertex(periodic_vert)%next
      end do

   endif
endif

! copy the element list to the result

allocate(elements(nelem),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in get_vertex_elements")
   stop
endif

elements = loc_elements(1:nelem)

contains

subroutine merge_periodic_vert
integer :: ii
do ii=1,size(periodic_elements)
   if (.not. any(loc_elements(1:nelem) == periodic_elements(ii))) then
      nelem = nelem + 1
      if (nelem > MAX_NEIGH) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("too many elements in get_vertex_elements")
         stop
      endif
      loc_elements(nelem) = periodic_elements(ii)
   endif
end do
end subroutine merge_periodic_vert

end subroutine get_vertex_elements

!                    -----------------
recursive subroutine get_edge_elements(grid,edge,elements,include_periodic)
!                    -----------------

!----------------------------------------------------
! This routine returns a list of all the elements around edge edge.
! If include_periodic is present and true, it also includes elements
! around periodic mates.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: edge
integer, pointer :: elements(:)
logical, intent(in), optional :: include_periodic
!----------------------------------------------------
! Local variables:

! assume there aren't more than 100 neighbors
integer, parameter :: MAX_NEIGH = 100
integer :: nelem, loc_elements(MAX_NEIGH), neighbors(NEIGHBORS_PER_ELEMENT), &
           neigh, ineigh, i, next, astat, periodic_edge
logical :: loc_periodic
integer, pointer :: periodic_elements(:)
!----------------------------------------------------
! Begin executable code

if (present(include_periodic)) then
   loc_periodic = include_periodic
else
   loc_periodic = .false.
endif

! begin with associated element

loc_elements(1) = grid%edge(edge)%assoc_elem
nelem = 1
next = 1

do while (next <= nelem)

! get the neighbors of the next element in the list

   neighbors = get_neighbors(loc_elements(next),grid)
   next = next + 1

! for each neighbor ...

   outer: do ineigh=1,NEIGHBORS_PER_ELEMENT
      neigh = neighbors(ineigh)
      if (neigh == BOUNDARY) cycle

! if it is already on the list, move on

      do i=1,nelem
         if (loc_elements(i) == neigh) cycle outer
      end do

! if any edge is edge, add it to the list

      if (any(grid%element(neigh)%edge == edge)) then
         nelem = nelem + 1
         if (nelem > MAX_NEIGH) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("too many elements in get_edge_elements")
            stop
         endif
         loc_elements(nelem) = neigh
      endif

   end do outer

end do

! if the edge is periodic and periodic elements are to be included,
! include the elements around the periodic mate

if (loc_periodic) then
   if (is_periodic_edge(edge,grid)) then

      periodic_edge = grid%edge(edge)%next
      call get_edge_elements(grid,periodic_edge,periodic_elements)
      call merge_periodic_edge
      deallocate(periodic_elements)

   endif
endif

! copy the element list to the result

allocate(elements(nelem),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in get_edge_elements")
   stop
endif

elements = loc_elements(1:nelem)

contains

subroutine merge_periodic_edge
integer :: ii
do ii=1,size(periodic_elements)
   if (.not. any(loc_elements(1:nelem) == periodic_elements(ii))) then
      nelem = nelem + 1
      if (nelem > MAX_NEIGH) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("too many elements in get_edge_elements")
         stop
      endif
      loc_elements(nelem) = periodic_elements(ii)
   endif
end do
end subroutine merge_periodic_edge

end subroutine get_edge_elements

!          -----------------
subroutine get_face_elements(grid,face,elements,include_periodic)
!          -----------------

!----------------------------------------------------
! This routine would return a list of all the elements around face face,
! except there are no faces in 2D.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: face
integer, pointer :: elements(:)
logical, intent(in), optional :: include_periodic
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call warning("get_face_elements called in a 2D compilation")
nullify(elements)

end subroutine get_face_elements

!          -------------
subroutine get_grid_info(grid,procs,still_sequential,tag,nelem,nlev,nvert, &
                         nvert_own,nelem_own,nelem_leaf,nelem_leaf_own,     &
                         dof,dof_own,total_nvert,total_nelem_leaf,total_dof, &
                         max_nlev,mindeg,maxdeg,no_master)
!          -------------

!----------------------------------------------------
! This routine returns information about the grid
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
type (proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: tag
integer, intent(out), optional :: nelem,nlev,nvert,nvert_own, &
                                  nelem_own,nelem_leaf,nelem_leaf_own, &
                                  dof,dof_own,total_nvert,total_nelem_leaf, &
                                  total_dof,max_nlev,mindeg,maxdeg
logical, intent(in), optional :: no_master
!----------------------------------------------------

integer :: proc,ni,nr,astat,lev,elem,min_mindeg,max_maxdeg
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
real(my_real) :: no_reals(1)
logical :: send_to_master

!----------------------------------------------------
! Begin executable code

if (present(no_master)) then
   if (no_master) then
      send_to_master = .false.
   else
      send_to_master = .true.
   endif
else
   send_to_master = .true.
endif

if (my_proc(procs)==MASTER) then

   if (present(nelem)) nelem = 0
   if (present(nlev)) nlev = 0
   if (present(nvert)) nvert = 0
   if (present(nvert_own)) nvert_own = 0
   if (present(nelem_own)) nelem_own = 0
   if (present(nelem_leaf)) nelem_leaf = 0
   if (present(nelem_leaf_own)) nelem_leaf_own = 0
   if (present(dof)) dof = 0
   if (present(dof_own)) dof_own = 0
   if (present(total_nvert)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+4)
      total_nvert = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if
   if (present(total_nelem_leaf)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+6)
      total_nelem_leaf = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if
   if (present(total_dof)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+9)
      total_dof = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if
   if (present(max_nlev)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+8)
      max_nlev = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if
   if (present(mindeg)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+11)
      mindeg = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if
   if (present(maxdeg)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+13)
      maxdeg = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if

else ! not processor 0

   if (present(nelem)) nelem = grid%nelem
   if (present(nlev)) nlev = grid%nlev
   if (present(nvert)) nvert = grid%nvert
   if (present(nvert_own)) nvert_own = grid%nvert_own
   if (present(nelem_own)) nelem_own = 0 ! was grid%nelem_own
   if (present(nelem_leaf)) nelem_leaf = grid%nelem_leaf
   if (present(nelem_leaf_own)) nelem_leaf_own = grid%nelem_leaf_own
   if (present(dof)) dof = grid%dof
   if (present(dof_own)) dof_own = grid%dof_own

   if (present(total_nvert)) then
      if (still_sequential) then
         total_nvert = grid%nvert_own
      else
         total_nvert = phaml_global_sum(procs,grid%nvert_own,10*tag+1)
      endif
      if (my_proc(procs)==1 .and. send_to_master) then
         call phaml_send(procs,MASTER,(/total_nvert/),1,no_reals,0,10*tag+4)
      endif
   end if
   if (present(total_nelem_leaf)) then
      if (still_sequential) then
         total_nelem_leaf = grid%nelem_leaf_own
      else
         total_nelem_leaf = phaml_global_sum(procs,grid%nelem_leaf_own,10*tag+3)
      endif
      if (my_proc(procs)==1 .and. send_to_master) then
         call phaml_send(procs,MASTER,(/total_nelem_leaf/),1,no_reals,0,10*tag+6)
      endif
   end if
   if (present(total_dof)) then
      if (still_sequential) then
         total_dof = grid%dof_own
      else
         total_dof = phaml_global_sum(procs,grid%dof_own,10*tag+2)
      endif
      if (my_proc(procs)==1 .and. send_to_master) then
         call phaml_send(procs,MASTER,(/total_dof/),1,no_reals,0,10*tag+9)
      endif
   end if
   if (present(max_nlev)) then
      if (still_sequential) then
         max_nlev = grid%nlev
      else
         max_nlev = phaml_global_max(procs,grid%nlev,10*tag+7)
      endif
      if (my_proc(procs)==1 .and. send_to_master) then
         call phaml_send(procs,MASTER,(/max_nlev/),1,no_reals,0,10*tag+8)
      endif
   endif
   if (present(mindeg) .or. present(maxdeg)) then
      if (present(mindeg)) mindeg = huge(0)
      if (present(maxdeg)) maxdeg = 0
      do lev=1,grid%nlev
         elem = grid%head_level_elem(lev)
         do while (elem /= END_OF_LIST)
            if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
               if (present(mindeg)) mindeg=min(mindeg,grid%element(elem)%degree)
               if (present(maxdeg)) maxdeg=max(maxdeg,grid%element(elem)%degree)
            endif
            elem = grid%element(elem)%next
         end do
      end do
      if (send_to_master) then
         if (present(mindeg)) then
            min_mindeg = phaml_global_min(procs,mindeg,10*tag+10)
            if (my_proc(procs)==1) then
               call phaml_send(procs,MASTER,(/min_mindeg/),1,no_reals,0,10*tag+11)
            endif
         endif
         if (present(maxdeg)) then
            max_maxdeg = phaml_global_max(procs,maxdeg,10*tag+12)
            if (my_proc(procs)==1) then
               call phaml_send(procs,MASTER,(/max_maxdeg/),1,no_reals,0,10*tag+13)
            endif
         endif
      endif
   endif

endif

end subroutine get_grid_info

!          -------------------------
subroutine get_vertices_and_solution(grid,x,y,u,comp,eigen,z)
!          -------------------------

!----------------------------------------------------
! This routine returns the coordinates of the vertices owned by this processor
! and the (comp,eigen) solution at them.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), pointer :: x(:), y(:), u(:)
integer, intent(in) :: comp, eigen
real(my_real), pointer, optional :: z(:)
!----------------------------------------------------
! Local variables:

integer :: i, lev, vert, astat, elem, ivert
logical(small_logical) :: visited_vert(grid%biggest_vert)
!----------------------------------------------------
! Begin executable code

if (present(z)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("get_vertices_and_solution: z present in call to 2D code")
   stop
endif

allocate(x(grid%nvert_own),y(grid%nvert_own),u(grid%nvert_own),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in get_vertices_and_solution")
   stop
endif

i = 0
visited_vert = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do ivert=1,VERTICES_PER_ELEMENT
         vert = grid%element(elem)%vertex(ivert)
         if (visited_vert(vert)) cycle
         visited_vert(vert) = .true.
         if (grid%element(grid%vertex(vert)%assoc_elem)%iown) then
            i = i+1
            if (i > grid%nvert_own) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("insufficient space in get_vertices_and_solution")
               stop
            endif
            x(i) = grid%vertex(vert)%coord%x
            y(i) = grid%vertex(vert)%coord%y
            u(i) = grid%vertex_solution(vert,comp,eigen)
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

if (i /= grid%nvert_own) then
   call warning("number of vertices mismatch in get_verticies_and_solution")
endif

end subroutine get_vertices_and_solution

!          ---------------
subroutine check_pde_terms(grid,has_first_order,has_cross,laplacian)
!          ---------------

!----------------------------------------------------
! This routine attempts to determine certain characteristics of the pde
! by examining the terms at the midpoints of the initial elements
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
logical, intent(out) :: has_first_order, has_cross, laplacian

!----------------------------------------------------
! Local variables:

integer :: elem
real(my_real) :: x,y,z
real(my_real) :: cxx(grid%system_size,grid%system_size), &
                 cyy(grid%system_size,grid%system_size), &
                 czz(grid%system_size,grid%system_size), &
                 cxy(grid%system_size,grid%system_size), &
                 cxz(grid%system_size,grid%system_size), &
                 cyz(grid%system_size,grid%system_size), &
                 cx(grid%system_size,grid%system_size),  &
                 cy(grid%system_size,grid%system_size),  &
                 cz(grid%system_size,grid%system_size),  &
                 c(grid%system_size,grid%system_size),   &
                 rs(grid%system_size)
!----------------------------------------------------
! Begin executable code

has_first_order = .false.
has_cross = .false.
laplacian = .true.

elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   x = sum(grid%vertex(grid%element(elem)%vertex)%coord%x)/VERTICES_PER_ELEMENT
   y = sum(grid%vertex(grid%element(elem)%vertex)%coord%y)/VERTICES_PER_ELEMENT
   z = sum(zcoord(grid%vertex(grid%element(elem)%vertex)%coord))/VERTICES_PER_ELEMENT
   call my_pdecoefs(elem,x,y,z,cxx,cyy,czz,cxy,cxz,cyz,cx,cy,cz,c,rs)
   if (any(cx /= 0.0_my_real) .or. any(cy /= 0.0_my_real) .or. &
       any(cz /= 0.0_my_real)) then
      has_first_order = .true.
      laplacian = .false.
      if (has_cross) exit
   endif
   if (any(cxy /= 0.0_my_real) .or. any(cxz /= 0.0_my_real) .or. &
       any(cyz /= 0.0_my_real)) then
      has_cross = .true.
      laplacian = .false.
      if (has_first_order) exit
   endif
   if (any(cxx /= 1.0_my_real) .or. any(cyy /= 1.0_my_real) .or. &
       any(c   /= 0.0_my_real)) then
      laplacian = .false.
      if (has_cross .and. has_first_order) exit
   endif
   elem = grid%element(elem)%next
end do

end subroutine check_pde_terms

!          ---------------------
subroutine check_triangle_shapes(grid,iso_right_tri)
!          ---------------------

!----------------------------------------------------
! This routine checks to see if all triangles on level 1 are isoceles right
! triangles with the third vertex at the right angle.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
logical, intent(out) :: iso_right_tri
!----------------------------------------------------
! Local variables:

real(my_real), parameter :: small = 100*epsilon(0.0_my_real)
integer :: elem
type(point) :: vec1, vec2
!----------------------------------------------------
! Begin executable code

iso_right_tri = .true.

elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)

! Compute the vectors from vertex 3 to vertex 1 and from vertex 3 to vertex 2

   vec1 = grid%vertex(grid%element(elem)%vertex(1))%coord - &
          grid%vertex(grid%element(elem)%vertex(3))%coord
   vec2 = grid%vertex(grid%element(elem)%vertex(2))%coord - &
          grid%vertex(grid%element(elem)%vertex(3))%coord

! Compute the inner product of the vectors.  If it is not (close to) 0, it is
! not a right angle.

   if (abs(dot_point(vec1,vec2)) > small) then
      iso_right_tri = .false.
      exit
   endif

! Compute the (squares of) the length of the vectors.  If they are not (close
! to) equal, it is not isoceles.

   if (abs(dot_point(vec1,vec1)-dot_point(vec2,vec2)) > small) then
      iso_right_tri = .false.
      exit
   endif

   elem = grid%element(elem)%next
end do

end subroutine check_triangle_shapes

!        ------------------
function get_next_free_vert(grid,errcode)
!        ------------------

!----------------------------------------------------
! This routine returns the next available local ID for a vertex.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(out) :: errcode
integer :: get_next_free_vert
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

errcode = 0

! next_free_vert always has the next one

get_next_free_vert = grid%next_free_vert

! if the next one is bigger than what is available, increase allocation

if (grid%next_free_vert > size(grid%vertex)) then
   call more_verts(grid,errcode)
   if (errcode /= 0) return
endif

! if the next one is bigger than ever used, the next next one follows it

if (grid%next_free_vert > grid%biggest_vert) then
   grid%next_free_vert = grid%next_free_vert + 1

! otherwise, it is the next one in the linked list

else
   grid%next_free_vert = grid%vertex(grid%next_free_vert)%next

endif

! check for new biggest ID

if (get_next_free_vert > grid%biggest_vert) then
   grid%biggest_vert = get_next_free_vert
endif

end function get_next_free_vert

!        ------------------
function get_next_free_edge(grid,errcode)
!        ------------------

!----------------------------------------------------
! This routine returns the next available local ID for an edge.
! For comments, see get_next_free_vert.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(out) :: errcode
integer :: get_next_free_edge
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

errcode = 0
get_next_free_edge = grid%next_free_edge
if (grid%next_free_edge > size(grid%edge)) then
   call more_edges(grid,errcode)
   if (errcode /= 0) return
endif
if (grid%next_free_edge > grid%biggest_edge) then
   grid%next_free_edge = grid%next_free_edge + 1
else
   grid%next_free_edge = grid%edge(grid%next_free_edge)%next
endif
if (get_next_free_edge > grid%biggest_edge) then
   grid%biggest_edge = get_next_free_edge
endif

end function get_next_free_edge

!        ------------------
function get_next_free_elem(grid,errcode,elist,reftype,new_p,numhref,numpref, &
                            desired_level,desired_degree,elem_list)
!        ------------------

!----------------------------------------------------
! This routine returns the next available local ID for an element.
! For comments, see get_next_free_vert.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(out) :: errcode
type(errind_list), optional, intent(inout) :: elist
character(len=*), optional, pointer :: reftype(:)
integer, optional, pointer :: numhref(:), numpref(:), new_p(:,:), &
                              desired_level(:), desired_degree(:), elem_list(:)
integer :: get_next_free_elem
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

errcode = 0
get_next_free_elem = grid%next_free_elem
if (grid%next_free_elem > size(grid%element)) then
   call more_elements(grid,errcode,elist,reftype,new_p,numhref,numpref, &
                         desired_level,desired_degree,elem_list)
   if (errcode /= 0) return
endif
if (grid%next_free_elem > grid%biggest_elem) then
   grid%next_free_elem = grid%next_free_elem + 1
else
   grid%next_free_elem = grid%element(grid%next_free_elem)%next
endif
if (get_next_free_elem > grid%biggest_elem) then
   grid%biggest_elem = get_next_free_elem
endif

end function get_next_free_elem

!          ------------------
subroutine put_next_free_vert(grid,vert)
!          ------------------

!----------------------------------------------------
! This routine adds vert to the list of available vertex IDs
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: vert
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! add vert to the beginning of the linked list of free space

grid%vertex(vert)%next = grid%next_free_vert
grid%next_free_vert = vert

end subroutine put_next_free_vert

!          ------------------
subroutine put_next_free_edge(grid,edge)
!          ------------------

!----------------------------------------------------
! This routine adds edge to the list of available edge IDs
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: edge
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! add edge to the beginning of the linked list of free space

grid%edge(edge)%next = grid%next_free_edge
grid%next_free_edge = edge

end subroutine put_next_free_edge

!          ------------------
subroutine put_next_free_elem(grid,elem)
!          ------------------

!----------------------------------------------------
! This routine adds elem to the list of available element IDs
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! add elem to the beginning of the linked list of free space

grid%element(elem)%next = grid%next_free_elem
grid%next_free_elem = elem

end subroutine put_next_free_elem

!          -----------
subroutine extend_nlev(grid)
!          -----------

!----------------------------------------------------
! This routine increases the number of refinement levels supported in
! the grid data structure
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid

!----------------------------------------------------
! Local variables:

integer, pointer :: old_array(:)
integer :: oldlev, newlev, astat, dstat

!----------------------------------------------------
! Begin executable code

oldlev = size(grid%head_level_elem)
newlev = 2*oldlev
old_array => grid%head_level_elem
allocate(grid%head_level_elem(newlev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in extend_nlev")
   grid%head_level_elem => old_array
   return
endif
grid%head_level_elem(1:oldlev) = old_array
grid%head_level_elem(oldlev+1:newlev) = END_OF_LIST
deallocate(old_array,stat=dstat)

end subroutine extend_nlev

!          -------------
subroutine more_elements(grid,errcode,elist,reftype,new_p,numhref,numpref, &
                         desired_level,desired_degree,elem_list,size_wanted)
!          -------------

!----------------------------------------------------
! This routine increases the space for the number of elements
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(out) :: errcode
type(errind_list), optional, intent(inout) :: elist
character(len=*), optional, pointer :: reftype(:)
integer, optional, pointer :: numhref(:), numpref(:), new_p(:,:), &
                              desired_level(:), desired_degree(:), elem_list(:)
integer, intent(in), optional :: size_wanted
!----------------------------------------------------
! Local variables:

type(element_t), pointer :: temp_elem(:)
integer, pointer :: temp_eprev(:), temp_enext(:)
character(len=1), pointer :: temp_reftype(:)
integer, pointer :: temp_numref(:), temp_new_p(:,:), temp_desired(:)
real(my_real), pointer :: temp_errind(:,:)
integer :: allocstat, oldsize, newsize, i
!----------------------------------------------------
! Begin executable code

oldsize = size(grid%element)
if (present(size_wanted)) then
   newsize = size_wanted
else
   newsize = int(1.5*oldsize)
endif

errcode = 0
temp_elem => grid%element
nullify(grid%element)
allocate(grid%element(newsize),stat=allocstat)
if (allocstat /= 0) then
   errcode = 1
   call fatal("increased allocation failed")
   grid%element => temp_elem
   return
endif

temp_errind => grid%element_errind
nullify(grid%element_errind)
allocate(grid%element_errind(newsize,max(1,grid%num_eval)),stat=allocstat)
if (allocstat /= 0) then
   errcode = 1
   call fatal("increased allocation failed")
   grid%element => temp_elem
   grid%element_errind => temp_errind
   return
endif

if (present(elist)) then
   temp_eprev => elist%prev_errind
   temp_enext => elist%next_errind
   nullify(elist%prev_errind, elist%next_errind)
   allocate(elist%prev_errind(newsize),elist%next_errind(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      elist%prev_errind => temp_eprev
      elist%next_errind => temp_enext
      deallocate(grid%element,stat=allocstat)
      grid%element => temp_elem
      grid%element_errind => temp_errind
      return
   endif
   elist%prev_errind(1:oldsize) = temp_eprev
   elist%next_errind(1:oldsize) = temp_enext
   deallocate(temp_eprev,temp_enext,stat=allocstat)
endif

if (present(reftype)) then
   temp_reftype => reftype
   nullify(reftype)
   allocate(reftype(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      stop
   endif
   reftype(1:oldsize) = temp_reftype
   deallocate(temp_reftype,stat=allocstat)
endif

if (present(new_p)) then
   temp_new_p => new_p
   nullify(new_p)
   allocate(new_p(2,newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      stop
   endif
   new_p(:,1:oldsize) = temp_new_p
   deallocate(temp_new_p,stat=allocstat)
endif

if (present(numhref)) then
   temp_numref => numhref
   nullify(numhref)
   allocate(numhref(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      stop
   endif
   numhref(1:oldsize) = temp_numref
   deallocate(temp_numref,stat=allocstat)
endif

if (present(numpref)) then
   temp_numref => numpref
   nullify(numpref)
   allocate(numpref(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      stop
   endif
   numpref(1:oldsize) = temp_numref
   deallocate(temp_numref,stat=allocstat)
endif

grid%element(1:oldsize) = temp_elem
grid%element_errind(1:oldsize,:) = temp_errind

deallocate(temp_elem,temp_errind,stat=allocstat)
grid%next_free_elem = oldsize+1

if (present(desired_level)) then
   temp_desired => desired_level
   nullify(desired_level)
   allocate(desired_level(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      stop
   endif
   desired_level(1:oldsize) = temp_desired
   deallocate(temp_desired,stat=allocstat)
endif

if (present(desired_degree)) then
   temp_desired => desired_degree
   nullify(desired_degree)
   allocate(desired_degree(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      stop
   endif
   desired_degree(1:oldsize) = temp_desired
   deallocate(temp_desired,stat=allocstat)
endif

if (present(elem_list)) then
   temp_desired => elem_list
   nullify(elem_list)
   allocate(elem_list(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      stop
   endif
   elem_list(1:oldsize) = temp_desired
   deallocate(temp_desired,stat=allocstat)
endif

!call count_grid_memory(grid)

end subroutine more_elements

!          ----------
subroutine more_edges(grid,errcode,size_wanted)
!          ----------

!----------------------------------------------------
! This routine increases the space for the number of edges
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(out) :: errcode
integer, intent(in), optional :: size_wanted
!----------------------------------------------------
! Local variables:

type(edge_t), pointer :: temp_edge(:)
integer, pointer :: temp_edge_type(:,:)
integer :: allocstat, oldsize, newsize, i
!----------------------------------------------------
! Begin executable code

oldsize = size(grid%edge)
if (present(size_wanted)) then
   newsize = size_wanted
else
   newsize = int(1.5*oldsize)
endif

errcode = 0
temp_edge => grid%edge
temp_edge_type => grid%edge_type
nullify(grid%edge,grid%edge_type)
allocate(grid%edge(newsize),grid%edge_type(newsize,grid%system_size), &
         stat=allocstat)
if (allocstat /= 0) then
   errcode = 1
   call fatal("increased allocation failed")
   grid%edge => temp_edge
   grid%edge_type => temp_edge_type
   return
endif
grid%edge(1:oldsize) = temp_edge
grid%edge_type(1:oldsize,:) = temp_edge_type

deallocate(temp_edge,temp_edge_type,stat=allocstat)
grid%next_free_edge = oldsize+1

!call count_grid_memory(grid)

end subroutine more_edges

!          ----------
subroutine more_verts(grid,errcode,size_wanted)
!          ----------

!----------------------------------------------------
! This routine increases the space for the number of vertices
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(out) :: errcode
integer, intent(in), optional :: size_wanted
!----------------------------------------------------
! Local variables:

type(vertex_t), pointer :: temp_vertex(:)
integer, pointer :: temp_vertex_type(:,:)
real(my_real), pointer :: temp_vertex_solution(:,:,:), &
                          temp_vertex_exact(:,:,:), temp_vertex_oldsoln(:,:,:)
integer :: astat, astat2, astat3, oldsize, newsize, i
!----------------------------------------------------
! Begin executable code

oldsize = size(grid%vertex)
if (present(size_wanted)) then
   newsize = size_wanted
else
   newsize = int(1.5*oldsize)
endif

errcode = 0
temp_vertex => grid%vertex
temp_vertex_type => grid%vertex_type
temp_vertex_solution => grid%vertex_solution
temp_vertex_exact => grid%vertex_exact
temp_vertex_oldsoln => grid%vertex_oldsoln
nullify(grid%vertex,grid%vertex_type)
allocate(grid%vertex(newsize),grid%vertex_type(newsize,grid%system_size), &
         grid%vertex_solution(newsize,grid%system_size,max(1,grid%num_eval)), &
         stat=astat)
if (grid%have_true) then
   allocate(grid%vertex_exact(newsize,grid%system_size,max(1,grid%num_eval)), &
            stat=astat2)
else
   astat2 = 0
endif
if (associated(grid%vertex_oldsoln)) then
   allocate(grid%vertex_oldsoln(newsize,grid%system_size,max(1,grid%num_eval)), &
            stat=astat3)
else
   astat3 = 0
endif
if (astat /= 0 .or. astat2 /= 0 .or. astat3 /= 0) then
   errcode = 1
   call fatal("increased allocation failed")
   grid%vertex => temp_vertex
   grid%vertex_type => temp_vertex_type
   grid%vertex_solution => temp_vertex_solution
   grid%vertex_exact => temp_vertex_exact
   grid%vertex_oldsoln => temp_vertex_oldsoln
   return
endif
grid%vertex(1:oldsize) = temp_vertex
grid%vertex_type(1:oldsize,:) = temp_vertex_type
grid%vertex_solution(1:oldsize,:,:) = temp_vertex_solution
if (grid%have_true) then
   grid%vertex_exact(1:oldsize,:,:) = temp_vertex_exact
endif
if (associated(grid%vertex_oldsoln)) then
   grid%vertex_oldsoln(1:oldsize,:,:) = temp_vertex_oldsoln
endif

deallocate(temp_vertex,temp_vertex_type,temp_vertex_solution,stat=astat)
if (associated(temp_vertex_exact)) deallocate(temp_vertex_exact,stat=astat)
if (associated(temp_vertex_oldsoln)) deallocate(temp_vertex_oldsoln,stat=astat)
grid%next_free_vert = oldsize+1

!call count_grid_memory(grid)

end subroutine more_verts

!          -----------------
subroutine count_grid_memory(grid)
!          -----------------

!----------------------------------------------------
! This routine adds up the amount of memory used in the grid data structure.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables:

integer :: isize, rsize, lsize, lssize, csize, p0size, p1size, p2size, &
           p3size, p4size, hsize
integer :: mem, i
!----------------------------------------------------
! Begin executable code

isize = bit_size(1)/8

! TEMP need a good portable way to determine the number of bytes in real.
! I think the following should always work, at least on the computers of 2006.
! The standard intrinsic function digits returns the number of binary digits in
! the real model for the type of real given.  The number of binary digits cannot
! possibly be bigger than the number of bits, so we look for the first bit
! size that can hold that number of binary digits.  I assume the number of
! bytes is a power of 2, it's at least 4 and no more than 64.

rsize = digits(1.0_my_real)
if (rsize <= 32) then
   rsize = 4
elseif (rsize <= 64) then
   rsize = 8
elseif (rsize <= 128) then
   rsize = 16
elseif (rsize <= 256) then
   rsize = 32
else
   rsize = 64
endif

! TEMP assume a default logical is 4 bytes and a small logical is 1 byte

lsize = 4
lssize = 1

! TEMP assume a default character is 1 byte

csize = 1

! The size of pointer descriptors varies among compilers.  There is also
! lost memory for small pointer arrays because of alignment requirements;
! this can be very large (lf95 can loose up to 192 bytes for every allocated
! pointer).  I do not include the memory lost to alignment in the memory
! count in this routine.  The following byte counts for pointer descriptors
! and alignment requirements were determined experimentally for some Fortran
! compilers on Linux.
!          Lahey   Intel   g95   Absoft   Nag
! Rank 0     8       4      4      24      4
! Rank 1    28      36     28      36     20
! Rank 2    44      48     40      48     32
! Rank 3    60      60     52      60     44
! Rank 4    76      72     64      72     56
! alignment 192     16     32      96      8

p0size = 8
p1size = 28
p2size = 44
p3size = 60
p4size = 76

! hash key size

hsize = KEY_SIZE*isize

mem = 0

! entries in grid_type

mem = mem + 3*p1size  & ! pointers for element, edge and vertex
          + 3*0       & ! hash tables ! TEMP look for size
          + 2*2*rsize & ! bounding box
          + 3*p3size  & ! solutions
          + 2*p2size  & ! error indicators
          + p1size    & ! pointer for eigenvalue
          + 2*rsize   & ! eigen linsys resids
          + 8*p1size  & ! pointers for eigenprob and errests
          + rsize     & ! max_blen
          + 2*p1size  & ! pointers for bp
          + 2*p2size  & ! pointers for type
          + p2size    & ! pointer for initial neighbor
          + 2*p1size  & ! pointer for linked list heads
          + 23*isize  & ! whole mess of integers
          + 2*lsize   & ! 2 logicals
          + FN_LEN*csize

! allocated pointers in grid_type, except element, edge and vertex

if (associated(grid%vertex_solution)) mem = mem + rsize*size(grid%vertex_solution)
if (associated(grid%vertex_exact)) mem = mem + rsize*size(grid%vertex_exact)
if (associated(grid%vertex_oldsoln)) mem = mem + rsize*size(grid%vertex_oldsoln)
if (associated(grid%element_errind)) mem = mem + rsize*size(grid%element_errind)
if (associated(grid%eigen_results%eigenvalue)) mem = mem + &
   rsize*size(grid%eigen_results%eigenvalue)
if (associated(grid%eigen_results%eigensolver_errbound)) mem = mem + &
   rsize*size(grid%eigen_results%eigensolver_errbound)
if (associated(grid%errest_energy)) mem = mem + rsize*size(grid%errest_energy)
if (associated(grid%errest_Linf)) mem = mem + rsize*size(grid%errest_Linf)
if (associated(grid%errest_L2)) mem = mem + rsize*size(grid%errest_L2)
if (associated(grid%errest_eigenvalue)) mem = mem + &
   rsize*size(grid%errest_eigenvalue)
if (associated(grid%bp_start)) mem = mem + rsize*size(grid%bp_start)
if (associated(grid%bp_finish)) mem = mem + rsize*size(grid%bp_finish)
if (associated(grid%edge_type)) mem = mem + isize*size(grid%edge_type)
if (associated(grid%vertex_type)) mem = mem + isize*size(grid%vertex_type)
if (associated(grid%initial_neighbor)) mem = mem + &
   isize*size(grid%initial_neighbor)
if (associated(grid%head_level_elem)) mem = mem + &
   isize*size(grid%head_level_elem)

! elements

if (associated(grid%element)) then
   mem = mem + size(grid%element)*( &
                  hsize*2                    & ! gid and mate
                + rsize                      & ! weight
                + p2size*3                   & ! pointers for solutions
                + p1size*2                   & ! pointers for error indicators
                + rsize*2                    & ! works
                + isize*VERTICES_PER_ELEMENT & ! vertex
                + isize*EDGES_PER_ELEMENT    & ! edge
                + isize*(4+MAX_CHILD+2)      & ! bunch of integers
                + lssize*5                   & ! small logicals
               )

   do i=1,grid%biggest_elem
      if (associated(grid%element(i)%solution)) mem = mem + rsize*size(grid%element(i)%solution)
      if (associated(grid%element(i)%exact)) mem = mem + rsize*size(grid%element(i)%exact)
      if (associated(grid%element(i)%oldsoln)) mem = mem + rsize*size(grid%element(i)%oldsoln)
   end do

endif

! edges

if (associated(grid%edge)) then
   mem = mem + size(grid%edge)*( &
                  hsize &    ! gid
                + isize*5 &  ! integers
                + p1size &   ! pointer for type
                + p2size*3 & ! pointers for solutions
                + isize    & ! next
               )

   do i=1,grid%biggest_edge
      if (associated(grid%edge(i)%solution)) mem = mem + rsize*size(grid%edge(i)%solution)
      if (associated(grid%edge(i)%exact)) mem = mem + rsize*size(grid%edge(i)%exact)
      if (associated(grid%edge(i)%oldsoln)) mem = mem + rsize*size(grid%edge(i)%oldsoln)
   end do

endif

! vertices

if (associated(grid%vertex)) then
   mem = mem + size(grid%vertex)*( &
                  hsize &    ! gid
                + rsize*2 &  ! coord
                + rsize &    ! bparam
                + p1size &   ! pointer for type
                + isize*4 &  ! integers
               )
endif

write(outunit,"(A,I12,A)") "Memory for grid data structure ",mem," bytes"

end subroutine count_grid_memory

!          -----------
subroutine phier2nodal(solution,xvert,yvert,degree)
!          -----------

!----------------------------------------------------
! This routine converts the values in solution from coefficients of a
! p-hierarchical basis to coefficients of a nodal basis, over four triangles
! that are siblings in an h refinement with the outer vertices given in
! (xvert,yvert).  The order of the coefficients (using cubics as an example) is

!                  vertex 1
!                      1
!                    / | \
!                   10 8  20
!                  / 12|22 \
!                 11   9    21
!                /     |      \
!    vertex 3   3-5-6--7-19-18-4  vertex 4
!                \     14     /
!                 16   |    24
!                  \ 17|25 /
!                   15 13 23
!                    \ | /
!                      2
!                  vertex 2

! For the p-hierarchical basis, those on the same edge or bubble functions of
! the same element are in order of degree, as in the basis function routines.
! All elements and edges must have the same degree.  The third and fourth
! triangles can be omitted by setting the 4th vertex to be huge(0.0_my_real).

! If the third and fourth triangles are omitted, the indices remain the same.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: solution(:,:)
real(my_real), intent(in) :: xvert(:),yvert(:)
integer, intent(in) :: degree
!----------------------------------------------------
! Local variables:

real(my_real) :: xnode(size(solution,dim=1)),ynode(size(solution,dim=1)), &
                 basis(((degree+1)*(degree+2))/2,size(solution,dim=1)), &
                 new_soln(size(solution,dim=1),size(solution,dim=2))
real(my_real) :: rdeg,xmid,ymid,fact,fact1,fact2,fact3
integer :: i,j,isub,nsub,nbasis,last
integer :: jsub(((degree+1)*(degree+2))/2), ksub(((degree+1)*(degree+2))/2)
!----------------------------------------------------
! Begin executable code


rdeg = real(degree,my_real)
nbasis = ((degree+1)*(degree+2))/2

! midpoint of the edge between vertices 1 and 2

xmid = (xvert(1)+xvert(2))/2
ymid = (yvert(1)+yvert(2))/2

! the coefficients of the vertex bases (1-4 and the midpoint 4+degree) are
! the same in both bases

new_soln = 0
new_soln(1:4,:) = solution(1:4,:)
new_soln(4+degree,:) = solution(4+degree,:)

! define the coordinates of the as yet unaddressed nodes in the first triangle,
! and set jsub to the index corresponding to those nodes,
! and set ksub to the index of the basis functions

ksub(1) = 1
ksub(2) = 3
ksub(3) = 4+degree
isub = 1
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xmid + (1-fact)*xvert(3)
   ynode(isub) = fact*ymid + (1-fact)*yvert(3)
   jsub(isub) = 4+i
   ksub(isub+3) = jsub(isub)
   isub = isub+1
end do
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xmid + (1-fact)*xvert(1)
   ynode(isub) = fact*ymid + (1-fact)*yvert(1)
   jsub(isub) = 4+degree+i
   ksub(isub+3) = jsub(isub)
   isub = isub+1
end do
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xvert(3) + (1-fact)*xvert(1)
   ynode(isub) = fact*yvert(3) + (1-fact)*yvert(1)
   jsub(isub) = 3+2*degree+i
   ksub(isub+3) = jsub(isub)
   isub = isub+1
end do
do i=1,degree-2
   fact3 = i/rdeg
   do j=1,degree-1-i
      fact2 = (1-fact3)*j/(rdeg-1)
      fact1 = 1-fact3-fact2
      xnode(isub) = fact1*xvert(1) + fact2*xvert(3) + fact3*xmid
      ynode(isub) = fact1*yvert(1) + fact2*yvert(3) + fact3*ymid
      jsub(isub) = 5+isub
      ksub(isub+3) = jsub(isub)
      isub = isub+1
   end do
end do

nsub = 3*(degree-1) + element_dof(degree)
if (nsub == 0) then
   last = 0
else
   last = jsub(nsub)
endif

! evaluate the p-hierarchical basis functions at the nodes

if (nsub /= 0) then
   call p_hier_basis_func(xnode(1:nsub),ynode(1:nsub), &
                          (/xvert(1),xvert(3),xmid/), &
                          (/yvert(1),yvert(3),ymid/), &
                          (/degree,degree,degree,degree/),"a",basis(:,1:nsub))
endif

! evaluate the solution at the nodes

do i=1,nsub
   new_soln(jsub(i),:) = 0
   do j=1,nbasis
      new_soln(jsub(i),:) = new_soln(jsub(i),:) + solution(ksub(j),:)*basis(j,i)
   end do
end do

! same process for the second triangle

ksub(1) = 2
ksub(2) = 3
ksub(3) = 4+degree
do i=1,degree-1
   ksub(i+3) = i+4
end do
isub = 1
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xmid + (1-fact)*xvert(2)
   ynode(isub) = fact*ymid + (1-fact)*yvert(2)
   jsub(isub) = last + isub
   ksub(isub+degree+2) = jsub(isub)
   isub = isub+1
end do
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xvert(3) + (1-fact)*xvert(2)
   ynode(isub) = fact*yvert(3) + (1-fact)*yvert(2)
   jsub(isub) = last + isub
   ksub(isub+degree+2) = jsub(isub)
   isub = isub+1
end do
do i=1,degree-2
   fact3 = i/rdeg
   do j=1,degree-1-i
      fact2 = (1-fact3)*j/(rdeg-1)
      fact1 = 1-fact3-fact2
      xnode(isub) = fact1*xvert(2) + fact2*xvert(3) + fact3*xmid
      ynode(isub) = fact1*yvert(2) + fact2*yvert(3) + fact3*ymid
      jsub(isub) = last + isub
      ksub(isub+degree+2) = jsub(isub)
      isub = isub+1
   end do
end do
nsub = 2*(degree-1) + element_dof(degree)
if (nsub == 0) then
   last = 0
else
   last = jsub(nsub)
endif
if (nsub /= 0) then
   call p_hier_basis_func(xnode(1:nsub),ynode(1:nsub), &
                          (/xvert(2),xvert(3),xmid/), &
                          (/yvert(2),yvert(3),ymid/), &
                          (/degree,degree,degree,degree/),"a",basis(:,1:nsub))
endif
do i=1,nsub
   new_soln(jsub(i),:) = 0
   do j=1,nbasis
      new_soln(jsub(i),:) = new_soln(jsub(i),:) + solution(ksub(j),:)*basis(j,i)
   end do
end do

! same process for the third triangle if it exists

if (xvert(4) /= huge(0.0_my_real)) then
   ksub(1) = 1
   ksub(2) = 4
   ksub(3) = 4+degree
   isub = 1
   do i=1,degree-1
      fact = i/rdeg
      xnode(isub) = fact*xmid + (1-fact)*xvert(4)
      ynode(isub) = fact*ymid + (1-fact)*yvert(4)
      jsub(isub) = last + isub
      ksub(isub+3) = jsub(isub)
      isub = isub+1
   end do
   do i=1,degree-1
      ksub(isub+2+i) = 4+degree+i
   end do
   do i=1,degree-1
      fact = i/rdeg
      xnode(isub) = fact*xvert(4) + (1-fact)*xvert(1)
      ynode(isub) = fact*yvert(4) + (1-fact)*yvert(1)
      jsub(isub) = last + isub
      ksub(isub+degree+2) = jsub(isub)
      isub = isub+1
   end do
   do i=1,degree-2
      fact3 = i/rdeg
      do j=1,degree-1-i
         fact2 = (1-fact3)*j/(rdeg-1)
         fact1 = 1-fact3-fact2
         xnode(isub) = fact1*xvert(1) + fact2*xvert(4) + fact3*xmid
         ynode(isub) = fact1*yvert(1) + fact2*yvert(4) + fact3*ymid
         jsub(isub) = last + isub
         ksub(isub+degree+2) = jsub(isub)
         isub = isub+1
      end do
   end do
   nsub = 2*(degree-1) + element_dof(degree)
   if (nsub == 0) then
      last = 0
   else
      last = jsub(nsub)
   endif
   if (nsub /= 0) then
      call p_hier_basis_func(xnode(1:nsub),ynode(1:nsub), &
                             (/xvert(1),xvert(4),xmid/), &
                             (/yvert(1),yvert(4),ymid/), &
                            (/degree,degree,degree,degree/),"a",basis(:,1:nsub))
   endif
   do i=1,nsub
      new_soln(jsub(i),:) = 0
      do j=1,nbasis
         new_soln(jsub(i),:) = new_soln(jsub(i),:) + solution(ksub(j),:)*basis(j,i)
      end do
   end do

! same process for the fourth triangle if it exists

   ksub(1) = 2
   ksub(2) = 4
   ksub(3) = 4+degree
   do i=1,degree-1
      ksub(i+3) = (degree+1)**2 + 1 + i
   end do
   do i=1,degree-1
      ksub(degree+2+i) = 2 + element_dof(degree)
   end do
   isub = 1
   do i=1,degree-1
      fact = i/rdeg
      xnode(isub) = fact*xvert(4) + (1-fact)*xvert(2)
      ynode(isub) = fact*yvert(4) + (1-fact)*yvert(2)
      jsub(isub) = last + isub
      ksub(isub+2*degree+1) = jsub(isub)
      isub = isub+1
   end do
   do i=1,degree-2
      fact3 = i/rdeg
      do j=1,degree-1-i
         fact2 = (1-fact3)*j/(rdeg-1)
         fact1 = 1-fact3-fact2
         xnode(isub) = fact1*xvert(2) + fact2*xvert(4) + fact3*xmid
         ynode(isub) = fact1*yvert(2) + fact2*yvert(4) + fact3*ymid
         jsub(isub) = last + isub
         ksub(isub+2*degree+1) = jsub(isub)
         isub = isub+1
      end do
   end do
   nsub = (degree-1) + element_dof(degree)
   if (nsub /= 0) then
      call p_hier_basis_func(xnode(1:nsub),ynode(1:nsub), &
                             (/xvert(2),xvert(4),xmid/), &
                             (/yvert(2),yvert(4),ymid/), &
                            (/degree,degree,degree,degree/),"a",basis(:,1:nsub))
   endif
   do i=1,nsub
      new_soln(jsub(i),:) = 0
      do j=1,nbasis
         new_soln(jsub(i),:) = new_soln(jsub(i),:) + solution(ksub(j),:)*basis(j,i)
      end do
   end do
endif

! copy new coefficients to return variable

solution = new_soln

end subroutine phier2nodal

!          -----------
subroutine nodal2phier(solution,xvert,yvert,degree)
!          -----------

!----------------------------------------------------
! This routine converts the values in solution from coefficients of a
! nodal basis to coefficients of a p-hierarchical basis, over two triangles
! that are mates with the outer vertices given in (xvert,yvert).  The order
! of the coefficients (using cubics as an example) is

!              vertex 1
!                  1
!                / | \
!               7  |  14
!              /   |   \
!             8    9    15
!            /     |     \
! vertex 3  3  11  |  16  4  vertex 4
!            \     |     /
!             6    10   13
!              \   |   /
!               5  |  12
!                \ | /
!                  2
!              vertex 2

! For the p-hierarchical basis, those on the same edge or bubble functions of
! the same element are in order of degree, as in the basis function routines.
! All elements and edges must have the same degree.  The second triangle
! can be omitted by setting the 4th vertex to be huge(0.0_my_real).

! If the second triangle is omitted, then all the indices after 3 are
! decremented by 1.

!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: solution(:,:)
real(my_real), intent(in) :: xvert(:),yvert(:)
integer, intent(in) :: degree
!----------------------------------------------------
! Local variables:

real(my_real) :: xnode((degree+1)**2),ynode((degree+1)**2), &
                 basis((degree+1)**2,(degree+1)**2), &
                 basistemp(((degree+1)*(degree+2))/2,((degree+1)*(degree+2))/2)
real(my_real) :: rdeg,fact,fact1,fact2,fact3
integer :: index1(((degree+1)*(degree+2))/2),index2(((degree+1)*(degree+2))/2),&
           ipiv((degree+1)**2)
integer :: i,j,isub,jsub,nnode,info
!----------------------------------------------------
! Begin executable code

rdeg = real(degree,my_real)

! define the nodes

! nodes that are vertices (will fix 4 later if only one triangle)

xnode(1:4) = xvert(1:4)
ynode(1:4) = yvert(1:4)
index1(1:3) = (/1,2,3/)
index2(1:3) = (/1,2,4/)

! nodes in first triangle

isub = 5
jsub = 4
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xvert(3) + (1-fact)*xvert(2)
   ynode(isub) = fact*yvert(3) + (1-fact)*yvert(2)
   index1(jsub) = isub
   jsub = jsub+1
   isub = isub+1
end do
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xvert(3) + (1-fact)*xvert(1)
   ynode(isub) = fact*yvert(3) + (1-fact)*yvert(1)
   index1(jsub) = isub
   jsub = jsub+1
   isub = isub+1
end do
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xvert(2) + (1-fact)*xvert(1)
   ynode(isub) = fact*yvert(2) + (1-fact)*yvert(1)
   index1(jsub) = isub
   jsub = jsub+1
   isub = isub+1
end do
do i=1,degree-2
   fact3 = i/rdeg
   do j=1,degree-1-i
      fact2 = (1-fact3)*j/(rdeg-1)
      fact1 = 1-fact3-fact2
      xnode(isub) = fact1*xvert(1) + fact2*xvert(2) + fact3*xvert(3)
      ynode(isub) = fact1*yvert(1) + fact2*yvert(2) + fact3*yvert(3)
      index1(jsub) = isub
      jsub = jsub+1
      isub = isub+1
   end do
end do

! if only one triangle, done, but remove vertex 4

if (xvert(4) == huge(0.0_my_real)) then
   nnode = isub-2
   xnode(4:nnode) = xnode(5:nnode+1)
   ynode(4:nnode) = ynode(5:nnode+1)
   index1(4:nnode) = index1(4:nnode)-1

! otherwise, nodes in second triangle

else

   jsub = 4
   do i=1,degree-1
      fact = i/rdeg
      xnode(isub) = fact*xvert(4) + (1-fact)*xvert(2)
      ynode(isub) = fact*yvert(4) + (1-fact)*yvert(2)
      index2(jsub) = isub
      jsub = jsub+1
      isub = isub+1
   end do
   do i=1,degree-1
      fact = i/rdeg
      xnode(isub) = fact*xvert(4) + (1-fact)*xvert(1)
      ynode(isub) = fact*yvert(4) + (1-fact)*yvert(1)
      index2(jsub) = isub
      jsub = jsub+1
      isub = isub+1
   end do
   do i=1,degree-1
      index2(jsub) = 2 + 2*degree + i
      jsub = jsub+1
   end do
   do i=1,degree-2
      fact3 = i/rdeg
      do j=1,degree-1-i
         fact2 = (1-fact3)*j/(rdeg-1)
         fact1 = 1-fact3-fact2
         xnode(isub) = fact1*xvert(1) + fact2*xvert(2) + fact3*xvert(4)
         ynode(isub) = fact1*yvert(1) + fact2*yvert(2) + fact3*yvert(4)
         index2(jsub) = isub
         jsub = jsub+1
         isub = isub+1
      end do
   end do
   nnode = isub-1

endif

! evaluate the p-hierarchical basis functions at the nodes

basis = 0
call p_hier_basis_func(xnode(index1),ynode(index1),xvert(1:3),yvert(1:3), &
                       (/degree,degree,degree,degree/),"a",basistemp)
basis(index1,index1) = basistemp
if (xvert(4) /= huge(0.0_my_real)) then
   call p_hier_basis_func(xnode(index2),ynode(index2), &
                          (/xvert(1),xvert(2),xvert(4)/), &
                          (/yvert(1),yvert(2),yvert(4)/), &
                          (/degree,degree,degree,degree/),"a",basistemp)
   basis(index2,index2) = basistemp
endif

! if n is the vector of nodal basis coefficients, p is the vector of
! p-hierarchical basis coefficients, and A is the matrix with A_ij having
! the value of the jth p-hierarchical basis at the ith node, then
! n = Ap.  Therefore, to get p from n, compute p = A^-1 n.  Note however that
! the matrix basis is A^T.

basis = transpose(basis)

if (my_real == kind(0.0)) then
   call sgesv(nnode,size(solution,dim=2),basis,size(basis,dim=1),ipiv, &
              solution,size(solution,dim=1),info)
else
   call dgesv(nnode,size(solution,dim=2),basis,size(basis,dim=1),ipiv, &
              solution,size(solution,dim=1),info)
endif

if (info /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("lapack failed in nodal2phier",intlist=(/info/))
   stop
endif

end subroutine nodal2phier

!          -----------
subroutine nodal2hhier(solution,xvert,yvert,degree)
!          -----------

!----------------------------------------------------
! This routine converts the values in solution from coefficients of a
! nodal basis to coefficients of a h-hierarchical basis, over four triangles
! that are siblings in an h refinement with the outer vertices given in
! (xvert,yvert).  The order of the coefficients is the same as in phier2nodal.
! All elements and edges must have the same degree.  The third and fourth
! triangles can be omitted by setting the 4th vertex to be huge(0.0_my_real).
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: solution(:,:)
real(my_real), intent(in) :: xvert(:),yvert(:)
integer, intent(in) :: degree
!----------------------------------------------------
! Local variables:

integer :: nred1, nblack1, nbasis1, nred2, nblack2, nbasis2, nred4, nblack4, &
           nbasis4, astat, i, j, k, rsub, bsub, ind
logical :: red
real(my_real) :: xvert5, yvert5, rdeg, w1, w2, w3
real(my_real), allocatable :: s(:,:),prod(:,:)
real(my_real), allocatable :: xnode(:), ynode(:), basis(:,:)
integer, allocatable :: ind_black(:), ind_red(:)
!----------------------------------------------------
! Begin executable code

! Allocate memory.  These used to be automatic arrays, but
! ifort (IFORT) 11.1 20100414 with OpenMP got a seg fault on s=0 when I used
! coef_decay with examples/parabolic and the degree hit 16 while computing
! the h-coarsen error indicator. (8/05/10)

allocate(s(size(solution,dim=1),size(solution,dim=1)), &
         prod(size(solution,dim=1),size(solution,dim=2)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in nodal2hhier")
   return
endif

! some useful constants

rdeg = real(degree,my_real)
xvert5 = (xvert(1)+xvert(2))/2
yvert5 = (yvert(1)+yvert(2))/2

! number of basis functions / nodes in 1 triangle and all 2 or 4 triangles

nbasis1 = ((degree+1)*(degree+2))/2
if (2*(degree/2) == degree) then
   nred1 = (degree/2)*(degree/2+1)
else
   nred1 = ((degree+1)/2)**2
endif
nblack1 = nbasis1 - nred1
nred2 = (degree*(degree+1))/2
nblack2 = nbasis1
nbasis2 = nred2 + nblack2
nred4 = degree**2
nblack4 = (degree+1)**2
nbasis4 = nred4 + nblack4
if (nbasis4 /= size(solution,dim=1)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("size mismatch in nodal2hhier")
   stop
endif

! allocate memory based on number of bases/nodes

allocate(xnode(nred2),ynode(nred2),basis(nblack2,nred2),ind_black(nblack2), &
         ind_red(nred2),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in nodal2hhier")
   return
endif

! The (i,j)th entry of the matrix that converts nodal coefficients to
! hierarchical coefficients is the negative of the jth h-hierarchical basis
! function at the ith node if j is black and i is red, and identity on the
! diagonal.  The black h-hierarchical bases are the nodal bases of the parent
! triangles.  Go through the 1 or 2 parents evaluating the nodal bases at
! the red nodes and inserting the value in the right place.

s = 0
do i=1,nbasis4
   s(i,i) = 1
end do

! parent

! define the red nodes, red node indices, and black basis indices; go along
! lines parallel to (vert1,vert2)

rsub = 1
bsub = 1
! vertex 1
ind_black(bsub) = 1
bsub = bsub + 1
! nodes along line between child 1 and child 3
red = .true.
ind = 5 + degree-1 + 1
do i=1,degree-1
   if (red) then
      w3 = i/rdeg
      xnode(rsub) = (1-w3)*xvert(1) + w3*xvert5
      ynode(rsub) = (1-w3)*yvert(1) + w3*yvert5
      ind_red(rsub) = ind
      rsub = rsub + 1
   else
      ind_black(bsub) = ind
      bsub = bsub + 1
   endif
   ind = ind + 1
   red = .not. red
end do
! central node
if (red) then
   xnode(rsub) = xvert5
   ynode(rsub) = yvert5
   ind_red(rsub) = degree + 4
   rsub = rsub + 1
else
   ind_black(bsub) = degree + 4
   bsub = bsub + 1
endif
red = .not. red
! nodes along edge between child 2 and child 4, in reverse order
ind = 5 + 4*(degree-1) + element_dof(degree)
do i=1,degree-1
   if (red) then
      w1 = i/rdeg
      xnode(rsub) = w1*xvert(2) + (1-w1)*xvert5
      ynode(rsub) = w1*yvert(2) + (1-w1)*yvert5
      ind_red(rsub) = ind
      rsub = rsub + 1
   else
      ind_black(bsub) = ind
      bsub = bsub + 1
   endif
   ind = ind - 1
   red = .not. red
end do
! vertex 2
ind_black(bsub) = 2
bsub = bsub + 1
! for each line parallel to (vert1, vert2)
do k=1,degree-1
   w2 = k/rdeg
! node on edge between vert1 and vert3
   ind_black(bsub) = 5 + 2*(degree-1) + k
   bsub = bsub + 1
! interior nodes in child 1
   red = .true.
   ind = 5 + 3*(degree-1) + k
   do j=1,degree-1-k
      w3 = (1-w2)*j/real(degree-k,my_real)
      w1 = 1-w2-w3
      if (red) then
         xnode(rsub) = w1*xvert(1) + w2*xvert(3) + w3*xvert5
         ynode(rsub) = w1*yvert(1) + w2*yvert(3) + w3*yvert5
         ind_red(rsub) = ind
         rsub = rsub + 1
      else
         ind_black(bsub) = ind
         bsub = bsub + 1
      endif
      ind = ind + degree - j - 1
      red = .not. red
   end do
! line between child 1 and child 2
   if (red) then
      xnode(rsub) = w2*xvert(3) + (1-w2)*xvert5
      ynode(rsub) = w2*yvert(3) + (1-w2)*yvert5
      ind_red(rsub) = 5 + degree-1 - k
      rsub = rsub + 1
   else
      ind_black(bsub) = 5 + degree-1 - k
      bsub = bsub + 1
   endif
   red = .not. red
! interior nodes in child 2
   ind = 5 + 5*(degree-1) + 2*element_dof(degree) - ((k-1)*k)/2
   do j=1,degree-1-k
      w1 = (1-w2)*j/real(degree-k,my_real)
      w3 = 1-w2-w1
      if (red) then
         xnode(rsub) = w1*xvert(2) + w2*xvert(3) + w3*xvert5
         ynode(rsub) = w1*yvert(2) + w2*yvert(3) + w3*yvert5
         ind_red(rsub) = ind
         rsub = rsub + 1
      else
         ind_black(bsub) = ind
         bsub = bsub + 1
      endif
      ind = ind - (k+j)
      red = .not. red
   end do
! node on line between vert2 and vert3
   ind_black(bsub) = 5 + 4*(degree-1) + element_dof(degree) + k
   bsub = bsub + 1
end do ! lines parallel to (vert1, vert2)
! vertex 3
ind_black(bsub) = 3
bsub = bsub + 1

! evaluate the nodal basis functions of the parent at the red nodes

call nodal_basis_func(xnode,ynode,(/xvert(1),xvert(2),xvert(3)/), &
                      (/yvert(1),yvert(2),yvert(3)/),degree,"a",basis)

! copy the values into the matrix

do j=1,nblack2
   do i=1,nred2
      s(ind_red(i),ind_black(j)) = -basis(j,i)
   end do
end do

if (xvert(4) /= huge(0.0_my_real)) then

! mate

! define the red nodes, red node indices, and black basis indices; go along
! lines parallel to (vert1,vert2)

   rsub = 1
   bsub = 1
! vertex 1
   ind_black(bsub) = 1
   bsub = bsub + 1
! nodes along line between child 1 and child 3
   red = .true.
   ind = 5 + degree-1 + 1
   do i=1,degree-1
      if (red) then
         w3 = i/rdeg
         xnode(rsub) = (1-w3)*xvert(1) + w3*xvert5
         ynode(rsub) = (1-w3)*yvert(1) + w3*yvert5
         ind_red(rsub) = ind
         rsub = rsub + 1
      else
         ind_black(bsub) = ind
         bsub = bsub + 1
      endif
      ind = ind + 1
      red = .not. red
   end do
! central node
   if (red) then
      xnode(rsub) = xvert5
      ynode(rsub) = yvert5
      ind_red(rsub) = degree + 4
      rsub = rsub + 1
   else
      ind_black(bsub) = degree + 4
      bsub = bsub + 1
   endif
   red = .not. red
! nodes along edge between child 2 and child 4, in reverse order
   ind = 5 + 4*(degree-1) + element_dof(degree)
   do i=1,degree-1
      if (red) then
         w1 = i/rdeg
         xnode(rsub) = w1*xvert(2) + (1-w1)*xvert5
         ynode(rsub) = w1*yvert(2) + (1-w1)*yvert5
         ind_red(rsub) = ind
         rsub = rsub + 1
      else
         ind_black(bsub) = ind
         bsub = bsub + 1
      endif
      ind = ind - 1
      red = .not. red
   end do
! vertex 2
   ind_black(bsub) = 2
   bsub = bsub + 1
! for each line parallel to (vert1, vert2)
   do k=1,degree-1
      w2 = k/rdeg
! node on edge between vert1 and vert4
      ind_black(bsub) = 5 + 6*(degree-1) + 2*element_dof(degree) + k
      bsub = bsub + 1
! interior nodes in child 3
      red = .true.
      ind = 5 + 7*(degree-1) + 2*element_dof(degree) + k
      do j=1,degree-1-k
         w3 = (1-w2)*j/real(degree-k,my_real)
         w1 = 1-w2-w3
         if (red) then
            xnode(rsub) = w1*xvert(1) + w2*xvert(4) + w3*xvert5
            ynode(rsub) = w1*yvert(1) + w2*yvert(4) + w3*yvert5
            ind_red(rsub) = ind
            rsub = rsub + 1
         else
            ind_black(bsub) = ind
            bsub = bsub + 1
         endif
         ind = ind + degree - j - 1
         red = .not. red
      end do
! line between child 3 and child 4
      if (red) then
         xnode(rsub) = w2*xvert(4) + (1-w2)*xvert5
         ynode(rsub) = w2*yvert(4) + (1-w2)*yvert5
         ind_red(rsub) = 6 + 6*(degree-1) + 2*element_dof(degree) - k
         rsub = rsub + 1
      else
         ind_black(bsub) = 6 + 6*(degree-1) + 2*element_dof(degree) - k
         bsub = bsub + 1
      endif
      red = .not. red
! interior nodes in child 4
      ind = 5 + 8*(degree-1) + 4*element_dof(degree) - ((k-1)*k)/2
      do j=1,degree-1-k
         w1 = (1-w2)*j/real(degree-k,my_real)
         w3 = 1-w2-w1
         if (red) then
            xnode(rsub) = w1*xvert(2) + w2*xvert(4) + w3*xvert5
            ynode(rsub) = w1*yvert(2) + w2*yvert(4) + w3*yvert5
            ind_red(rsub) = ind
            rsub = rsub + 1
         else
            ind_black(bsub) = ind
            bsub = bsub + 1
         endif
         ind = ind - (k+j)
         red = .not. red
      end do
! node on line between vert2 and vert4
      ind_black(bsub) = 5 + 7*(degree-1) + 3*element_dof(degree) + k
      bsub = bsub + 1
   end do ! lines parallel to (vert1, vert2)
! vertex 4
   ind_black(bsub) = 4
   bsub = bsub + 1

! evaluate the nodal basis functions of the mate at the red nodes

   call nodal_basis_func(xnode,ynode,(/xvert(1),xvert(2),xvert(4)/), &
                         (/yvert(1),yvert(2),yvert(4)/),degree,"a",basis)

! copy the values into the matrix

   do j=1,nblack2
      do i=1,nred2
         s(ind_red(i),ind_black(j)) = -basis(j,i)
      end do
   end do

endif ! xvert(4) /= huge

! convert values in solution to h-hierarchical

if (my_real == kind(1.0)) then
   call sgemm("N","N",nbasis4,size(solution,dim=2),nbasis4,1.0, &
              s,nbasis4,solution,nbasis4,0.0,prod,nbasis4)
elseif (my_real == kind(1.0d0)) then
   call dgemm("N","N",nbasis4,size(solution,dim=2),nbasis4,1.0d0, &
              s,nbasis4,solution,nbasis4,0.0d0,prod,nbasis4)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither single nor double precision. Can't call GEMM")
   stop
endif

! copy result back to solution

solution = prod

deallocate(xnode,ynode,basis,ind_black,ind_red)
deallocate(s,prod)

end subroutine nodal2hhier

!        ----------------
function element_diameter(grid,elem)
!        ----------------

!----------------------------------------------------
! This routine computes the diameter of element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real) :: element_diameter
!----------------------------------------------------
! Local variables:

integer :: i, v1, v2
!----------------------------------------------------
! Begin executable code

! take the diameter to be the longest side length

element_diameter = -1
do i=1,EDGES_PER_ELEMENT
   v1 = grid%edge(grid%element(elem)%edge(i))%vertex(1)
   v2 = grid%edge(grid%element(elem)%edge(i))%vertex(2)

   element_diameter = max(element_diameter, &
      sqrt((grid%vertex(v1)%coord%x - grid%vertex(v2)%coord%x)**2 + &
           (grid%vertex(v1)%coord%y - grid%vertex(v2)%coord%y)**2))
end do

end function element_diameter

!        --------------
function element_volume(elem,grid)
!        --------------

!----------------------------------------------------
! This routine computes the area of triangular element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
type(grid_type), intent(in) :: grid
real(my_real) :: element_volume
!----------------------------------------------------
! Local variables:

real(my_real) :: A(2,2)
integer :: i
!----------------------------------------------------
! Begin executable code

do i=1,2
   A(1,i) = grid%vertex(grid%element(elem)%vertex(i))%coord%x - &
            grid%vertex(grid%element(elem)%vertex(3))%coord%x
   A(2,i) = grid%vertex(grid%element(elem)%vertex(i))%coord%y - &
            grid%vertex(grid%element(elem)%vertex(3))%coord%y
end do

element_volume = abs(A(1,1)*A(2,2) - A(2,1)*A(1,2))/2

end function element_volume

!          -------------
subroutine list_elements(grid,list,nelem,level,own,leaf,unowned_neigh, &
                         bound_vert,owned_neigh)
!          -------------

!----------------------------------------------------
! This routine creates a list of elements in the beginning of list.
! nelem is the number of returned elements.
! If level is present, it is the h-refinement level to use; otherwise all levels
! If own is present and true, only owned elements are listed.
! If leaf is present and true, only leaves are listed.
! If unowned_neigh is present and true, only elements with an unowned neighbor
! are listed.
! If bound_vert is present and true, only elements that have a vertex that
! touches the partition boundary are listed.  Currently this requires
! own=leaf=true and unowned_neigh=owned_neigh=false.
! If owned_neigh is present and true, unowned elements that have an owned
! neighbor are included (despite own being true).
! list must be large enough to hold the list.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(out) :: list(:), nelem
integer, intent(in), optional :: level
logical, intent(in), optional :: own, leaf, unowned_neigh, bound_vert, &
                                 owned_neigh
!----------------------------------------------------
! Local variables:

integer :: lolim, hilim, lev, elem, next_elem
integer :: loc_level, i, neigh(NEIGHBORS_PER_ELEMENT), astat
logical :: loc_own, loc_leaf, loc_unowned_neigh, unowned_neigh_satisfied, &
           loc_bound_vert, loc_owned_neigh, owned_neigh_satisfied
!----------------------------------------------------
! Begin executable code

if (present(own)) then
   loc_own = own
else
   loc_own = .false.
endif

if (present(leaf)) then
   loc_leaf = leaf
else
   loc_leaf = .false.
endif

if (present(unowned_neigh)) then
   loc_unowned_neigh = unowned_neigh
else
   loc_unowned_neigh = .false.
endif

if (present(bound_vert)) then
   loc_bound_vert = bound_vert
else
   loc_bound_vert = .false.
endif

if (present(owned_neigh)) then
   loc_owned_neigh = owned_neigh
else
   loc_owned_neigh = .false.
endif

if (present(level)) then
   lolim = level
   hilim = level
else
   lolim = 1
   hilim = grid%nlev
endif

! I had a hard enough time getting the logic right before I added the bound_vert
! option, so for now I'm just going to handle that one separately.

if (loc_bound_vert) then
   call bound_vert_list
   return
endif

nelem = 0
do lev=lolim,hilim
   next_elem = grid%head_level_elem(lev)
   do while (next_elem /= END_OF_LIST)
      elem = next_elem
      next_elem = grid%element(elem)%next
      if (loc_leaf .and. .not. grid%element(elem)%isleaf) cycle
      if (loc_own .and. .not. grid%element(elem)%iown .and. &
          (.not. loc_owned_neigh)) cycle
      if (loc_unowned_neigh .or. loc_owned_neigh) then
         neigh = get_neighbors(elem,grid)
      endif
      if (loc_unowned_neigh) then
         unowned_neigh_satisfied = .false.
         do i=1,NEIGHBORS_PER_ELEMENT
            if (neigh(i) /= BOUNDARY) then
               if (.not. grid%element(neigh(i))%iown) then
                  unowned_neigh_satisfied = .true.
               endif
            endif
         end do
      else
         unowned_neigh_satisfied = .false.
      endif
      if (loc_owned_neigh) then
         owned_neigh_satisfied = .false.
         do i=1,NEIGHBORS_PER_ELEMENT
            if (neigh(i) /= BOUNDARY) then
               if ((.not. grid%element(elem)%iown) .and. &
                   grid%element(neigh(i))%iown) then
                  owned_neigh_satisfied = .true.
               endif
            endif
         end do
      else
         owned_neigh_satisfied = .false.
      endif
      if (loc_unowned_neigh .and. .not. unowned_neigh_satisfied) cycle
      if (loc_own .and. .not. grid%element(elem)%iown .and. &
          loc_owned_neigh .and. .not. owned_neigh_satisfied) cycle
      nelem = nelem + 1
      list(nelem) = elem
   end do
end do

contains

subroutine bound_vert_list
logical, allocatable :: isboundary(:)

if (.not. loc_own .or. .not. loc_leaf .or. loc_unowned_neigh .or. loc_owned_neigh) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("list_elements: bound_vert=T requires own=T, leaf=T, unowned_neigh=T, owned_neigh=T")
   stop
endif

! mark each vertex as being on the partition boundary or not by going through
! the elements and marking each vertex as being on the boundary if I own the
! element or vertex but not the other

allocate(isboundary(grid%biggest_vert),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in list_elements")
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

! make a list of owned leaves that touch the partition boundary by going
! through the elements and finding those that have a marked vertex.

nelem = 0
do lev=lolim,hilim
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
         do i=1,VERTICES_PER_ELEMENT
            if (isboundary(grid%element(elem)%vertex(i))) then
               nelem = nelem + 1
               list(nelem) = elem
               exit
            endif
         end do
      endif
      elem = grid%element(elem)%next
   end do
end do

deallocate(isboundary)

end subroutine bound_vert_list

end subroutine list_elements

!          -----------------------
subroutine list_edges_without_rule(grid,refine_control,element_list,nelem, &
                                   edge_list,nedge)
!          -----------------------

!----------------------------------------------------
! This routine makes a list of the edges of elements in element_list(1:nelem)
! that do not satisfy the edge rule.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: element_list(:), nelem
integer, intent(out) :: edge_list(:), nedge
!----------------------------------------------------
! Local variables:

integer :: i, j, k, elem, edge, need_deg
integer, pointer :: edge_elements(:)
logical(small_logical) :: visited_edge(grid%biggest_edge)
!----------------------------------------------------
! Begin executable code

visited_edge = .false.

nedge = 0
do i=1,nelem
   elem = element_list(i)
   do j=1,EDGES_PER_ELEMENT
      edge = grid%element(elem)%edge(j)
      if (visited_edge(edge)) cycle
      visited_edge(edge) = .true.
      call get_edge_elements(grid,edge,edge_elements)
      if (refine_control%edge_rule == MINIMUM_RULE) then
         need_deg = huge(0)
      else ! (refine_control%edge_rule == MAXIMUM_RULE)
         need_deg = 0
      endif
      do k=1,size(edge_elements)
         if (refine_control%edge_rule == MINIMUM_RULE) then
            need_deg = min(grid%element(edge_elements(k))%degree,need_deg)
         else ! (refine_control%edge_rule == MAXIMUM_RULE)
            need_deg = max(grid%element(edge_elements(k))%degree,need_deg)
         endif
      end do
      if (grid%edge(edge)%degree /= need_deg) then
         nedge = nedge + 1
         edge_list(nedge) = edge
      endif
      deallocate(edge_elements)
   end do
end do

end subroutine list_edges_without_rule

!        -----------
function element_dof(deg)
!        -----------

!----------------------------------------------------
! This routine returns the number of degrees of freedom on the interior of
! an element of degree deg
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: deg
integer :: element_dof
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

element_dof = ((deg-1)*(deg-2))/2

end function element_dof

!        --------
function face_dof(deg)
!        --------

!----------------------------------------------------
! This routine returns the number of degrees of freedom on a face, which is
! 0 because there are no faces in 2D
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: deg
integer :: face_dof
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

face_dof = 0

end function face_dof

!        ---------
function count_dof(grid,procs,just_owned)
!        ---------

!----------------------------------------------------
! This routine counts the degrees of freedom in grid.  If procs is present,
! it returns the dof of the global grid.  If just_owned is present and true, it
! returns the owned dof of this processor's grid.  Otherwise it returns all
! the dof of this processor's grid.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in), optional :: procs
logical, intent(in), optional :: just_owned
integer :: count_dof
!----------------------------------------------------
! Local variables:

logical(small_logical) :: visited_vert(grid%biggest_vert), &
                          visited_edge(grid%biggest_edge)
logical :: only_owned
integer :: lev, ivert, vert, elem, deg, iedge, edge
!----------------------------------------------------
! Begin executable code

! can't have both procs and just_owned

if (present(procs) .and. present(just_owned)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Can't give both procs and just_owned in count_dof.")
   stop
endif

! determine if the count is of all dof or only the owned dof

if (present(procs)) then
   only_owned = .true.
elseif (present(just_owned)) then
   only_owned = just_owned
else
   only_owned = .false.
endif

! count the vertices; each has system_size dofs

count_dof = 0
visited_vert = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do ivert=1,VERTICES_PER_ELEMENT
         vert = grid%element(elem)%vertex(ivert)
         if (visited_vert(vert)) cycle
         visited_vert(vert) = .true.
         if (grid%element(grid%vertex(vert)%assoc_elem)%iown .or. &
             .not. only_owned) then
            count_dof = count_dof + count(grid%vertex_type(vert,:)/=PERIODIC_SLAVE)
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

! go through the leaf elements counting the degree-dependent dofs and counting
! the dofs on the edges that have not yet been visited

visited_edge = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf .and. &
          (grid%element(elem)%iown .or. .not. only_owned)) then
         deg = grid%element(elem)%degree
         if (deg >= 3) then
            count_dof = count_dof + grid%system_size*element_dof(deg)
         endif
         do iedge=1,EDGES_PER_ELEMENT
            edge = grid%element(elem)%edge(iedge)
            if (visited_edge(edge)) cycle
            visited_edge(edge) = .true.
            if (grid%element(grid%edge(edge)%assoc_elem)%iown .or. &
                 .not. only_owned) then
               count_dof = count_dof + (grid%edge(edge)%degree-1)* &
                           count(grid%edge_type(edge,:)/=PERIODIC_SLAVE)
            endif
         end do
      endif
      elem = grid%element(elem)%next
   end do
end do

! sum over all processors

if (present(procs)) then
   count_dof = phaml_global_sum(procs,count_dof,3501)
endif

end function count_dof

!        -------------------------
function compute_global_max_errind(grid,procs,still_sequential,elem_list,nelem)
!        -------------------------

!----------------------------------------------------
! This routine computes the global maximum error indicator
! If elem_list and nelem are present, it computes it OpenMP parallel.
! elem_list must contain at least all owned leaf elements, but can contain more.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
real(my_real) :: compute_global_max_errind
integer, intent(in), optional :: elem_list(:)
integer, intent(in), optional :: nelem
!----------------------------------------------------
! Local variables:

integer :: lev, elem, i
real(my_real) :: max_errind
!----------------------------------------------------
! Begin executable code

! OpenMP version

if (present(elem_list)) then
   if (.not. present(nelem)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("compute_global_max_errind must contain neither or both of elem_list and nelem")
      stop
   endif

   max_errind = 0.0_my_real

!$omp parallel do &
!$omp  default(shared) &
!$omp  private(i,elem) &
!$omp  reduction(max : max_errind)

   do i=1,nelem
      elem = elem_list(i)
      if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
         max_errind = max(max_errind, &
                 maxval(grid%element_errind(elem,:))/grid%element(elem)%work)
      endif
   end do
!$omp end parallel do

   compute_global_max_errind = max_errind

else

! sequential version

   compute_global_max_errind = 0.0_my_real
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
            compute_global_max_errind = max(compute_global_max_errind, &
                    maxval(grid%element_errind(elem,:))/grid%element(elem)%work)
         endif
         elem = grid%element(elem)%next
      end do
   end do

endif

if (.not. still_sequential) then
   compute_global_max_errind = phaml_global_max(procs, &
                                                compute_global_max_errind,3540)
endif

end function compute_global_max_errind

!        -------------
function zcoord_scalar(pt)
!        -------------

!----------------------------------------------------
! This routine returns 0.  In grid_util3D.f90 it returns the z coordinate of pt
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(point), intent(in) :: pt
real(my_real) :: zcoord_scalar

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

zcoord_scalar = 0.0_my_real

end function zcoord_scalar

!        ------------
function zcoord_array(pt)
!        ------------

!----------------------------------------------------
! This routine returns 0.  In grid_util3D.f90 it returns the z coordinate of pt
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(point), intent(in) :: pt(:)
real(my_real) :: zcoord_array(size(pt))

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

zcoord_array = 0.0_my_real

end function zcoord_array

!          ----------
subroutine set_zcoord(pt,val)
!          ----------

!----------------------------------------------------
! This routine does nothing; the 3D counterpart sets the z coordinate of pt.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(point), intent(in) :: pt
real(my_real), intent(in) :: val
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

end subroutine set_zcoord

!----------------------------------------------------
! The following routines act as wrappers for the user supplied routines with
! point coordinate arguments, so the user doesn't have to make z an optional
! argument.
!----------------------------------------------------

subroutine my_pdecoefs(elem,x,y,z,cxx,cyy,czz,cxy,cxz,cyz,cx,cy,cz,c,rs)
use global
integer, intent(in) :: elem
real(my_real), intent(in) :: x,y,z
real(my_real), intent(out) :: cxx(:,:),cyy(:,:),czz(:,:),cxy(:,:),cxz(:,:), &
                              cyz(:,:),cx(:,:),cy(:,:),cz(:,:),c(:,:),rs(:)
secret_elem = elem
call pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
czz=0; cxz=0; cyz=0; cz=0
end subroutine my_pdecoefs

subroutine my_bconds(x,y,z,bmark,itype,c,rs)
use global
real(my_real), intent(in) :: x,y,z
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:),rs(:)
call bconds(x,y,bmark,itype,c,rs)
end subroutine my_bconds

function my_iconds(x,y,z,comp,eigen)
use global
real(my_real), intent(in) :: x,y,z
integer, intent(in) :: comp, eigen
real(my_real) :: my_iconds
my_iconds = iconds(x,y,comp,eigen)
end function my_iconds

function my_trues(x,y,z,comp,eigen)
use global
real (my_real), intent(in) :: x,y,z
integer, intent(in) :: comp, eigen
real (my_real) :: my_trues
my_trues = trues(x,y,comp,eigen)
end function my_trues

function my_truexs(x,y,z,comp,eigen)
use global
real (my_real), intent(in) :: x,y,z
integer, intent(in) :: comp, eigen
real (my_real) :: my_truexs
my_truexs = truexs(x,y,comp,eigen)
end function my_truexs

function my_trueys(x,y,z,comp,eigen)
use global
real (my_real), intent(in) :: x,y,z
integer, intent(in) :: comp, eigen
real (my_real) :: my_trueys
my_trueys = trueys(x,y,comp,eigen)
end function my_trueys

function my_truezs(x,y,z,comp,eigen)
use global
real (my_real), intent(in) :: x,y,z
integer, intent(in) :: comp, eigen
real (my_real) :: my_truezs
my_truezs = 0.0_my_real
end function my_truezs

subroutine my_boundary_point(ipiece,s,x,y)
use global
integer, intent(in) :: ipiece
real(my_real), intent(in) :: s
real(my_real), intent(out) :: x,y
call boundary_point(ipiece,s,x,y)
end subroutine my_boundary_point

function my_boundary_npiece(hole)
integer, intent(in) :: hole
integer :: my_boundary_npiece
my_boundary_npiece = boundary_npiece(hole)
end function my_boundary_npiece

function my_phaml_integral_kernel(kernel,x,y,z)
use global
integer, intent(in) :: kernel
real(my_real), intent(in) :: x,y,z
real(my_real) :: my_phaml_integral_kernel
my_phaml_integral_kernel = phaml_integral_kernel(kernel,x,y)
end function my_phaml_integral_kernel

function my_regularity(x,y,z)
use global
real(my_real), intent(in) :: x(*),y(*),z(*)
real(my_real) :: my_regularity
my_regularity = regularity(x,y)
end function my_regularity

end module grid_util
