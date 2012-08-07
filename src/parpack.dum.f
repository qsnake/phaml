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

! This file contains dummy routines to satisfy PARPACK external references
! when PARPACK is not needed.  The code was lifted from the PARPACK source code.
! These call the equivalent nonparallel routines from ARPACK, because it
! appears SLEPc calls these parallel routines even when built for sequential.

      subroutine psnaupd
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &     iparam, ipntr, workd, workl, lworkl, info )
      use message_passing
      integer    comm
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Real 
     &           tol
      integer    iparam(11), ipntr(14)
      Real 
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call snaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &     iparam, ipntr, workd, workl, lworkl, info )
      end


      subroutine pdnaupd 
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &     iparam, ipntr, workd, workl, lworkl, info )
      use message_passing
      integer    comm
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision 
     &           tol
      integer    iparam(11), ipntr(14)
      Double precision 
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call dnaupd 
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &     iparam, ipntr, workd, workl, lworkl, info )
      end

      subroutine psneupd 
     &         (comm , rvec , howmny, select, dr    , di  ,
     &          z    , ldz  , sigmar, sigmai, workev, bmat,
     &          n    , which, nev   , tol   , resid ,
     &          ncv  , v    , ldv   , iparam, ipntr ,
     &          workd, workl, lworkl, info  )
      use message_passing
      integer   comm
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Real      
     &           sigmar, sigmai, tol
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Real 
     &           dr(nev+1)    , di(nev+1)    , resid(n)  ,
     &           v(ldv,ncv)   , z(ldz,*)     , workd(3*n),
     &           workl(lworkl), workev(3*ncv)
      call sneupd 
     &         (rvec , howmny, select, dr    , di  ,
     &          z    , ldz  , sigmar, sigmai, workev, bmat,
     &          n    , which, nev   , tol   , resid ,
     &          ncv  , v    , ldv   , iparam, ipntr ,
     &          workd, workl, lworkl, info  )
      end

      subroutine pdneupd  
     &         (comm , rvec , howmny, select, dr    , di  ,
     &          z    , ldz  , sigmar, sigmai, workev, bmat,
     &          n    , which, nev   , tol   , resid ,
     &          ncv  , v    , ldv   , iparam, ipntr ,
     &          workd, workl, lworkl, info  )
      use message_passing
      integer   comm
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision      
     &           sigmar, sigmai, tol
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision 
     &           dr(nev+1)    , di(nev+1)    , resid(n)  ,
     &           v(ldv,ncv)   , z(ldz,*)     , workd(3*n),
     &           workl(lworkl), workev(3*ncv)
      call dneupd  
     &         (rvec , howmny, select, dr    , di  ,
     &          z    , ldz  , sigmar, sigmai, workev, bmat,
     &          n    , which, nev   , tol   , resid ,
     &          ncv  , v    , ldv   , iparam, ipntr ,
     &          workd, workl, lworkl, info  )
      end

      subroutine pdsaupd
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &     iparam, ipntr, workd, workl, lworkl, info )
c
      use message_passing
      character  bmat*1, which*2
      integer    comm, ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
      integer    iparam(11), ipntr(11)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call dsaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &     iparam, ipntr, workd, workl, lworkl, info )
      end

      subroutine pdseupd
     &    (comm, rvec, howmny, select, d, z, ldz, sigma, bmat,
     &     n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat, howmny, which*2
      logical    rvec
      integer    comm, info, ldz, ldv, lworkl, n, ncv, nev
      Double precision
     &           sigma, tol
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Double precision
     &           d(nev), resid(n), v(ldv,ncv), z(ldz, nev),
     &           workd(2*n), workl(lworkl)
      call dseupd
     &    (rvec, howmny, select, d, z, ldz, sigma, bmat,
     &     n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      end

      subroutine pssaupd
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &     iparam, ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat*1, which*2
      integer    comm, ido, info, ldv, lworkl, n, ncv, nev
      Real
     &           tol
      integer    iparam(11), ipntr(11)
      Real
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      call ssaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &     iparam, ipntr, workd, workl, lworkl, info )
      end

      subroutine psseupd
     &    (comm, rvec, howmny, select, d, z, ldz, sigma, bmat,
     &     n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      use message_passing
      character  bmat, howmny, which*2
      logical    rvec
      integer    comm, info, ldz, ldv, lworkl, n, ncv, nev
      Real
     &           sigma, tol
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Real
     &           d(nev), resid(n), v(ldv,ncv), z(ldz, nev),
     &           workd(2*n), workl(lworkl)
      call sseupd
     &    (rvec, howmny, select, d, z, ldz, sigma, bmat,
     &     n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      end

      subroutine pznaupd
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &     iparam, ipntr, workd, workl, lworkl, rwork, info )
      use message_passing
      integer    comm
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
      integer    iparam(11), ipntr(14)
      Complex*16
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      Double precision
     &           rwork(ncv)
      call znaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &     iparam, ipntr, workd, workl, lworkl, rwork, info )
      end

      subroutine pzneupd
     &         ( comm , rvec  , howmny, select, d    ,
     &           z    , ldz   , sigma , workev, bmat ,
     &           n    , which , nev   , tol   , resid,
     &           ncv  , v     , ldv   , iparam, ipntr,
     &           workd, workl , lworkl, rwork , info )
      use message_passing
      integer   comm
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Complex*16
     &           sigma
      Double precision
     &           tol
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision
     &           rwork(ncv)
      Complex*16
     &           d(nev)     , resid(n)  , v(ldv,ncv)   ,
     &           z(ldz, nev), workd(3*n), workl(lworkl),
     &           workev(2*ncv)
      call zneupd
     &         ( rvec  , howmny, select, d    ,
     &           z    , ldz   , sigma , workev, bmat ,
     &           n    , which , nev   , tol   , resid,
     &           ncv  , v     , ldv   , iparam, ipntr,
     &           workd, workl , lworkl, rwork , info )
      end
