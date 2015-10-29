SUBROUTINE cbcg_solve(h_psi, cg_psi, e, d0psi, dpsi, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol, cw, tprec, itol)
!
!-----------------------------------------------------------------------
!
!   Iterative solution of the linear system:
!
!                 ( h - e + w + i * eta ) * x = b
!                 ( h - cw + i * eta ) * G(G,G') = -\delta(G,G')
!
!   where h is a complex hermitian matrix, e, w, and eta are
!   real scalar, x and b are complex vectors

USE kinds,       ONLY: DP
!USE control_gw,  ONLY: maxter_green
USE mp_global,      ONLY : intra_pool_comm
USE mp,             ONLY : mp_sum
USE io_global,            ONLY : stdout, ionode

implicit none

! first I/O variables

logical :: tprec

integer ::   ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             kter, & ! output: counter on iterations
             nbnd, & ! input: the number of bands
             npol, & ! input: number of components of the wavefunctions
             ik      ! input: the k point

real(DP) :: &
             anorm,   & ! output: the norm of the error in the solution
             ethr,    & ! input: the required precision
             h_diag(ndmx*npol,nbnd) ! input: an estimate of ( H - \epsilon )
             !h_diag(ndmx*npol,nbnd) ! input: an estimate of ( H - \epsilon )

!COMPLEX(DP) :: h_diag(ndmx*npol,nbnd) ! input: an estimate of ( H - \epsilon )

  complex(DP) :: &
             dpsi    (ndmx*npol, nbnd), & ! output: the solution of the linear syst
             d0psi   (ndmx*npol, nbnd)    ! input: the known term
!             dpsitil (ndmx*npol, nbnd)   ! output: the conjugate solution!

  logical :: conv_root ! output: if true the root is converged

  external h_psi       ! input: the routine computing h_psi

  external cg_psi      ! input: the routine computing cg_psi

  !
  !  here the local variables
  !

!HL upping iterations to get convergence with green_linsys?

  integer, parameter :: maxiter = 200
  !integer, parameter :: maxter = 600
  !the maximum number of iterations
  integer :: iter, ibnd, lbnd, itol
  ! counters on iteration, bands
  integer , allocatable :: conv (:)
  ! if 1 the root is converged
  !HL NB GWTC: g t h hold -> SGW: r q p pold
  complex(DP), allocatable :: g (:,:), t (:,:), h (:,:), hold (:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space

  COMPLEX(DP)  :: cw

! HL need to introduce gt tt ht htold for BICON
! also gp grp for preconditioned systems

  complex(DP), allocatable :: gt (:,:), tt (:,:), ht (:,:), htold (:,:)
  complex(DP), allocatable :: gp (:,:), gtp (:,:)
  complex(DP) ::  dcgamma, dclambda, alpha, beta
  !  the ratio between rho
  !  step length
  complex(DP), external :: zdotc
  !HL (eigenvalue + iw) 

  complex(DP) :: e(nbnd), eu(nbnd)

  ! the scalar product
  real(DP), allocatable :: rho (:), a(:), c(:), astar(:), cstar(:), b(:)
  ! the residue
  ! auxiliary for h_diag
  real(DP) :: kter_eff
  ! account the number of iterations with b
  ! coefficient of quadratic form
  allocate ( g(ndmx*npol,nbnd), t(ndmx*npol,nbnd), h(ndmx*npol,nbnd), &
             hold(ndmx*npol ,nbnd) )
  allocate ( gt(ndmx*npol,nbnd), tt(ndmx*npol,nbnd), ht(ndmx*npol,nbnd), &
             htold(ndmx*npol, nbnd) )
  allocate ( gp(ndmx*npol,nbnd), gtp(ndmx*npol,nbnd))
  allocate (a(nbnd), c(nbnd))
  allocate (conv ( nbnd))
  allocate (rho(nbnd))
  if(itol==1)allocate(b(nbnd))
 !WRITE(6,*) g,t,h,hold
 !Initialize
  eu(:) = dcmplx(0.0d0,0.0d0)
  kter_eff = 0.d0

  do ibnd = 1, nbnd
     conv (ibnd) = 0
  enddo
 
  conv_root = .false.

  g=(0.d0,0.d0)
  t=(0.d0,0.d0)
  h=(0.d0,0.d0)
  hold=(0.d0,0.d0)

  gt = (0.d0,0.d0)
  tt = (0.d0,0.d0)
  ht = (0.d0,0.d0)
  htold = (0.d0,0.d0)
  gp(:,:) = (0.d0, 0.0d0)
  gtp(:,:) = (0.d0, 0.0d0)

  call start_clock ('cbcgsolve')

  do iter = 1, maxiter
    ! kter = kter + 1
    ! g    = (-PcDv\Psi) - (H \Delta\Psi)   ! AX
    ! gt   = conjg( g)
    ! r    = b - Ax 
    ! rt   = conjg ( r )
     if (iter .eq. 1) then
        !r = b - A* x
        !rt = conjg (r) 
        call h_psi (ndim, dpsi, g, e, cw, ik, nbnd)  !g=Ax   ! (H-e+\alpha |p><p|)dpsi left hand side of the linear equation
        do ibnd = 1, nbnd
           call zaxpy (ndmx*npol, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1)
           IF(itol==1) b(ibnd) = abs(ZDOTC (ndmx*npol, d0psi(1,ibnd), 1, d0psi(1,ibnd), 1)) !calculate |b|^2      
        end do 

#ifdef __MPI
call mp_sum(  b(1:nbnd) , intra_pool_comm )
#endif       
           ! After this step, g=Ax-b

        do ibnd = 1, nbnd
           call zscal (ndmx*npol, (-1.0d0, 0.0d0), g(1,ibnd), 1) !g=r_0=b-Ax
!           gt(:,ibnd) = dconjg ( g(:,ibnd) )       ! obtain r_0^*
           gt(:,ibnd) =  g(:,ibnd)   
           ! HL: for convergence with spin polarized along y
           ! KC: This should not change the results as along as it converges
        !  p   =  inv(M) * r
        !  pt  =  conjg ( p )
           call zcopy (ndmx*npol, g (1, ibnd), 1, h (1, ibnd), 1)  !Copy g to h
           if(tprec) call cg_psi(ndmx, ndim, 1, h(1,ibnd), h_diag(1,ibnd) )
!           ht(:,ibnd) = dconjg( h(:,ibnd) )  
           ht(:,ibnd) =  h(:,ibnd)    
           ! HL: for convergence with spin polarized along y
        enddo
     endif

!HL: Convergence check... 
     lbnd = 0
     do ibnd = 1, nbnd
        if (conv (ibnd).eq.0) then
            lbnd = lbnd+1
            rho(lbnd) = abs(ZDOTC (ndmx*npol, g(1,ibnd), 1, g(1,ibnd), 1))
            if(itol==1)rho(lbnd)=rho(lbnd)/b(ibnd)
             
         !!!!!KC: Here the convergence check is different from cg, why?
         !!!!! This is not a relative convergence threshold 
        endif
     enddo

     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)

! sum within pools
#ifdef __MPI
call mp_sum(  rho(1:lbnd) , intra_pool_comm )
#endif

     do ibnd = nbnd, 1, -1
        if (conv(ibnd).eq.0) then
            rho(ibnd) = rho(lbnd)
            lbnd = lbnd-1
            anorm = sqrt(rho(ibnd))
            if (anorm.lt.ethr) conv (ibnd) = 1
        endif
     enddo

     conv_root = .true.
     do ibnd = 1, nbnd
        conv_root = conv_root.and.(conv (ibnd).eq.1)
     enddo

     if (conv_root) goto 100
    
 if (iter.eq.maxiter .and. .not.conv_root) then
   do ibnd=1, nbnd
      if(conv(ibnd)/=1)then
      WRITE( stdout, '(5x,"kpoint",i4," ibnd",i4, &
                &              " solve_linter: root not converged ",e10.3)') &
                &              ik , ibnd, anorm
      end if
   end do
 end if 
! we only apply hamiltonian to unconverged bands.
     lbnd = 0 
     do ibnd = 1, nbnd
        if (conv(ibnd).eq.0) then 
            lbnd = lbnd + 1
            call zcopy(ndmx*npol, h(1,ibnd),  1, hold(1,  lbnd), 1)
            call zcopy(ndmx*npol, ht(1,ibnd), 1, htold(1, lbnd), 1)
            eu(lbnd) = e(ibnd)
        endif
     enddo

!****************** THIS IS THE MOST EXPENSIVE PART**********************!
! compute t = A*h and tt=(A^H)*ht

     call h_psi (ndim, hold, t, eu(1), cw, ik, lbnd)  ! t=A*h
     call h_psi (ndim, htold, tt, eu(1), dconjg(cw), ik, lbnd)  ! tt=(A^H)*ht

    lbnd=0
    do ibnd = 1, nbnd
       if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1
!alpha = <rt|rp>/<pt|q>
           call ZCOPY (ndmx*npol, g  (1, ibnd), 1, gp  (1, ibnd), 1)  !copy g to gp
           if (tprec) call cg_psi (ndmx, ndim, 1, gp(1,ibnd), h_diag(1,ibnd) )
           a(lbnd) = ZDOTC (ndmx*npol, gt(1,ibnd), 1, gp(1,ibnd), 1)  ! a=<r^|r>
           c(lbnd) = ZDOTC (ndmx*npol, ht(1,ibnd), 1, t (1,lbnd), 1)  ! c=<h|Ah>
       endif
    enddo

!when you have ZDOTC which corresponds to a sum within each band, you need intra_pool_comm operation

#ifdef __MPI          
     call mp_sum(  a(1:lbnd), intra_pool_comm )
     call mp_sum(  c(1:lbnd), intra_pool_comm )
#endif

     lbnd=0
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1 
           alpha = a(lbnd) / c(lbnd)
! x  = x  + alpha        * p
           call ZAXPY (ndmx*npol,  alpha,  h(1,ibnd), 1, dpsi(1,ibnd), 1)
! r  = r  - alpha        * q
! rt = rt - conjg(alpha) * qt
           call ZAXPY (ndmx*npol, -alpha,        t  (1, lbnd), 1, g  (1,ibnd), 1)
           call ZAXPY (ndmx*npol, -dconjg(alpha), tt (1, lbnd), 1, gt (1,ibnd), 1)
! rp  = inv(M) * r
! rtp = inv(M) * rt
           call ZCOPY (ndmx*npol, g  (1, ibnd), 1, gp  (1, ibnd), 1)
           call ZCOPY (ndmx*npol, gt (1, ibnd), 1, gtp (1, ibnd), 1)
           if (tprec) call cg_psi (ndmx, ndim, 1, gp  (1,ibnd), h_diag(1,ibnd) )
           if (tprec) call cg_psi (ndmx, ndim, 1, gtp (1,ibnd), h_diag(1,ibnd) )
! beta = - <qt|rp>/<pt|q>
           a(lbnd) = ZDOTC (ndmx*npol, tt(1,lbnd), 1, gp(1,ibnd), 1)

        end if
     end do

#ifdef __MPI          
     call mp_sum(  a(1:lbnd), intra_pool_comm )
#endif 
  
     lbnd = 0
     do ibnd =1, nbnd
        if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1         
           beta = - a(lbnd) / c(lbnd)
! pold  = p
! ptold = pt
! Extra Copy
           call ZCOPY (ndmx*npol, h  (1, ibnd), 1, hold  (1, ibnd), 1)
           call ZCOPY (ndmx*npol, ht (1, ibnd), 1, htold (1, ibnd), 1)
! p  = rp  +       beta  * pold
! pt = rtp + conjg(beta) * ptold
           call ZCOPY (ndmx*npol, gp  (1, ibnd), 1, h  (1, ibnd), 1)
           call ZCOPY (ndmx*npol, gtp (1, ibnd), 1, ht (1, ibnd), 1)
           call ZAXPY (ndmx*npol,       beta,  hold  (1,ibnd), 1, h (1,ibnd), 1)
           call ZAXPY (ndmx*npol, dconjg(beta), htold (1,ibnd), 1, ht(1,ibnd), 1)
        endif
     enddo
  enddo

100 continue
  kter = kter_eff
  deallocate (rho)
  deallocate (conv)
  deallocate (a,c)
  deallocate (g, t, h, hold)
  deallocate (gt, tt, ht, htold)
  deallocate (gtp, gp)
  call stop_clock ('cbcgsolve')
  return
END SUBROUTINE cbcg_solve
 
