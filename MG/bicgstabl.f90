! This file is copied and modified from QUANTUM ESPRESSO
! Kun Cao, Henry Lambert, Feliciano Giustino
 
SUBROUTINE cbicgstabl(h_psi, cg_psi, e, d0psi, dpsi, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol, cw, lmresloc, tprec, itol)
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
USE mp_global,   ONLY: intra_pool_comm
USE mp,          ONLY: mp_sum
!USE control_gw,  ONLY: maxter_green
!USE units_gw,    ONLY: iunresid, lrresid, iunalphabeta, lralphabeta
USE mp_global,   ONLY: inter_pool_comm, intra_pool_comm, mp_global_end, mpime, &
                        nproc_pool, nproc, me_pool, my_pool_id, npool
USE mp,          ONLY : mp_barrier, mp_bcast, mp_sum
USE io_global,   ONLY: stdout

implicit none

!tprec true or false conditioning?
logical :: tprec

integer ::   ndmx, &   ! input: the maximum dimension of the vectors
             ndim, &   ! input: the actual dimension of the vectors
             kter, &   ! output: counter on iterations
             nbnd, &   ! input: the number of bands
             npol, &   ! input: number of components of the wavefunctions
             ik,   &   ! input: the k point
             nrec, &    ! for composite rec numbers
             lmresloc, &    ! set number of loops over minimal residuals...
             itol

real(DP) :: &
             anorm,   &        ! output: the norm of the error in the solution
             ethr,    &        ! input: the required precision
             h_diag(ndmx*npol,nbnd) ! input: an estimate of ( H - \epsilon )


!COMPLEX(DP) :: alphabeta(2), beta_old

complex(DP) :: &
             dpsi (ndmx*npol, nbnd), & ! output: the solution of the linear syst
             d0psi (ndmx*npol, nbnd)   ! input: the known term

!GMRES part
  complex(DP) :: tau(lmresloc ,lmresloc)
  complex(DP) :: gmgamma(lmresloc), gmgammat(lmresloc)

logical :: conv_root ! output: if true the root is converged
external h_psi       ! input: the routine computing h_psi
external cg_psi      ! input: the routine computing cg_psi
external inv

!
!  here the local variables
!

  !HL upping iterations to get convergence with green_linsys?
  integer, parameter :: maxiter = 200
  !integer, parameter :: maxter = 600
  !the maximum number of iterations
  integer :: iter, ibnd, lbnd, iterj, iteri
  ! counters on iteration, bands
  integer , allocatable :: conv (:)
  ! if 1 the root is converged
  !HL NB GWTC: g t h hold -> SGW: r q p pold
  complex(DP), allocatable :: g (:,:,:),  &
                            h (:,:,:), hold (:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space
  COMPLEX(DP)  :: cw
  !HL need to introduce gt tt ht htold for BICON
  ! also gp grp for preconditioned systems

  complex(DP), allocatable :: t(:,:), gt (:,:)
  complex(DP), allocatable :: alpha (:), omega (:)
  complex(DP) ::  dcgamma,  beta
  !  the ratio between rho
  !  step length
  complex(DP), external :: zdotc
  !HL (eigenvalue + iw) 

  complex(DP) :: e(nbnd), eu(nbnd)

  ! the scalar product
  complex(DP), allocatable :: rho(:), rhoold(:)
  ! the residue
  ! auxiliary for h_diag
  real(DP) :: kter_eff
  real(DP), allocatable:: rho0(:), b(:)
  ! account the number of iterations with b
  ! coefficient of quadratic form
  !
  call start_clock ('cgsolve')

! WRITE(stdout, *)'before allocation'

  allocate ( g(ndmx*npol,nbnd, 0:lmresloc), t(ndmx*npol,nbnd), & 
         h(ndmx*npol,nbnd, 0:lmresloc), hold(ndmx*npol, nbnd))
           
  allocate (rho(nbnd),rhoold(nbnd), omega(nbnd), alpha(nbnd))
  allocate (gt(ndmx*npol,nbnd))
             
!  allocate ( gp(ndmx*npol,nbnd), gtp(ndmx*npol,nbnd))
!  allocate (a(nbnd), c(nbnd))
  allocate (conv ( nbnd))
  allocate (rho0(nbnd),b(nbnd))

!WRITE(stdout, *)'after allocation'

  kter_eff = 0.d0

  do ibnd = 1, nbnd
     conv (ibnd) = 0
  enddo
  conv_root = .false.

  g        = dcmplx(0.d0,0.d0)
  t        = dcmplx(0.d0,0.d0)
  h        = dcmplx(0.d0,0.d0)
  hold     = dcmplx(0.d0,0.d0)
  gt       = dcmplx(1.d0,1.d0)
!  gp(:,:)  = dcmplx(0.d0,0.0d0)
!  gtp(:,:) = dcmplx(0.d0,0.0d0)
  rho      = dcmplx(1.0d0,0.0d0)
  rhoold   = dcmplx(1.0d0,0.0d0)
  omega    = dcmplx(1.0d0,0.0d0)
  alpha    = dcmplx(1.0d0,0.0d0)
gmgamma  = dcmplx(0.0d0,0.0d0)
tau = dcmplx(0.0d0,0.0d0)

  do iter = 1, maxiter
  ! r    = b - Ax 
  ! rt   = conjg ( r )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!cBICG PART!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! iterj allows for the bicgstab(2,3,4,...) generalizations.
  ! a choice of lmresloc = 1 corresponds to the standard bicgstab.
  !H. A. Van Der Corst Siam J. Sci. Stat. Comput.
  !Vol. 13, No. 2, pp. 631-644
!            if (mod(iter, 20) .eq. 1) then
             if (iter .eq. 1) then                           
               !initialize iter 0 stuff
               !r = b - A* x
               !rt = conjg (r) 
               call h_psi (ndim, dpsi, g(:, :, 0), e, cw, ik, nbnd)

               do ibnd = 1, nbnd
!                  IF(conv(ibnd) .eq. 0)then

!############################# for restart ################################
!                  h(:,ibnd,0)= (0.d0, 0.d0)
!                  alpha(ibnd)=(1.d0, 0.d0)
!                  rhoold(ibnd)=(1.d0, 0.d0)
!                  omega(ibnd)=(1.d0, 0.d0)
!##########################################################################                  
               !initial residual should be r = b
               ! I/O for multishift.
               ! call davcio (d0psi(:,1), lrresid, iunresid, iter, +1)
                  call zaxpy (ndmx*npol, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd, 0), 1)
                  call zscal (ndmx*npol, (-1.0d0, 0.0d0), g(1,ibnd,0), 1)
                  if(tprec) call cg_psi(ndmx, ndim, 1, g(1,ibnd,0), h_diag(1,ibnd))

!                  IF(itol==1)THEN
                      call zcopy(ndmx*npol, d0psi(:,ibnd),  1, hold(:,  ibnd), 1)
                      IF(tprec) call cg_psi(ndmx, ndim, 1, hold(1,ibnd),h_diag(1,ibnd))
                      b(ibnd)= ABS(ZDOTC (ndmx*npol, hold(1,ibnd), 1, hold(1,ibnd), 1))
                      IF(ik==1)write(stdout,'(i5, e10.3)')ibnd, b(ibnd)
!                  END IF
    ! copy r -> u, i.e. u = r_{0}
    ! call zcopy (ndim, g (1, ibnd), 1, h (1, ibnd), 1)
    ! set \tilde{r} = r^{*}
    ! gt(:,ibnd)    = conjg ( g(:,ibnd) )
                     gt(:,ibnd) = g(:,ibnd,0)
    !for preconditioning we solve the system MAx = Mb
    !x shouldn't change...
    !            gt = Mb
!                 END IF
               enddo
            endif

  
!WRITE(stdout, *)'after initilization'
!####################### begin convergence check #############################
     lbnd = 0
     do ibnd = 1, nbnd
        if (conv (ibnd).eq.0) then
            lbnd = lbnd+1
!            rho(lbnd) = abs(ZDOTC (ndmx*npol, g(1,ibnd, 0), 1, g(1,ibnd, 0), 1))
!            if(itol==1)rho(lbnd)=rho(lbnd)/b(ibnd)
            rho0(ibnd) = abs(ZDOTC (ndmx*npol, g(1,ibnd, 0), 1, g(1,ibnd, 0), 1))
            IF(itol==1)rho0(ibnd)=rho0(ibnd)/b(ibnd)
            anorm = sqrt(rho0(ibnd))

            IF (anorm.lt.ethr) THEN
               conv (ibnd) = 1
            ENDIF
         !!h!!KC: Here the convergence check is different from cg, why?
         !!!!! This is not a relative convergence threshold 
        endif
     enddo

     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)

! sum within pools
!#ifdef __MPI
!call mp_sum(  rho(1:lbnd) , intra_pool_comm )
!#endif

!     do ibnd = nbnd, 1, -1
!        if (conv(ibnd).eq.0) then
!            rho0(ibnd) = rho(lbnd)
!            lbnd = lbnd-1
!            anorm = sqrt(rho0(ibnd))

!            IF (anorm.lt.ethr) THEN
!               conv (ibnd) = 1
!            ENDIF

!        endif
!     enddo

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
                &              ik , ibnd, sqrt(rho0(ibnd))
      end if
   end do
 end if

!WRITE(stdout, *)'after convergence check'
!!##################### END convergence check ###########################
!SUBROUTINE SAXPY(N,A,X,INCX,Y,INCY)
!Y = A * X + Y

rhoold(:)  = -omega(:)*rhoold(:)

Do iterj=0, lmresloc-1

!     rhoold(:)  = -omega(:)*rhoold(:)

 Do ibnd=1,nbnd
   IF(conv(ibnd)/=1)THEN

!     rhoold(ibnd)  = -omega(ibnd)*rhoold(ibnd)
!        do iterj = 0, lmresloc-1
                rho(ibnd) = ZDOTC (ndmx*npol, g(1,ibnd, iterj), 1, gt(1,ibnd), 1)
                beta    = alpha(ibnd)*(rho(ibnd)/rhoold(ibnd))
                rhoold(ibnd)  = rho(ibnd)

                do iteri = 0, iterj
          !\hat{u}_{i} =  \hat{r}_{i} - \beta \hat{u}_{i}
!                   hold(:,:) = dcmplx(0.0d0, 0.d0)
!                   call ZCOPY (ndmx*npol,  g(1, iteri),      1, hold (1, iteri), 1)
                   call ZSCAL (ndmx*npol, -beta, h(1,ibnd, iteri), 1)
                   call ZAXPY (ndmx*npol,  (1.d0, 0.d0), g(1,ibnd, iteri), 1, h(1,ibnd, iteri), 1)
!                   call ZCOPY (ndmx*npol,  hold(1, iteri),   1, h (1, iteri), 1)
                enddo
!       enddo

   END IF
 END DO

     !****************** THIS IS THE MOST EXPENSIVE PART**********************!
     ! HLstab
     ! t = u_{j+1} = A u_{j}

    lbnd = 0
    do ibnd = 1, nbnd
      if (conv(ibnd).eq.0) then
           lbnd = lbnd + 1
           call zcopy(ndmx*npol, h(:,ibnd,iterj),  1, hold(:,  lbnd), 1)
           eu(lbnd) = e(ibnd)
       endif
    enddo

    call h_psi (ndim, hold, t, eu(1), cw, ik, lbnd)


   lbnd=0
 DO ibnd = 1, nbnd
   IF (conv (ibnd) .eq.0) THEN
          lbnd=lbnd+1

          if(tprec) call cg_psi(ndmx, ndim, 1, t(1,ibnd), h_diag(1,ibnd) )

          call zcopy(ndmx*npol, t(1,lbnd),  1, h(1,  ibnd, iterj+1), 1)
        
          dcgamma = ZDOTC (ndmx*npol, h(:,ibnd,iterj+1), 1, gt(:,ibnd), 1)

          alpha(ibnd) = rhoold(ibnd) / dcgamma
       

       do iteri = 0, iterj
       call ZAXPY (ndmx*npol, -alpha(ibnd), h(1,ibnd,iteri+1), 1, g(1,ibnd,iteri), 1)
       enddo

 !update the solution x_{0}  = x_{0}  + alpha * u_{0}

   call ZAXPY (ndmx*npol, alpha(ibnd), h(1,ibnd, 0), 1, dpsi (1,ibnd), 1)

   END IF
 END DO


lbnd = 0       
    DO ibnd = 1, nbnd                                               
       IF (conv(ibnd).eq.0) then
           lbnd = lbnd + 1
           call zcopy(ndmx*npol, g(1,ibnd,iterj),  1, hold(1,  lbnd), 1)

           eu(lbnd) = e(ibnd)
       endif
   enddo

   call h_psi (ndim, hold, t , eu(1), cw, ik, lbnd)

   lbnd =0 
   DO ibnd = 1, nbnd
      IF (conv(ibnd).eq.0) then
      lbnd = lbnd + 1

      if(tprec) call cg_psi(ndmx, ndim, 1, t(1,ibnd), h_diag(1,ibnd) )

      call zcopy(ndmx*npol, t(:,lbnd),  1, g(:,  ibnd, iterj+1), 1)

      END IF ! conv
   END DO   ! ibnd


END DO   ! iterj


!WRITE(stdout, *)'after BICG part'
!call flush(6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!Min Residual Part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!IF(1==0)THEN

  DO ibnd = 1, nbnd

!################remember to try a version with all the coeff ibnd dependent
!################## to test consistency ###############!!!!!!!!!!!!!!!!

    IF(conv(ibnd)/=1)THEN

!  IF(1==0)THEN

      do iteri = 1, lmresloc
          
          do iterj = 1, iteri
            
            tau(iteri, iterj) = ZDOTC(ndmx*npol, g(1,ibnd, iterj), 1, g(1,ibnd, iteri), 1)
!orthogonalizes each residual:
            IF(iteri/=iterj)tau(iterj, iteri) = dconjg(tau(iteri, iterj))
!            IF(iteri/=iterj)tau(iterj, iteri) = tau(iteri, iterj)
          end do

       gmgamma(iteri) = ZDOTC(ndmx*npol, g(1,ibnd,0), 1, g(1,ibnd,iteri), 1)

     enddo  !iterj

! calculate tau^{-1}
     call inv(tau, lmresloc)

     gmgammat(:)= gmgamma(:)

     call ZGEMV('N',lmresloc,lmresloc,(1.d0, 0.d0),tau,lmresloc, gmgammat, 1, (0.d0, 0.d0) , gmgamma, 1)
     
     omega(ibnd)=gmgamma(lmresloc)
     


        do iterj = 1, lmresloc
!      \hat{u_0}  = \hat{u_0} + gmgamma(lmresloc)*u_{j}
        call ZAXPY (ndmx*npol, -(gmgamma(iterj)),   h  (1,ibnd,iterj), 1, h (1,ibnd, 0), 1)
!      x_{0}  = x_{0} + gmgamma''_{j}  * r_{j}
        call ZAXPY (ndmx*npol,  gmgamma(iterj), g(1,ibnd,iterj-1),  1, dpsi  (1,ibnd), 1)
!      r_{0}  = r_{0}  - gamma'_{j}*r_{j}
        call ZAXPY (ndmx*npol, -(gmgamma(iterj)), g(1,ibnd, iterj), 1, g  (1,ibnd,0), 1)
        enddo


     endif !conv(ibnd)

   enddo !ibnd


enddo!maxter

100  continue



  kter   =  kter_eff
  deallocate (conv)
  deallocate (rho, rhoold, omega, alpha, rho0, b)
  deallocate (g, t, h, hold)
  deallocate (gt)
!  deallocate (gtp, gp)
  call stop_clock ('cgsolve')
  return
END SUBROUTINE cbicgstabl
 

! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
Subroutine inv(A, n)
  implicit none
  integer::n, info
  complex(kind=8), intent(inout) :: A(n,n)
!  complex(kind=8), dimension(size(A,1),size(A,2)) :: Ainv

  complex(kind=8), allocatable :: work(:)  ! work array for LAPACK
  integer, allocatable :: ipiv(:)   ! pivot indices

  ! External procedures defined in LAPACK
  external ZGETRF
  external ZGETRI

  allocate(work(n),ipiv(n))
  ! Store A in Ainv to prevent it from being overwritten by LAPACK

 ! write(*,*)'A', A
 ! Write(*,*)'n', n
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call ZGETRF(n, n, A, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call ZGETRI(n, A, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end subroutine inv
