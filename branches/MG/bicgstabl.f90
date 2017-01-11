SUBROUTINE cbicgstabl(h_psi, cg_psi, e, d0psi, dpsi, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol, cw, lmresloc, tprec)
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
             niters, &  ! number of iterations for this BiCG min
             nrec, &    ! for composite rec numbers
             lmresloc     ! set number of loops over minimal residuals...

real(DP) :: &
             anorm,   &        ! output: the norm of the error in the solution
             ethr,    &        ! input: the required precision
             h_diag(ndmx*npol,nbnd) ! input: an estimate of ( H - \epsilon )


COMPLEX(DP) :: alphabeta(2), beta_old

complex(DP) :: &
             dpsi (ndmx*npol, nbnd), & ! output: the solution of the linear syst
             d0psi (ndmx*npol, nbnd)   ! input: the known term

!GMRES part
  complex(DP) :: tau(0:lmresloc ,0:lmresloc)
  complex(DP) :: sigma(0:lmresloc)
  complex(DP) :: gmgamma(0:lmresloc)
  complex(DP) :: gmgammapp(0:lmresloc )
  complex(DP) :: gmgammap(0:lmresloc)
  complex(DP) :: gschmidt

logical :: conv_root ! output: if true the root is converged
external h_psi       ! input: the routine computing h_psi
external cg_psi      ! input: the routine computing cg_psi

!
!  here the local variables
!

  !HL upping iterations to get convergence with green_linsys?
  integer, parameter :: maxiter = 1
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
  complex(DP) ::  dcgamma, dclambda, beta
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
  ! account the number of iterations with b
  ! coefficient of quadratic form
  !
  call start_clock ('cgsolve')

 WRITE(stdout, *)'before allocation'

  allocate ( g(ndmx*npol,nbnd, 0:lmresloc), t(ndmx*npol,nbnd), & 
         h(ndmx*npol,nbnd, 0:lmresloc), hold(ndmx*npol, nbnd))
           
  allocate (rho(nbnd),rhoold(nbnd), omega(nbnd), alpha(nbnd))
  allocate (gt(ndmx*npol,nbnd))
             
!  allocate ( gp(ndmx*npol,nbnd), gtp(ndmx*npol,nbnd))
!  allocate (a(nbnd), c(nbnd))
  allocate (conv ( nbnd))

WRITE(stdout, *)'after allocation'

  kter_eff = 0.d0

  do ibnd = 1, nbnd
     conv (ibnd) = 0
  enddo
  conv_root = .false.

  g        = dcmplx(0.d0,0.d0)
  t        = dcmplx(0.d0,0.d0)
  h        = dcmplx(0.d0,0.d0)
  hold     = dcmplx(0.d0,0.d0)
  gt       = dcmplx(0.d0,0.d0)
!  gp(:,:)  = dcmplx(0.d0,0.0d0)
!  gtp(:,:) = dcmplx(0.d0,0.0d0)
  rho      = dcmplx(1.0d0,0.0d0)
  rhoold   = dcmplx(1.0d0,0.0d0)
  omega    = dcmplx(1.0d0,0.0d0)
  alpha    = dcmplx(0.0d0,0.0d0)
!  sigma    = dcmplx(1.0d0, 0.0d0)

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
            if (iter .eq. 1) then
               !initialize iter 0 stuff
               !r = b - A* x
               !rt = conjg (r) 
!               call h_psi (ndim, dpsi, g(:, :, 0), e, cw, ik, nbnd)

               do ibnd = 1, nbnd
               !initial residual should be r = b
               ! I/O for multishift.
               ! call davcio (d0psi(:,1), lrresid, iunresid, iter, +1)
                  call zaxpy (ndmx*npol, (1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd, 0), 1)
!                  call zscal (ndmx*npol, (-1.0d0, 0.0d0), g(1,ibnd,0), 1)


!           if(tprec) call cg_psi(ndmx, ndim, 1, g(1,ibnd), h_diag(1,ibnd) )                                                                        

    ! copy r -> u, i.e. u = r_{0}
    ! call zcopy (ndim, g (1, ibnd), 1, h (1, ibnd), 1)
    ! set \tilde{r} = r^{*}
    ! gt(:,ibnd)    = conjg ( g(:,ibnd) )
                  gt(:,ibnd) = g(:,ibnd,0)
    !for preconditioning we solve the system MAx = Mb
    !x shouldn't change...
    !            gt = Mb
               enddo
            endif

WRITE(stdout, *)'after initilization'
!####################### begin convergence check #############################
     lbnd = 0
     do ibnd = 1, nbnd
        if (conv (ibnd).eq.0) then
            lbnd = lbnd+1
            rho(lbnd) = abs(ZDOTC (ndmx*npol, g(1,ibnd, 0), 1, g(1,ibnd, 0), 1))
!            if(itol==1)rho(lbnd)=rho(lbnd)/b(ibnd)

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

 !if (iter.eq.maxiter .and. .not.conv_root) then
 !  do ibnd=1, nbnd
 !     if(conv(ibnd)/=1)then
 !     WRITE( stdout, '(5x,"kpoint",i4," ibnd",i4, &
 !               &              " solve_linter: root not converged ",e10.3)') &
 !               &              ik , ibnd, sqrt(rho(ibnd))
 !     end if
 !  end do
 !end if

WRITE(stdout, *)'after convergence check'
!!##################### END convergence check ###########################
!SUBROUTINE SAXPY(N,A,X,INCX,Y,INCY)
!Y = A * X + Y

Do iterj=0, lmresloc-1

!     rhoold(:)  = -omega(:)*rhoold(:)

 Do ibnd=1,nbnd
   IF(conv(ibnd)/=1)THEN

     rhoold(ibnd)  = -omega(ibnd)*rhoold(ibnd)
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

!           call h_psi (ndim, h(:,ibnd), hn(:,ibnd), e(1), cw, ik, nbnd)
!           if(tprec) call cg_psi(ndmx, ndim, 1, hn(1,ibnd), h_diag(1,ibnd))
    lbnd = 0                                                                  
    do ibnd = 1, nbnd                                                         
       if (conv(ibnd).eq.0) then
           lbnd = lbnd + 1
           call zcopy(ndmx*npol, h(:,ibnd,iterj),  1, hold(:,  lbnd), 1)
!           call zcopy(ndmx*npol, ht(1,ibnd), 1, htold(1, lbnd), 1)            
           eu(lbnd) = e(ibnd)                                                 
       endif
    enddo

    call h_psi (ndim, hold, t, eu(1), cw, ik, lbnd)

!    lbnd = 0
!    DO ibnd = 1,nbnd
!      IF (conv(ibnd).eq.0) then
!         lbnd = lbnd + 1
!         call zcopy(ndmx*npol, t(1,lbnd),  1, h(1,  ibnd, iterj+1), 1)
!      END IF 
!    END DO

   lbnd=0
 DO ibnd = 1, nbnd
   IF (conv (ibnd) .eq.0) THEN
          lbnd=lbnd+1

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

      call zcopy(ndmx*npol, t(:,lbnd),  1, g(:,  ibnd, iterj+1), 1)

      END IF ! conv
   END DO   ! ibnd


END DO   ! iterj


WRITE(stdout, *)'after BICG part'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!Min Residual Part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO ibnd = 1, nbnd

!################remember to try a version with all the coeff ibnd dependent
!################## to test consistency ###############!!!!!!!!!!!!!!!!

    IF(conv(ibnd)/=1)THEN



     do iterj = 1, lmresloc

       IF(iterj>1)THEN
          do iteri = 1, iterj-1
            tau(iteri, iterj) = (1.0d0/sigma(iteri))*ZDOTC(ndmx*npol, g(1,ibnd, iterj), g(1,ibnd, iteri), 1)
!orthogonalizes each residual:
            call ZAXPY (ndmx*npol, -tau(iteri, iterj),  g(1,ibnd,iteri), 1, g(1,ibnd, iterj), 1)
           enddo
       END IF
    
!should correspond to van der vorst's:sleipjen:SGW
!w=(t,s)(t,t):(\sigma_{1} = (\r_{1},\r_{1})), \gammap=(1/\sigma_{j})(\r_{0}, \r_{j}): 
!sigma=(t,t) gmgammap=(1/sigma(j))*(ro,t)
!sigma(iterj)   = ZDOTC(ndmx*npol, g(1,iterj), 1, g(1,iterj), 1)
        sigma(iterj)   = ZDOTC(ndmx*npol, g(1,ibnd,iterj), 1, g(1,ibnd, iterj), 1) 
!should be r_{1}! i.e. t 
!gmgammap(iterj) = (1.0d0/sigma(iterj)) * ZDOTC(ndim, g(1,ibnd), 1, g(1,iterj), 1)
!       gmgammap(iterj) = (1.0d0/sigma(iterj)) * ZDOTC(ndim, g(1,ibnd), 1, t(1,ibnd), 1)
        gmgammap(iterj) = (1.0d0/sigma(iterj)) * ZDOTC(ndmx*npol, g(1,ibnd,0), 1, g(1,ibnd,iterj), 1)

     enddo  !iterj


!w = (t,s)(t,t) = (\sigma_{j} = (\r_{j}, \r_{j})), \gammap = (1/\sigma_{j})(\r_{0}, \r_{0})
     gmgamma(lmresloc) = gmgammap(lmresloc)
     omega(ibnd)      = gmgamma(lmresloc)


!write(600,'(6f12.7)') gmgamma(:), omega
!ok what's goin on here??
     IF(lmresloc>1)THEN

         do iterj = lmresloc-1, 1, -1
            gschmidt = dcmplx(0.0d0, 0.0d0)
         do iteri = iterj + 1, lmresloc
            gschmidt = gschmidt + tau(iterj, iteri)*gmgamma(iteri)
         enddo
            gmgamma(iterj) = gmgammap(iterj) - gschmidt
         enddo


        do iterj = 1, lmresloc-1
           gschmidt = dcmplx(0.0d0, 0.0d0)
           do iteri = iterj + 1, lmresloc-1
              gschmidt = gschmidt + tau(iterj, iteri)*gmgamma(iteri + 1)
           enddo
              gmgammapp(iterj) = gmgamma(iterj+1) + gschmidt
        enddo

     END IF
!The second update phase:
!   do ibnd = 1, nbnd
!      x  = x  + gmgamma(1)  * \hat{r}_{0}
       call ZAXPY (ndmx*npol,  gmgamma(1), g(1,ibnd, 0),  1, dpsi(1,ibnd), 1)
!     \hat{u}_{0}  = \hat{u}_{0}  +  beta  * u_old
       call ZAXPY (ndmx*npol, -gmgamma(lmresloc), h(1,ibnd, lmresloc), 1, h(1,ibnd, 0), 1)
!     \hat{r}_{0}  = \hat{r}_{0}  - gammap_{l}*r_{lmresloc}
       call ZAXPY (ndmx*npol, -gmgammap(lmresloc), g(1,ibnd,lmresloc), 1, g(1,ibnd,0), 1)
      !and again here we get to a point and then stop progressing...
      !call ZAXPY (ndmx*npol, (gmgamma(lmresloc+1)), hn(1,ibnd), 1, h(1,ibnd), 1)
      IF(lmresloc.gt.1) then

        do iterj = 1, lmresloc-1
!      \hat{u_0}  = \hat{u_0} + gmgamma(lmresloc)*u_{j}
        call ZAXPY (ndmx*npol, -(gmgamma(iterj)),   h  (1,ibnd,iterj), 1, h (1,ibnd, 0), 1)
!      x_{0}  = x_{0} + gmgamma''_{j}  * r_{j}
        call ZAXPY (ndmx*npol,  gmgammapp(iterj), g(1,ibnd,iterj),  1, dpsi  (1,ibnd), 1)
!      r_{0}  = r_{0}  - gamma'_{j}*r_{j}
        call ZAXPY (ndmx*npol, -(gmgammap(iterj)), g(1,ibnd, iterj), 1, g  (1,ibnd,0), 1)
        enddo

      END IF !lmresloc > 1
 
     endif !conv(ibnd)

   enddo !ibnd

enddo!maxter

100 continue

  kter   =  kter_eff
  deallocate (conv)
  deallocate (rho, rhoold, omega, alpha)
  deallocate (g, t, h, hold)
  deallocate (gt)
!  deallocate (gtp, gp)
  call stop_clock ('cgsolve')
  return
END SUBROUTINE cbicgstabl
 
