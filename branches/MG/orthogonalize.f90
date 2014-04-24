!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq)
!------------------------------------------------------------------------
  !
  ! This routine ortogonalizes dvpsi to the valence states: ps = <evq|dvpsi>
  ! It should be quite general. It works for metals and insulators, with
  ! NC as well as with US PP, both SR or FR.
  ! Note that on output it changes sign. So it applies -P^+_c.
  !
  ! NB: IN/OUT is dvpsi ; dpsi is used as work_space
  !
USE kinds, ONLY : DP
USE klist, ONLY : lgauss, degauss, ngauss
USE noncollin_module, ONLY : noncolin, npol
USE wvfct, ONLY : npwx, nbnd, et
USE ener, ONLY : ef
USE control_ph,  ONLY : alpha_pv, nbnd_occ
USE becmod,      ONLY : bec_type, becp, calbec
USE uspp,        ONLY : vkb, okvan
USE mp_global,   ONLY : intra_pool_comm
USE mp,          ONLY : mp_sum
USE control_flags, ONLY : gamma_only
USE realus,      ONLY : npw_k
USE gvect,       ONLY : gstart
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ikk, ikq   ! the index of the k and k+q points
INTEGER, INTENT(IN) :: npwq       ! the number of plane waves for q
COMPLEX(DP), INTENT(IN) :: evq(npwx*npol,nbnd)
COMPLEX(DP), INTENT(INOUT) :: dvpsi(npwx*npol,nbnd)
COMPLEX(DP), INTENT(INOUT) :: dpsi(npwx*npol,nbnd) ! work space allocated by
                                                   ! the calling routine

COMPLEX(DP), ALLOCATABLE :: ps(:,:)
REAL(DP), ALLOCATABLE :: ps_r(:,:)

INTEGER :: ibnd, jbnd, nbnd_eff
REAL(DP) :: wg1, w0g, wgp, wwg, deltae, theta
REAL(DP), EXTERNAL :: w0gauss, wgauss
! functions computing the delta and theta function

CALL start_clock ('ortho')
IF (gamma_only) THEN
   ALLOCATE(ps_r(nbnd,nbnd))
   ps_r = 0.0_DP
ENDIF

ALLOCATE(ps(nbnd,nbnd))
ps = (0.0_DP, 0.0_DP)

!
if (lgauss) then
   !
   IF (gamma_only) CALL errore ('orthogonalize', "degauss with gamma &
        & point algorithms",1)
   !
   !  metallic case
   !
   IF (noncolin) THEN
      CALL zgemm( 'C', 'N', nbnd, nbnd_occ (ikk), npwx*npol, (1.d0,0.d0), &
                 evq, npwx*npol, dvpsi, npwx*npol, (0.d0,0.d0), ps, nbnd )
   ELSE
      CALL zgemm( 'C', 'N', nbnd, nbnd_occ (ikk), npwq, (1.d0,0.d0), &
                 evq, npwx, dvpsi, npwx, (0.d0,0.d0), ps, nbnd )
   END IF
   !
   DO ibnd = 1, nbnd_occ (ikk)
      wg1 = wgauss ((ef-et(ibnd,ikk)) / degauss, ngauss)
      w0g = w0gauss((ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
      DO jbnd = 1, nbnd
         wgp = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
         deltae = et (jbnd, ikq) - et (ibnd, ikk)
         theta = wgauss (deltae / degauss, 0)
         wwg = wg1 * (1.d0 - theta) + wgp * theta
         IF (jbnd <= nbnd_occ (ikq) ) THEN
            IF (abs (deltae) > 1.0d-5) THEN
                wwg = wwg + alpha_pv * theta * (wgp - wg1) / deltae
            ELSE
               !
               !  if the two energies are too close takes the limit
               !  of the 0/0 ratio
               !
               wwg = wwg - alpha_pv * theta * w0g
            ENDIF
         ENDIF
         !
         ps(jbnd,ibnd) = wwg * ps(jbnd,ibnd)
         !
      ENDDO
      IF (noncolin) THEN
         CALL dscal (2*npwx*npol, wg1, dvpsi(1,ibnd), 1)
      ELSE
         call dscal (2*npwq, wg1, dvpsi(1,ibnd), 1)
      END IF
   END DO
   nbnd_eff=nbnd
ELSE
   !
   !  insulators
   !
   IF (noncolin) THEN
      CALL zgemm( 'C', 'N',nbnd_occ(ikq), nbnd_occ(ikk), npwx*npol, &
             (1.d0,0.d0), evq, npwx*npol, dvpsi, npwx*npol, &
             (0.d0,0.d0), ps, nbnd )
   ELSEIF (gamma_only) THEN
            CALL dgemm( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), 2*npwq, &
             2.0_DP, evq, 2*npwx, dvpsi, 2*npwx, &
             0.0_DP, ps_r, nbnd )
            IF (gstart == 2 ) THEN
               CALL DGER( nbnd_occ(ikq), nbnd_occ (ikk), -1.0_DP, evq, &
                    & 2*npwq, dvpsi, 2*npwx, ps_r, nbnd )
            ENDIF
   ELSE
      CALL zgemm( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), npwq, &
             (1.d0,0.d0), evq, npwx, dvpsi, npwx, &
             (0.d0,0.d0), ps, nbnd )
   END IF
   nbnd_eff=nbnd_occ(ikk)
END IF
#ifdef __MPI
IF (gamma_only) THEN
   call mp_sum(ps_r(:,:),intra_pool_comm)
ELSE
   call mp_sum(ps(:,1:nbnd_eff),intra_pool_comm)
ENDIF
#endif
!
! dpsi is used as work space to store S|evc>
!
IF (okvan) CALL calbec ( npwq, vkb, evq, becp, nbnd_eff)
CALL s_psi (npwx, npwq, nbnd_eff, evq, dpsi)
!
! |dvspi> =  -(|dvpsi> - S|evq><evq|dvpsi>)
!
if (lgauss) then
   !
   !  metallic case
   !
   IF (noncolin) THEN
      CALL zgemm( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd, &
                (1.d0,0.d0), dpsi, npwx*npol, ps, nbnd, (-1.0d0,0.d0), &
                dvpsi, npwx*npol )
   ELSE
      CALL zgemm( 'N', 'N', npwq, nbnd_occ(ikk), nbnd, &
             (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
              dvpsi, npwx )
   END IF
ELSE
   !
   !  Insulators: note that nbnd_occ(ikk)=nbnd_occ(ikq) in an insulator
   !
   IF (noncolin) THEN
      CALL zgemm( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd_occ(ikk), &
                (1.d0,0.d0),dpsi,npwx*npol,ps,nbnd,(-1.0d0,0.d0), &
                dvpsi, npwx*npol )
   ELSEIF (gamma_only) THEN
      ps = CMPLX (ps_r,0.0_DP, KIND=DP)
      CALL ZGEMM( 'N', 'N', npwq, nbnd_occ(ikk), nbnd_occ(ikk), &
               (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
                dvpsi, npwx )
   ELSE
      CALL zgemm( 'N', 'N', npwq, nbnd_occ(ikk), nbnd_occ(ikk), &
             (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
              dvpsi, npwx )
   END IF
ENDIF

DEALLOCATE(ps)
CALL stop_clock ('ortho')
RETURN
END SUBROUTINE orthogonalize
