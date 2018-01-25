! This file is copied and modified from QUANTUM ESPRESSO
! Kun Cao, Henry Lambert, Feliciano Giustino
 
!
! Copyright (C) 2001-2008 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE magscf
  !-----------------------------------------------------------------------
  !
  !     This subroutine is the main driver of the self consistent cycle
  !     which gives as output the change of the wavefunctions and the
  !     change of the self-consistent potential due to a phonon of
  !     a fixed q or to an electric field.
  !
  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat
  USE lsda_mod, ONLY : nspin
  USE io_global,  ONLY : stdout, ionode
!  USE fft_base,   ONLY : dfftp
  USE cell_base, ONLY : omega
  USE uspp,  ONLY: okvan
  USE efield_mod, ONLY : zstarue0, zstarue0_rec
  USE control_ph, ONLY : zue, convt, rec_code, do_elec,lgamma, dbext, lrpa, dvext!, do_trans,
  USE partial,    ONLY : done_irr, comp_irr
  USE modes,      ONLY : nirr, npert, npertx
  USE phus,       ONLY : int3, int3_nc, int3_paw
  USE uspp_param, ONLY : nhm
  USE eqv,        ONLY : drhoscfs
  USE paw_variables, ONLY : okpaw
  USE noncollin_module, ONLY : noncolin, nspin_mag, npol
  USE recover_mod, ONLY : write_rec
  USE qpoint,          ONLY : xq
  USE fft_base,   ONLY: dfftp, dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE spin_orb, ONLY : domag

  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum
  USE freq_ph,       ONLY : fpol, fiu, nfs, nfsmax
  USE gvecs,         ONLY : doublegrid, nls  !KC

  IMPLICIT NONE

  INTEGER :: irr, irr1, imode0, npe, ig,iw, ir, i
  ! counter on the representations
  ! counter on the representations
  ! counter on the modes
  ! npert(irr)
  COMPLEX(DP) :: magtot_nc(1:3)
  REAL(DP) :: tcpu, get_clock
  ! timing variables

  LOGICAL :: exst
  ! used to test the recover file

  EXTERNAL get_clock
  ! the change of density due to perturbations

  CALL start_clock ('magscf')

!        ALLOCATE (drhoscfs( dfftp%nnr , nspin_mag))
         ALLOCATE (drhoscfs( dffts%nnr , nspin_mag)) !KC
        imode0 = 0

        WRITE( stdout, '(/,5x,"qpoint= ", 3f12.5)'), xq(1:3)
        WRITE( stdout, '(/,5x,"Self-consistent Calculation")')
        WRITE( stdout, *)'npol,nspin_mag, domag, lrpa', npol, nspin_mag, domag, lrpa
        WRITE( stdout, *)'dfftp,dffts, doublegrid,nls(1)', dfftp%nnr, dffts%nnr, doublegrid,nls(1)

        IF (okvan) THEN
           ALLOCATE (int3 ( nhm, nhm, 1, nat, nspin_mag))
           IF (okpaw) ALLOCATE (int3_paw (nhm, nhm, 1, nat, nspin_mag))
           IF (noncolin) ALLOCATE(int3_nc( nhm, nhm, 1, nat, nspin))
        ENDIF

!        do iw =1, nfs
do iw =1, nfs
!                CALL solve_linter (drhoscfs(1,1), iw)
                WRITE( stdout,'(/,5x,"frequency = ", 2f12.5)'), (fiu(iw))
                WRITE( stdout, '(/,5x,"qpoint= ", 3f12.5)'), xq(1:3)
                CALL solve_linter (drhoscfs(1,1), iw)
                WRITE( stdout, '(/,5x,"End of self-consistent calculation")')
                WRITE( stdout, '(/,5x,"qpoint= ", 3f12.5)'), xq(1:3)
                WRITE( stdout, '(/,5x,"X_[G](Gp)")')
                WRITE( stdout, '("charge density response ")')
                WRITE( stdout, *)

!  KC: output the magnetization in real space for test purpose only 
       if(noncolin)then
          magtot_nc = (0.D0, 0.d0)
!          absmag    = 0.D0
          !
          DO ir = 1,dffts%nnr
             !
             !mag = SQRT( drhoscfs(ir,2)**2 + &
             !            drhoscfs(ir,3)**2 + &
             !            drhoscfs(ir,4)**2 )
             !
             DO i = 1, 3
                !
                magtot_nc(i) = magtot_nc(i) + drhoscfs(ir,i+1)
                !
             END DO
             !
             !absmag = absmag + ABS( mag )
             !
          END DO
          !
          !CALL mp_sum( magtot_nc, intra_bgrp_comm )
          !CALL mp_sum( absmag, intra_bgrp_comm )
          !
          DO i = 1, 3
             !
             magtot_nc(i) = magtot_nc(i) * omega / (dffts%nr1*dffts%nr2*dffts%nr3 )
             !
          END DO
          !
!          absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
          write(stdout, '("change of the total magnetiztion in Bohr mag/cell")')
          write(stdout, '("mx",2f8.3)') magtot_nc(1)
          write(stdout, '("my",2f8.3)') magtot_nc(2)
          write(stdout, '("mz",2f8.3)') magtot_nc(3)
     end if
!  KC: end of the test

!      write(*,*)'do_elec, lgamma', do_elec, lgamma
end do    !iw

        tcpu = get_clock ('MAGNON')
        !
        DEALLOCATE (drhoscfs)
        IF (okvan) THEN
           DEALLOCATE (int3)
           IF (okpaw) DEALLOCATE (int3_paw)
           IF (noncolin) DEALLOCATE(int3_nc)
        ENDIF
  CALL stop_clock ('magscf')
  RETURN
END SUBROUTINE magscf
