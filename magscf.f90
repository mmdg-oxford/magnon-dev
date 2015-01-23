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
  USE uspp,  ONLY: okvan
  USE efield_mod, ONLY : zstarue0, zstarue0_rec
  USE control_ph, ONLY : zue, convt, rec_code, do_elec!, do_trans
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

  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum
  USE freq_ph,       ONLY : fpol, fiu, nfs, nfsmax
  USE gvecs,         ONLY : doublegrid

  IMPLICIT NONE

  INTEGER :: irr, irr1, imode0, npe, ig,iw
  ! counter on the representations
  ! counter on the representations
  ! counter on the modes
  ! npert(irr)

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
        WRITE( stdout, *)'npol,nspin_mag', npol, nspin_mag
        WRITE( stdout, *)'dfftp,dffts, doublegrid', dfftp%nnr, dffts%nnr, doublegrid

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
!                write(stdout,'(7f14.7)') (real(drhoscfs (ig,1)), ig = 1,7)
                do ig=1, nspin_mag
                   CALL fwfft ('Smooth', drhoscfs(:,ig), dffts)
                enddo
!                WRITE(stdout, *)
!                WRITE(stdout, '("magnetization density response" )')
!                WRITE(stdout, *)
!               write(stdout,'(7f14.7)') (real(drhoscfs (ig,2)), ig = 1,7)
!                write(stdout,'(7f14.7)') (real(drhoscfs (ig,3)), ig = 1,7)
!                write(stdout,'(7f14.7)') (real(drhoscfs (ig,4)), ig = 1,7)
!        enddo


!        WRITE(stdout, '("magnetization density response" )')
      if(do_elec) then
        write(stdout,*)'charge density response'
        write(stdout,'("eps^{-1}, "3f12.5,"  ",f14.7)') xq(:), (1.0d0 + real(drhoscfs (1,1)))
        write(stdout,'("eps, "3f12.5,"  ",f14.7)') xq(:), (1.0/(1.0d0 + real(drhoscfs (1,1))))
        write(stdout,'("real(eps) im(eps) "4f12.5,"  ",2f14.7)') xq(:),real(fiu(iw)), real((1.0/(1.0d0 + drhoscfs (1,1)))), &
                     aimag((1.0/(1.0d0 + drhoscfs (1,1))))
      end if

      IF (noncolin) THEN
        write(stdout,*)"density matrix response with G=G'=(0 0 0 )"
!       G=G'= (0 0 0)
        write(stdout,'("real drho, "3f12.5,"  ",4f14.7)') xq(:), (real(drhoscfs (1,ig)), ig = 1,4)
        write(stdout,'("im drho", 3f12.5,"  ",4f14.7)') xq(:), (aimag(drhoscfs (1,ig)), ig = 1,4)
!      endif

!      if(do_trans) then
!      write(stdout,*)'output2'
        WRITE(stdout, '("transverse magnetic response" )')
        write(stdout,'("w, chiq+-, "f12.5,"  ",2f14.7)') real(fiu(iw))*13600, real(drhoscfs(1,3)+drhoscfs(1,2)) &
            -aimag(drhoscfs(1,3)-drhoscfs(1,2)), aimag(drhoscfs(1,3)+drhoscfs(1,2))+real(drhoscfs(1,3)-drhoscfs(1,2))
      END IF                                                      
!      write(stdout,'("imchiq+-", 3f12.5,"  ",f14.7)') xq(:),drhoscfs(1,3)-drhoscfs(1,2)
!      end if
end do
        tcpu = get_clock ('MAGNON')
        !
        DEALLOCATE (drhoscfs)
  CALL stop_clock ('magscf')
  RETURN
END SUBROUTINE magscf
