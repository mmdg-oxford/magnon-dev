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
  USE fft_base,   ONLY : dfftp
  USE uspp,  ONLY: okvan
  USE efield_mod, ONLY : zstarue0, zstarue0_rec
  USE control_ph, ONLY : zue, convt, rec_code
  USE partial,    ONLY : done_irr, comp_irr
  USE modes,      ONLY : nirr, npert, npertx
  USE phus,       ONLY : int3, int3_nc, int3_paw
  USE uspp_param, ONLY : nhm
  USE eqv,        ONLY : drhoscfs
  USE paw_variables, ONLY : okpaw
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE recover_mod, ONLY : write_rec
  USE qpoint,          ONLY : xq
  USE fft_base,   ONLY: dfftp, dffts
  USE fft_interfaces, ONLY: fwfft, invfft

  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum

  IMPLICIT NONE

  INTEGER :: irr, irr1, imode0, npe, ig
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

        ALLOCATE (drhoscfs( dfftp%nnr , nspin_mag))
        imode0 = 0

        WRITE( stdout, '(/,5x,"qpoint= ", 3f12.5)'), xq(1:3)
        WRITE( stdout, '(/,5x,"Self-consistent Calculation")')

        CALL solve_linter (drhoscfs(1,1))

        
      
        WRITE( stdout, '(/,5x,"End of self-consistent calculation")')
        WRITE( stdout, '(/,5x,"qpoint= ", 3f12.5)'), xq(1:3)
        WRITE( stdout, '(/,5x,"X_[G](Gp)")')
        WRITE( stdout, '("charge density response ")')
        WRITE( stdout, *)


        write(stdout,'(7f14.7)') (real(drhoscfs (ig,1)), ig = 1,7)

        do ig=1, nspin_mag
          CALL fwfft ('Smooth', drhoscfs(:,ig), dffts)
        enddo

        WRITE( stdout, *)
        WRITE(stdout, '("magnetization density response" )')
        WRITE(stdout, *)
        write(stdout,'(7f14.7)') (real(drhoscfs (ig,2)), ig = 1,7)
        write(stdout,'(7f14.7)') (real(drhoscfs (ig,3)), ig = 1,7)
        write(stdout,'(7f14.7)') (real(drhoscfs (ig,4)), ig = 1,7)


        WRITE(stdout, '("magnetization density response" )')
        write(stdout,'("eps^{-1}, "3f12.5,"  ",4f14.7)') xq(:), (1.0d0 + real(drhoscfs (1,1)))
        write(stdout,'("eps, "3f12.5,"  ",4f14.7)') xq(:), (1.0/(1.0d0 + real(drhoscfs (1,1))))
!       G=G'= (0 0 0)
        write(stdout,'("rechiq, "3f12.5,"  ",4f14.7)') xq(:), (real(drhoscfs (1,ig)), ig = 1,4)
        write(stdout,'("imchiq", 3f12.5,"  ",4f14.7)') xq(:), (aimag(drhoscfs (1,ig)), ig = 1,4)
!        write(600,'("rechiq, "3f12.5,"  ",4f14.7)') xq(:), (real(drhoscfs (1,ig)), ig = 1,4)
!        write(601,'("imchiq", 3f12.5,"  ",4f14.7)') xq(:), (aimag(drhoscfs (1,ig)), ig = 1,4)

        tcpu = get_clock ('MAGNON')
        !
        DEALLOCATE (drhoscfs)
  CALL stop_clock ('magscf')
  RETURN
END SUBROUTINE magscf
