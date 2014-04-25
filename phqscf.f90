!
! Copyright (C) 2001-2008 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE phqscf
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

  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum

  IMPLICIT NONE

  INTEGER :: irr, irr1, imode0, npe
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

  CALL start_clock ('phqscf')
  !
  !    For each irreducible representation we compute the change
  !    of the wavefunctions
  !
        ALLOCATE (drhoscfs( dfftp%nnr , nspin_mag, npe))
        imode0 = 0
        IF (okvan) THEN
           ALLOCATE (int3 ( nhm, nhm, npe, nat, nspin_mag))
           IF (okpaw) ALLOCATE (int3_paw (nhm, nhm, npe, nat, nspin_mag))
           IF (noncolin) ALLOCATE(int3_nc( nhm, nhm, npe, nat, nspin))
        ENDIF


        WRITE( stdout, '(/,5x,"Self-consistent Calculation")')
        CALL solve_linter (drhoscfs)
        WRITE( stdout, '(/,5x,"End of self-consistent calculation")')
        !
        !   Add the contribution of this mode to the dynamical matrix
        !
        IF (okvan) THEN
           DEALLOCATE (int3)
           IF (okpaw) DEALLOCATE (int3_paw)
           IF (noncolin) DEALLOCATE(int3_nc)
        ENDIF
        tcpu = get_clock ('PHONON')
        !
        DEALLOCATE (drhoscfs)
  CALL stop_clock ('phqscf')
  RETURN
END SUBROUTINE phqscf
