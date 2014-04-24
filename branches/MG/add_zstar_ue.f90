!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine add_zstar_ue (imode0, npe)
  !-----------------------------------------------------------------------
  ! add the contribution of the modes imode0+1 -> imode+npe
  ! to the effective charges Z(Us,E) (Us=scf,E=bare)
  !
  ! trans =.true. is needed for this calculation to be meaningful
  !
  USE kinds, only : DP
  USE klist, ONLY : xk, wk
  USE uspp,  ONLY : vkb
  USE wvfct, ONLY : npwx, npw, igk
  USE wavefunctions_module,  ONLY: evc
  USE noncollin_module,      ONLY: noncolin
  USE io_files, ONLY: iunigk
  USE qpoint,   ONLY: npwq, nksq
  USE eqv,      ONLY: dpsi, dvpsi
  USE efield_mod, ONLY: zstarue0_rec
  USE control_ph, ONLY : nbnd_occ
  USE units_ph,   ONLY : iudwf, lrdwf, iuwfc, lrwfc

  implicit none

  integer, intent(in) :: imode0, npe

  integer :: ibnd, jpol, ipert, nrec, mode, ik
  ! counter on bands
  ! counter on polarization
  ! counter on pertubations
  ! counter on records
  ! counter on modes
  ! counter on k points

  real(DP) :: weight

  complex(DP), external :: zdotc

  call start_clock('add_zstar_ue')
  zstarue0_rec=(0.0_DP,0.0_DP)
  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq
     if (nksq.gt.1) read (iunigk) npw, igk
     npwq = npw
     weight = wk (ik)
     if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     do jpol = 1, 3
        !
        ! read/compute DeltaV*psi(bare) for electric field
        !
        call dvpsi_e (ik, jpol)
        !
        do ipert = 1, npe
           mode = imode0 + ipert
           nrec = (ipert - 1) * nksq + ik
           !
           ! read dpsi(scf)/du for phonon mode # mode
           !

           call davcio (dpsi, lrdwf, iudwf, nrec, -1)
           do ibnd = 1, nbnd_occ(ik)
              zstarue0_rec (mode, jpol) = zstarue0_rec (mode, jpol) - 2.d0 * weight * &
                   zdotc (npw, dpsi (1, ibnd), 1, dvpsi (1, ibnd), 1)
              IF (noncolin) &
                 zstarue0_rec(mode,jpol)=zstarue0_rec (mode, jpol) - 2.d0 * weight * &
                   zdotc (npw, dpsi (1+npwx, ibnd), 1, dvpsi (1+npwx, ibnd), 1)

           enddo
        enddo
     enddo

  enddo
  call stop_clock('add_zstar_ue')
  return
end subroutine add_zstar_ue
