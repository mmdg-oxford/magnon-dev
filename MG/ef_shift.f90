! This file is copied and modified from QUANTUM ESPRESSO
! Kun Cao, Henry Lambert, Feliciano Giustino
 
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
subroutine ef_shift (drhoscf, ldos, ldoss, dos_ef, irr, npe, flag)
  !-----------------------------------------------------------------------
  !    This routine takes care of the effects of a shift of Ef, due to the
  !    perturbation, that can take place in a metal at q=0
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE wavefunctions_module, ONLY : evc
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : gg, nl
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npw, npwx, et
  USE klist,                ONLY : degauss, ngauss, ngk
  USE ener,                 ONLY : ef
  USE noncollin_module,     ONLY : nspin_mag, nspin_lsda
! modules from phcom
  USE qpoint,               ONLY : nksq
  USE control_ph,           ONLY : nbnd_occ, lgamma_gamma
  USE noncollin_module,     ONLY : noncolin, npol
  USE units_ph,             ONLY : lrwfc, iuwfc, lrdwf
  USE eqv,                  ONLY : dpsi
  USE modes,                ONLY : npert
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum

  implicit none
  !
  ! input/output variables
  !
  integer :: npe
  ! input: the number of perturbation

  complex(DP) :: drhoscf(dfftp%nnr, nspin_mag), &
       ldos(dfftp%nnr,nspin_mag), ldoss(dffts%nnr,nspin_mag)
  ! inp/out:the change of the charge
  ! inp: local DOS at Ef
  ! inp: local DOS at Ef without augme
  real(DP) :: dos_ef
  ! inp: density of states at Ef
  integer :: irr
  ! inp: index of the current irr. rep.
  logical :: flag
  ! inp: if true the eigenfunctions are updated
  !
  ! local variables
  !
  !--> these quantities may be complex since perturbation may be

  complex(DP) :: delta_n, wfshift, def(3)
  ! the change in electron number
  ! the shift coefficient for the wavefunction
  ! the change of the Fermi energy for each pert.
  ! NB: def(3) should be def (npertx) but then it cannot be saved
  !     anyway at Gamma the dimension of irreps never exceeds 3

  real(DP), external :: w0gauss
  ! the smeared delta function

  integer :: ibnd, ik, is, ipert, nrec, ikrec
  ! counter on occupied bands
  ! counter on k-point
  ! counter on spin polarizations
  ! counter on perturbations
  ! record number
  ! record position of wfc at k
  ! auxiliary for spin
  save def
  !
  ! determines Fermi energy shift (such that each pertubation is neutral)
  !
  call start_clock ('ef_shift')
  if (.not.flag) then
!     do ipert = 1, npert (irr)
     do ipert = 1, 1
        delta_n = (0.d0, 0.d0)
        do is = 1, nspin_lsda
           CALL fwfft ('Dense', drhoscf(:,is), dfftp)
           if (gg(1).lt.1.0d-8) delta_n = delta_n + omega*drhoscf(nl(1),is)
           CALL invfft ('Dense', drhoscf(:,is), dfftp)
        enddo
        call mp_sum ( delta_n, intra_pool_comm )
        def (ipert) = - delta_n / dos_ef
     enddo
     !
     ! symmetrizes the Fermi energy shift
     !
!     if (.not.lgamma_gamma) call sym_def (def, irr)
     WRITE( stdout, '(5x,"Pert. #",i3,": Fermi energy shift (Ry) =", 2e15.4)') &
 !         (ipert, def (ipert) , ipert = 1, npert (irr) )
          (ipert, def (ipert) , ipert = 1, 1 )
     !
     ! corrects the density response accordingly...
     !
!     do ipert = 1, npert (irr)
     do ipert = 1, 1
        call zaxpy (dfftp%nnr*nspin_mag, def(ipert), ldos, 1, drhoscf(1,1), 1)
     enddo
  else
     !
     ! does the same for perturbed wfc
     !
     do is = 1, nspin_mag
        call zaxpy (dffts%nnr, def(ipert), ldoss(1,is), 1, drhoscf(1,is), 1)
     enddo
  endif
  call stop_clock ('ef_shift')
  return
end subroutine ef_shift

