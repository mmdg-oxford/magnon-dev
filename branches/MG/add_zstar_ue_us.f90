!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
subroutine add_zstar_ue_us(imode0,npe)
!----------===============-------------------------------
  ! add the contribution of the modes imode0+1 -> imode+npe
  ! to the effective charges Z(Us,E) (Us=scf,E=bare)
  !
  ! This subroutine is just for the USPP case
  !
  ! trans =.true. is needed for this calculation to be meaningful
  !

  USE kinds, ONLY : DP
  USE klist, ONLY : xk, wk
  USE uspp,  ONLY : nkb, vkb
  USE wvfct, ONLY : npwx, npw, nbnd, igk
  USE noncollin_module,   ONLY : npol
  USE wavefunctions_module,    ONLY : evc
  USE io_files, ONLY: iunigk

  USE qpoint,     ONLY : npwq, nksq
  USE efield_mod, ONLY: zstarue0_rec
  USE control_ph, ONLY : nbnd_occ
  USE eqv,        ONLY : dpsi, dvpsi
  USE modes,      ONLY : u
  USE units_ph,   ONLY : iucom, lrcom, iuwfc, lrwfc

  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum

  implicit none

  integer, intent(in) :: imode0, npe

  integer :: ik, jpol, nrec, mode, ipert, ibnd, jbnd, i,j

  real(DP) :: weight

  complex(DP), allocatable :: pdsp(:,:)
  complex(DP), allocatable :: dvkb(:,:,:)
  !  auxiliary space for <psi|ds/du|psi>

  call start_clock('add_zstar_us')
!  call compute_qdipol(dpqq)

  allocate (pdsp(nbnd,nbnd))
  allocate (dvkb(npwx,nkb,3))
  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq
     if (nksq.gt.1) read (iunigk) npw, igk
     npwq = npw
     weight = wk (ik)
     if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     call dvkb3(ik,dvkb)
     do ipert = 1, npe
        mode = imode0 + ipert
        do jpol = 1, 3
           dvpsi = (0.d0,0.d0)
           !
           ! read/compute the Commutator with the additional term
           call dvpsi_e(ik,jpol)
           !
           ! Calculate the matrix elements <psi_v'k|dS/du|psi_vk>
           ! Note: we need becp1
           !
           pdsp = (0.d0,0.d0)
           call psidspsi (ik, u (1, mode), pdsp )
#ifdef __MPI
           call mp_sum(pdsp, intra_pool_comm )
#endif
           !
           ! add the term of the double summation
           !
           do ibnd = 1, nbnd_occ(ik)
              do jbnd = 1, nbnd_occ(ik)
                 zstarue0_rec(mode,jpol)=zstarue0_rec(mode,jpol) +           &
                      weight *                                          &
                      dot_product(evc(1:npwx*npol,ibnd), &
                           dvpsi(1:npwx*npol,jbnd))*pdsp(jbnd,ibnd)
              enddo
           enddo

           dvpsi = (0.d0,0.d0)
           dpsi  = (0.d0,0.d0)
           !
           ! For the last part, we read the commutator from disc,
           ! but this time we calculate
           ! dS/du P_c [H-eS]|psi> + (dK(r)/du - dS/du)r|psi>
           !
           ! first we read  P_c [H-eS]|psi> and store it in dpsi
           !
           nrec = (jpol - 1) * nksq + ik
           call davcio (dpsi, lrcom, iucom, nrec, -1)
           !
           ! Apply the matrix dS/du, the result is stored in dvpsi
           !
           call add_for_charges(ik, u(1,mode))
           !
           ! Add  (dK(r)/du - dS/du) r | psi>
           !
           call add_dkmds(ik, u(1,mode),jpol, dvkb)
           !
           ! And calculate finally the scalar product
           !
           do ibnd = 1, nbnd_occ(ik)
              zstarue0_rec(mode,jpol)=zstarue0_rec(mode,jpol) - weight *   &
                   dot_product(evc(1:npwx*npol,ibnd),dvpsi(1:npwx*npol,ibnd))
           enddo
        enddo
     enddo
  enddo

  deallocate(dvkb)
  deallocate(pdsp)
  call stop_clock('add_zstar_us')

  return
end subroutine add_zstar_ue_us

