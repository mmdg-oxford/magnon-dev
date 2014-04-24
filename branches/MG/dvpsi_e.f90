!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvpsi_e (ik, ipol)
  !----------------------------------------------------------------------
  !
  ! On output: dvpsi contains P_c^+ x | psi_ik > in crystal axis
  !            (projected on at(*,ipol) )
  !
  ! dvpsi is READ from file if this_pcxpsi_is_on_file(ik,ipol)=.true.
  ! otherwise dvpsi is COMPUTED and WRITTEN on file (vkb,evc,igk must be set)
  !
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : tpiba2
  USE io_global,       ONLY : stdout
  USE klist,           ONLY : xk
  USE gvect,           ONLY : g
  USE wvfct,           ONLY : npw, npwx, nbnd, igk, g2kin, et
  USE wavefunctions_module, ONLY: evc
  USE noncollin_module,ONLY : noncolin, npol
  USE becmod,          ONLY : bec_type, becp, calbec, &
                              allocate_bec_type, deallocate_bec_type
  USE uspp,            ONLY : okvan, nkb, vkb
  USE uspp_param,      ONLY : nh, nhm
  USE ramanm,          ONLY : eth_rps
  USE eqv,             ONLY : dpsi, dvpsi, eprec
  USE phus,            ONLY : becp1
  USE qpoint,          ONLY : nksq, npwq
  USE units_ph,        ONLY : this_pcxpsi_is_on_file, lrcom, iucom, &
                              lrebar, iuebar
  USE control_ph,      ONLY : nbnd_occ

  implicit none
  !
  integer, intent(IN) :: ipol, ik
  !
  ! Local variables
  !
  integer :: ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0,  &
             nrec, is, js, ijs
  ! counters

  real(DP), allocatable  :: h_diag (:,:)
  ! the diagonal part of h_scf
  type(bec_type) :: becp2 ! the scalar products
  complex(DP), allocatable :: spsi(:,:)
  real(DP) ::   anorm, thresh
  ! preconditioning cut-off
  ! the desired convergence of linter
  logical :: conv_root
  ! true if convergence has been achieved

  external ch_psi_all, cg_psi
  !
  call start_clock ('dvpsi_e')
  dpsi=(0.d0, 0.d0)
  dvpsi=(0.d0, 0.d0)
  if (this_pcxpsi_is_on_file(ik,ipol)) then
     nrec = (ipol - 1)*nksq + ik
     call davcio(dvpsi, lrebar, iuebar, nrec, -1)
     call stop_clock ('dvpsi_e')
     return
  end if
  !
  call allocate_bec_type ( nkb, nbnd, becp2)

  ! calculate the commutator [H,x_ipol]  psi > and store it in dpsi
  ! dvpsi used as workspace
  call commutator_Hx_psi (ik, nbnd_occ(ik), becp1(ik), becp2, ipol, dpsi, dvpsi )
  !
  !    orthogonalize dpsi to the valence subspace: ps = <evc|dpsi>
  !    Apply -P^+_c
  !    NB it uses dvpsi as workspace
  !
  CALL orthogonalize(dpsi, evc, ik, ik, dvpsi, npwq)
  dpsi=-dpsi
  !
  !   dpsi contains P^+_c [H-eS,x] psi_v for the three crystal polarizations
  !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
  !

  do ig = 1, npw
     g2kin (ig) = SUM((xk(1:3,ik) +g (1:3, igk (ig)) ) **2) *tpiba2
  enddo
  allocate (h_diag( npwx*npol, nbnd))
  h_diag=0.d0
  do ibnd = 1, nbnd_occ (ik)
     do ig = 1, npw
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd,ik) )
     enddo
     IF (noncolin) THEN
        do ig = 1, npw
           h_diag (ig+npwx, ibnd) = 1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
        enddo
     END IF
  enddo
  !
  dvpsi(:,:) = (0.d0, 0.d0)
  !
  thresh = eth_rps
  call cgsolve_all (ch_psi_all, cg_psi, et (1, ik), dpsi, dvpsi, &
       h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, &
       nbnd_occ (ik), npol)

  if (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
       & " linter: root not converged ",e10.3)') &
       ik, ibnd, anorm
  !
  CALL flush_unit( stdout )
  deallocate (h_diag)
  !
  ! we have now obtained P_c x |psi>.
  ! In the case of USPP this quantity is needed for the Born
  ! effective charges, so we save it to disc
  !
  ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
  ! therefore we apply S again, and then subtract the additional term
  ! furthermore we add the term due to dipole of the augmentation charges.
  !
  if (okvan) then
     !
     ! for effective charges
     !
     nrec = (ipol - 1) * nksq + ik
     call davcio (dvpsi, lrcom, iucom, nrec, 1)
     !
     allocate (spsi ( npwx*npol, nbnd))
     CALL calbec (npw, vkb, dvpsi, becp )
     CALL s_psi(npwx,npw,nbnd,dvpsi,spsi)
     call dcopy(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
     deallocate (spsi)
     CALL adddvepsi_us(becp1(ik),becp2,ipol,ik,dvpsi)
  endif

  IF (nkb > 0) call deallocate_bec_type (becp2)

  nrec = (ipol - 1)*nksq + ik
  call davcio(dvpsi, lrebar, iuebar, nrec, 1)
  this_pcxpsi_is_on_file(ik,ipol) = .true.
  call stop_clock ('dvpsi_e')
  return
end subroutine dvpsi_e
