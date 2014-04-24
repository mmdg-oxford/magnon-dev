!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------
subroutine zstar_eu_us
!----------===========-----------------------------------------
  !
  ! Calculates the additional part of the Born effective charges
  ! in the case of USPP
  !
  !
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE klist,     ONLY : xk, wk
  USE gvecs,   ONLY : doublegrid
  USE fft_base, ONLY : dfftp, dffts
  USE lsda_mod,  ONLY : nspin, current_spin, isk, lsda
  USE io_files,  ONLY : iunigk
  USE uspp,      ONLY : okvan, nkb, vkb
  USE wvfct,     ONLY : nbnd, npw, npwx, igk
  USE paw_variables, ONLY : okpaw
  USE wavefunctions_module,    ONLY : evc
  USE uspp_param,       ONLY : upf, nhm, nh
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE qpoint,    ONLY : nksq, npwq
  USE control_ph, ONLY : nbnd_occ
  USE efield_mod, ONLY : zstareu0
  USE eqv,        ONLY : dvpsi, dpsi
  USE phus,       ONLY : becsumort, int3, int3_paw
  USE modes,      ONLY : u, npert, nirr
  USE nlcc_ph,    ONLY : nlcc_any
  USE units_ph,   ONLY : lrdwf, iucom, lrcom, lrebar, iuebar, lrdrhous, &
                         iudrhous, iudwf, lrwfc, iuwfc
  USE mp_global, ONLY : nproc_pool, npool

  !
  implicit none
  integer :: ibnd, jbnd, ipol, jpol, imode0, irr, imode, nrec, mode
  integer :: ik,  ig, ir, is, i, j, mu, ipert
  integer :: ih, jh, ijh
  integer :: iuhxc, lrhxc
  !
  real(DP) :: weight, fact
  !
  complex(DP), allocatable :: dbecsum(:,:,:,:), aux1 (:)
  COMPLEX(DP), ALLOCATABLE :: dbecsum_nc(:, :, :, :, :)
  ! the becsum with dpsi
  ! auxillary work space for fft
  complex(DP) , pointer ::      &
      dvscf(:,:,:)
  complex(DP), allocatable :: pdsp(:,:)
  complex(DP), allocatable :: drhoscfh (:,:)
  complex(DP), allocatable :: dvkb (:,:,:)
  integer :: npe, irr1, imode1, na, nt

#ifdef TIMINIG_ZSTAR_US
  call start_clock('zstar_eu_us')
  call start_clock('zstar_us_1')
#endif

  !  auxiliary space for <psi|ds/du|psi>
  allocate (dvscf( dfftp%nnr , nspin_mag, 3))
  allocate (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, 3))
  if (noncolin) allocate (dbecsum_nc( nhm, nhm, nat, nspin, 3))
  allocate (aux1(  dffts%nnr))
  allocate (pdsp(nbnd,nbnd))

  !
  ! Set the initial values to zero
  !
  pdsp    = (0.d0,0.d0)
  dvscf   = (0.d0,0.d0)
  dbecsum = (0.d0,0.d0)
  if (noncolin) dbecsum_nc=(0.d0,0.d0)
  !
  ! first we calculate the perturbed charge density and the perturbed
  ! Hartree and exchange and correlation potential , which we need later
  ! for the calculation of the Hartree and xc part
  !
  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq
     if (nksq.gt.1) read (iunigk) npw, igk
     npwq = npw
     if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
     if (lsda) current_spin = isk (ik)
     call init_us_2 (npw, igk, xk(1,ik), vkb)
     weight = wk (ik)
     do jpol = 1, 3
        nrec = (jpol - 1) * nksq + ik
        call davcio (dpsi, lrdwf, iudwf, nrec, - 1)
        if (noncolin) then
           call incdrhoscf_nc (dvscf(1,1,jpol),weight,ik, &
                              dbecsum_nc(1,1,1,1,jpol), dpsi)
        else
           call incdrhoscf (dvscf(1,current_spin,jpol),weight,ik, &
                            dbecsum(1,1,current_spin,jpol), dpsi)
        endif
     end do
  end do

#ifdef __MPI
     IF (noncolin) THEN
        call mp_sum ( dbecsum_nc, intra_pool_comm )
     ELSE
        call mp_sum ( dbecsum, intra_pool_comm )
     END IF
#endif
#ifdef TIMINIG_ZSTAR_US
  call stop_clock('zstar_us_1')
  call start_clock('zstar_us_2')
#endif

  if (doublegrid) then
     do is = 1, nspin_mag
        do ipol = 1, 3
           call cinterpolate(dvscf(1,is,ipol),dvscf(1,is,ipol), 1)
        end do
     end do
  end if

  IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 3)

  call addusddense (dvscf, dbecsum)

#ifdef __MPI
  call mp_sum ( dvscf, inter_pool_comm )
#endif

#ifdef TIMINIG_ZSTAR_US
  call stop_clock('zstar_us_2')
  call start_clock('zstar_us_3')
#endif

  if (nlcc_any) call addnlcc_zstar_eu_us (dvscf)

  do ipol = 1, 3
     !
     ! Instead of recalculating the perturbed charge density,
     ! it can also be read from file
     ! NB: Then the variable fildrho must be set
     !
     ! call davcio_drho(dvscf(1,1,ipol),lrdrho,iudrho,ipol,-1)
     !
     call dv_of_drho (0, dvscf (1, 1, ipol), .false.)
  enddo
#ifdef __MPI
  call psyme (dvscf)
#else
  call syme (dvscf)
#endif

#ifdef TIMINIG_ZSTAR_US
  call stop_clock('zstar_us_3')
  call start_clock('zstar_us_4')
#endif
!
! Calculate the parts with the perturbed Hartree and exchange and correlation
! potenial
!
  imode0 = 0
  allocate(drhoscfh(dfftp%nnr,nspin_mag))
  do irr = 1, nirr
     npe = npert(irr)
     do imode = 1, npe
        mode = imode0 + imode
        call davcio (drhoscfh, lrdrhous, iudrhous, mode, -1)
        do jpol = 1, 3
           do is=1,nspin_mag
              zstareu0(jpol,mode) =  zstareu0(jpol,mode) -                  &
                 dot_product(dvscf(1:dfftp%nnr,is,jpol),drhoscfh(1:dfftp%nnr,is)) &
               * omega / DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
           end do
        end do
     end do
     imode0 = imode0 + npe
  end do
  deallocate (drhoscfh)
#ifdef TIMINIG_ZSTAR_US
  call stop_clock('zstar_us_4')
  call start_clock('zstar_us_5')
#endif
  !
  !  Calculate the part with the position operator
  !
  allocate (dvkb(npwx,nkb,3))
  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq
     if (nksq.gt.1) read (iunigk) npw, igk
     npwq = npw
     weight = wk (ik)
     if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     call dvkb3(ik, dvkb)
     imode0 = 0
     do irr = 1, nirr
        do imode = 1, npert (irr)
           mode = imode+imode0
           do jpol = 1, 3
              dvpsi = (0.d0,0.d0)
              !
              ! read the Commutator+add. terms
              !
              nrec = (jpol - 1) * nksq + ik
              call davcio (dvpsi, lrebar, iuebar, nrec, - 1)
              !
              pdsp = (0.d0,0.d0)
              call psidspsi (ik, u (1, mode), pdsp )
#ifdef __MPI
              call mp_sum( pdsp, intra_pool_comm )
#endif
              !
              ! add the term of the double summation
              !
              do ibnd = 1, nbnd_occ(ik)
                 do jbnd = 1, nbnd_occ(ik)
                    zstareu0(jpol,mode)=zstareu0(jpol, mode) +           &
                         weight *                                        &
                         dot_product(evc(1:npwx*npol,ibnd),              &
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
              ! Apply the matrix dS/du
              !
              call add_for_charges(ik, u(1,mode))
              !
              ! Add  (dK(r)/du - dS/du) r | psi>
              !
              call add_dkmds(ik, u(1,mode), jpol, dvkb)
              !
              ! And calculate finally the scalar product
              !
              do ibnd = 1, nbnd_occ(ik)
                 zstareu0(jpol,mode)=zstareu0(jpol, mode) - weight *  &
                      dot_product(evc(1:npwx*npol,ibnd),dvpsi(1:npwx*npol,ibnd))
              enddo
           enddo
        enddo
        imode0 = imode0 + npert (irr)
     enddo
  enddo
  deallocate (dvkb)

  deallocate (pdsp)
  deallocate (dbecsum)
  if (noncolin) deallocate(dbecsum_nc)
  deallocate (dvscf)
  deallocate (aux1)

  fact=1.0_DP
#ifdef __MPI
  fact=1.0_DP/nproc_pool/npool
#endif
  IF (okpaw) THEN
     imode0 = 0
     do irr1 = 1, nirr
        do ipert = 1, npert (irr1)
           mode = imode0 + ipert
           do nt=1,ntyp
              if (upf(nt)%tpawp) then
                 ijh=0
                 do ih=1,nh(nt)
                    do jh=ih,nh(nt)
                       ijh=ijh+1
                       do na=1,nat
                          if (ityp(na)==nt) then
                             do jpol = 1, 3
                                do is=1,nspin_mag
                                 zstareu0(jpol,mode)=zstareu0(jpol,mode)  &
                                    -fact*int3_paw(ih,jh,jpol,na,is)* &
                                            becsumort(ijh,na,is,mode)
                                enddo
                             enddo
                          endif
                       enddo
                    enddo
                 enddo
              endif
           enddo
        enddo
        imode0 = imode0 + npert (irr1)
     enddo
  endif

#ifdef TIMINIG_ZSTAR_US
  call stop_clock('zstar_us_5')
  call stop_clock('zstar_eu_us')
#endif

  return
end subroutine zstar_eu_us
