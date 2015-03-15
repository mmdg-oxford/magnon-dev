!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvqpsi_mag_us (ik, addnlcc)
  !----------------------------------------------------------------------
  !
  ! This routine applies the external magnetic field to 
  ! the wave functions in the first iteration of the self-consistent process.
  !
  !
  USE kinds, only : DP
  USE ions_base, ONLY : nat, ityp
  USE cell_base, ONLY : tpiba
  USE fft_base,   ONLY: dfftp, dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvect,     ONLY : eigts1, eigts2, eigts3, mill, g, nl, &
                        ngm
  USE gvecs,   ONLY : ngms, doublegrid, nls
  USE lsda_mod,  ONLY : lsda, isk
  USE noncollin_module, ONLY : npol, noncolin
  use uspp_param,ONLY : upf
  USE wvfct,     ONLY : nbnd, npw, npwx, igk
  USE wavefunctions_module,  ONLY: evc
  USE nlcc_ph,    ONLY : nlcc_any, drc
  USE eqv,        ONLY : dvpsi, dmuxc, vlocq
  USE qpoint,     ONLY : npwq, igkq, xq, eigqts, ikks, dbext
  USE control_ph, ONLY : do_elec

  implicit none
  !
  !   The dummy variables
  !

  integer :: ik
  ! input: the k point
  complex(DP) :: uact (3 * nat)
  ! input: the pattern of displacements
  logical :: addnlcc
  !
  !   And the local variables
  !

  integer :: na, mu, ikk, ig, nt, ibnd, ir, is, ip
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! counter on G vectors
  ! the type of atom
  ! counter on bands
  ! counter on real mesh

  complex(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  complex(DP) , allocatable, target :: aux (:)
  complex(DP),  allocatable  :: psic (:,:)
  complex(DP) , allocatable :: aux1 (:), aux2 (:,:)
  complex(DP) , pointer :: auxs (:)
  ! work space

  call start_clock ('dvqpsi_us')
  if (nlcc_any.and.addnlcc) then
     allocate (aux( dfftp%nnr))
     if (doublegrid) then
        allocate (auxs(dffts%nnr))
     else
        auxs => aux
     endif
  endif

  allocate (aux1(dffts%nnr))
  allocate (aux2(dffts%nnr, npol))
  allocate (psic(dffts%nnr, npol))
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  ikk = ikks(ik)
  dvpsi(:,:) = (0.d0, 0.d0)
  aux1(:)    = (0.d0, 0.d0)
  aux2(:,:)  = (0.d0, 0.d0)

  !in a.u. \mu_{B} = 1/2
  ! \mu_{B} * FFT[B_{q}(G)]
  !scalar component of field.
  !just set to G= 0 0 0 
  !aux1 (nls(1)) = 0.5*4.255d-6
  !aux1 (nls(1)) = 1.0
  !for electric field
  aux1 (nls(1)) = dcmplx(1.d0, 0.0d0)
  !moved into real space:
  CALL invfft ('Smooth', aux1, dffts)
  !
  ! add NLCC when present
  !
  if (nlcc_any.and.addnlcc) then
     WRITE(6,'("WARNING NLCC NOT IMPLEMENTED.")')
  endif

  do ibnd = 1, nbnd
       psic = (0.d0, 0.d0)
       aux2 = (0.d0, 0.d0)
     
       do ip=1,npol
          if(ip==1)then
              do ig = 1, npw
              psic (nls (igk (ig)), 1 ) = evc (ig, ibnd)
              end do
           else
              do ig = 1, npw
              psic (nls (igk (ig)), 2 ) = evc (ig+npwx, ibnd)
              end do
           end if
        end do
!    enddo

    do ip=1,npol
       if(ip==1)then 
       CALL invfft ('Wave', psic(:,1), dffts)
       else
       CALL invfft ('Wave', psic(:,2), dffts)
       end if
     end do

!For electric field perturbation couples directly to wave functions.
   if(do_elec) then
      do ip=1,npol
        if(ip==1)then 
           do ir = 1, dffts%nnr
             aux2(ir,1) = aux2(ir,1) + aux1 (ir)*psic(ir,1)
           end do
        else
           do ir = 1, dffts%nnr
             aux2(ir,2) = aux2(ir,2) + aux1 (ir)*psic(ir,2)
           end do
        end if
      end do 
   else  if(noncolin)then
!HL initial test with B_{+-}

!  if(noncolin)then
!X: \mu_b \sigma_{x}B_{x}
    if(dbext(1).ne.0.d0) then
        do ir = 1, dffts%nnr
          aux2(ir,1) = aux2(ir,1) + aux1 (ir)*dbext(1)*psic(ir,2)
          aux2(ir,2) = aux2(ir,2) + aux1 (ir)*dbext(1)*psic(ir,1)
        enddo
    endif

!Y: \mu_b \sigma_{y}B_{y}
    if(dbext(2).ne.0.d0) then
        do ir = 1, dffts%nnr
           aux2(ir,1) = aux2(ir,1) + (0.0d0, -1.0d0)*aux1(ir)*dbext(2)*(0.d0,1.d0)*psic(ir,2)
           aux2(ir,2) = aux2(ir,2) + (0.0d0, 1.0d0)*aux1(ir)*dbext(2)*(0.d0,1.d0)*psic(ir,1)
        enddo
    endif

!Z: \mu_b \sigma_{z}B_{z}
    if(dbext(3).ne.0.d0) then
        do ir = 1, dffts%nnr
           aux2(ir,1) = aux2(ir,1) + aux1 (ir)*dbext(3)*psic(ir,1)
           aux2(ir,2) = aux2(ir,2) - aux1 (ir)*dbext(3)*psic(ir,2)
        enddo
    endif
   

   endif

   do ip=1,npol
      if(ip==1)then
      CALL fwfft ('Wave', aux2(:,1), dffts)
      else
      CALL fwfft ('Wave', aux2(:,2), dffts)
      end if
   enddo

      do ip = 1, npol
        if (ip==1) then
           do ig = 1, npwq
              dvpsi (ig, ibnd) = aux2 (nls (igkq (ig) ), 1)
           enddo
        else
           do ig = 1, npwq
             dvpsi (ig+npwx, ibnd) = aux2 (nls (igkq (ig) ), 2)
           enddo
        end if
      enddo!npol

  enddo!nbnd
!
  deallocate (psic)
  deallocate (aux1)
  deallocate (aux2)
  if (nlcc_any.and.addnlcc) then
     deallocate (aux)
     if (doublegrid) deallocate (auxs)
  endif
  !
  !HL additional term from contribution nonlocal potential
  !call dvqpsi_us_only (ik, uact)
  call stop_clock ('dvqpsi_us')
  return
end subroutine dvqpsi_mag_us
