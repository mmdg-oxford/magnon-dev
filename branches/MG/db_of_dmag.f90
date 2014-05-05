subroutine db_of_dmag (dvscf, add_nlcc)
  USE kinds,     ONLY : DP
  USE constants, ONLY : e2, fpi
  USE fft_base,  ONLY: dfftp
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvect,     ONLY : nl, ngm, g,nlm
  USE cell_base, ONLY : alat, tpiba2
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag, nspin_gga
  USE funct,     ONLY : dft_is_gradient
  USE scf,       ONLY : rho, rho_core
  USE eqv,       ONLY : dmuxc
  USE nlcc_ph,   ONLY : nlcc_any
  USE qpoint,    ONLY : xq
  USE gc_ph,     ONLY : grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  USE control_ph, ONLY : lrpa
  USE control_flags, only : gamma_only
  !OBM: gamma_only is disregarded for phonon calculations, TDDFPT purposes only

  implicit none

  complex(DP), intent(inout):: dvscf (dfftp%nnr, nspin_mag)
  complex(DP), intent(inout):: db_ext (dfftp%nnr, 3)

  integer :: ipol


  do ipol=1, 3
    call fwfft('Dense', dvscf(:, 1+ipol), dfftp)
    do ig = 1, ngm 
      dbg  = dvscf(nl(ig),1+ipol)
      qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
      db_ext(nl(ig), ipol) = -(0.50/3.0)( )

    enddo

    call fwfft('Dense', dvscf(:, 1+ipol), dfftp)
    call invfft('Dense', dvscf(:, 1+ipol), dfftp)
  enddo

end subroutine db_of_dmag

