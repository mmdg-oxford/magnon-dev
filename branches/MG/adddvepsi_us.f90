!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE adddvepsi_us(becp1,becp2,ipol,kpoint,dvpsi)
  ! This subdoutine adds to dvpsi the terms which depend on the augmentation
  ! charge. It assumes that the variable dpqq, has been set and it is in
  ! the crystal basis.
  ! It calculates the last two terms of Eq.10 in JCP 21, 9934 (2004).
  ! P^+_c is applied in solve_e.
  !

  USE kinds, only : DP
  USE spin_orb, ONLY : lspinorb
  USE uspp,  ONLY : nkb, vkb, qq, qq_so
  USE wvfct, ONLY : npwx, npw, nbnd
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp_param, only: nh
  USE phus,     ONLY : dpqq, dpqq_so
  USE becmod,   ONLY : bec_type
  USE control_ph, ONLY: nbnd_occ
  USE control_flags, ONLY : gamma_only

  implicit none

  integer, intent(in) :: ipol, kpoint
  TYPE(bec_type), intent(in) :: becp1
  TYPE(bec_type), intent(in) :: becp2
  COMPLEX(KIND=DP), INTENT(INOUT) :: dvpsi(npwx*npol,nbnd)

  complex(DP), allocatable :: ps(:), ps_nc(:,:)
  REAL(KIND=DP), ALLOCATABLE :: ps_r(:)
  integer:: ijkb0, nt, na, ih, jh, ikb, jkb, ibnd, is, js, ijs

  IF (noncolin) THEN
     allocate (ps_nc(nbnd,npol))
  ELSEIF (gamma_only) THEN
     ALLOCATE (ps_r(nbnd))
  ELSE 
     allocate (ps(nbnd))
  END IF

  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp(na).eq.nt) then
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              IF (noncolin) THEN
                 ps_nc = (0.d0,0.d0)
              ELSEIF (gamma_only) THEN
                 ps_r = 0.0_DP
              ELSE
                 ps = (0.d0,0.d0)
              END IF
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 do ibnd=1, nbnd_occ(kpoint)
                    IF (noncolin) THEN
                       IF (lspinorb) THEN
                          ijs=0
                          do is=1,npol
                             do js=1,npol
                                ijs=ijs+1
                                ps_nc(ibnd,is)=ps_nc(ibnd,is) +          &
                                    qq_so(ih,jh,ijs,nt)*                 &
                                    (0.d0,1.d0)*becp2%nc(jkb,js,ibnd)       &
                                  + becp1%nc(jkb,js,ibnd)*        &
                                    dpqq_so(ih,jh,ijs,ipol,nt)
                              enddo
                           enddo
                       ELSE
                          DO is=1,npol
                             ps_nc(ibnd,is)=ps_nc(ibnd,is)+           &
                                qq(ih,jh,nt)*becp2%nc(jkb,is,ibnd)*(0.d0,1.d0) &
                               + dpqq(ih,jh,ipol,nt)*  &
                                 becp1%nc(jkb,is,ibnd)
                          END DO
                       END IF
                    ELSEIF (gamma_only) THEN
                       ps_r(ibnd) = ps_r(ibnd)+qq(ih,jh,nt)*becp2&
                            &%r(jkb,ibnd) + dpqq(ih,jh,ipol,nt)*&
                            & becp1%r(jkb,ibnd) 
                    ELSE
                       ps(ibnd) = ps(ibnd)+qq(ih,jh,nt)*becp2%k(jkb,ibnd) &
                           *(0.d0,1.d0) +  &
                            dpqq(ih,jh,ipol,nt)* becp1%k(jkb,ibnd)
                    END IF
                 enddo
              enddo
              do ibnd = 1, nbnd_occ (kpoint)
                 IF (noncolin) THEN
                    CALL zaxpy(npw,ps_nc(ibnd,1),vkb(1,ikb),1, &
                                                     dvpsi(1,ibnd),1)
                    CALL zaxpy(npw,ps_nc(ibnd,2),vkb(1,ikb),1, &
                                                     dvpsi(1+npwx,ibnd),1)
                 ELSEIF (gamma_only) THEN
                    CALL zaxpy(npw,CMPLX(ps_r(ibnd),0.0_DP,KIND=DP)&
                         &,vkb(1,ikb),1,dvpsi(1,ibnd),1) 
                 ELSE
                    CALL zaxpy(npw,ps(ibnd),vkb(1,ikb),1,dvpsi(1,ibnd),1)
                 END IF
              enddo
           enddo
           ijkb0=ijkb0+nh(nt)
        endif
     enddo
  enddo
  if (jkb.ne.nkb) call errore ('adddvepsi_us', 'unexpected error', 1)

  IF (noncolin) THEN
     deallocate(ps_nc)
  ELSEIF (gamma_only) THEN
     DEALLOCATE(ps_r)
  ELSE 
     deallocate(ps)
  END IF


  RETURN
END SUBROUTINE adddvepsi_us
