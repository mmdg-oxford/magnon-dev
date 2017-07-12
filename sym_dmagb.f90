! This file is copied and modified from QUANTUM ESPRESSO
! Kun Cao, Henry Lambert, Feliciano Giustino
 
!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine sym_dmagb (dvsym)
  !---------------------------------------------------------------------
  !
  !     This routine symmetrize the change of the potential due to an
  !     electric field perturbation. It is assumed that the perturbations
  !     are on the basis of the crystal
  !
  !
  USE kinds, only : DP
  USE cell_base,only : at, bg
  USE fft_base, only : dfftp
  USE symm_base,only : nsym, sname, s, ftau, t_rev, invs
  USE lsda_mod, only : nspin
  USE modes, only : nsymq
  implicit none

  complex(DP) :: dvsym (dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, nspin)
  complex(DP), allocatable ::  aux (:,:,:,:)
  complex(DP) ::  dmags(3), mag(3), magrot(3)
  ! the potential to symmetrize
  ! auxiliary quantity

  integer :: is, ri, rj, rk, i, j, k, irot, ipol, jpol, kpol
  ! counter on spin polarization
  ! the rotated points
  ! the point
  ! counter on symmetries
  ! counter on polarizations

!  do is = 2,4
!     do ipol = 1, 3
!        dvsym(:,:,:,is,ipol) = CMPLX(DBLE(dvsym(:,:,:,is,ipol)),0.d0,kind=DP)
!     end do
!  end do
  if (nsymq == 1) return
  allocate (aux(dfftp%nr1x , dfftp%nr2x , dfftp%nr3x , 3))

  do is = 2, 4
        aux(:,:,:,is-1) = dvsym(:,:,:,is)
        dvsym(:,:,:,is) = (0.d0, 0.d0)
  enddo
  !
  !  symmmetrize
  !
  do k = 1, dfftp%nr3
     do j = 1, dfftp%nr2
        do i = 1, dfftp%nr1
           do irot = 1, nsymq
              call ruotaijk (s(1,1,irot), ftau(1,irot), i, j, k, &
                             dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
              !
              !    ruotaijk find the rotated of i,j,k with the inverse of S
              !
!              dmags=(0.d0,0.d0)
!                 do is=1,3
!                       dmags(is)= aux(ri,rj,rk,is)
!                 enddo
!                 do kpol = 1, 3
!                    mag(kpol)=bg(1,kpol)*dmags(1) + &
!                              bg(2,kpol)*dmags(2) + &
!                              bg(3,kpol)*dmags(3)
!                 enddo
! rotate the magnetic moment
!                 do kpol = 1, 3
!                    magrot(kpol) = s(1,kpol,invs(irot))*mag(1) + &
!                                   s(2,kpol,invs(irot))*mag(2) + &
!                                   s(3,kpol,invs(irot))*mag(3)
!                 enddo
!                 if (sname(irot)(1:3)=='inv') magrot=-magrot
             IF(t_rev(irot).eq.1) THEN
                 dvsym(i,j,k,2) = dvsym(i,j,k,2) + conjg(aux(ri,rj,rk,1)) !+ mag(1)
                 dvsym(i,j,k,3) = dvsym(i,j,k,3) + conjg(aux(ri,rj,rk,2)) !+ mag(2)
                 dvsym(i,j,k,4) = dvsym(i,j,k,4) + conjg(aux(ri,rj,rk,3)) !+ mag(3)

             ELSE IF(0 .eq. 1)THEN
!                 dvsym(i,j,k,2) = dvsym(i,j,k,2) + aux(ri,rj,rk,1) !+ mag(1)
!                 dvsym(i,j,k,3) = dvsym(i,j,k,3) + aux(ri,rj,rk,2) !+ mag(2)
!                 dvsym(i,j,k,4) = dvsym(i,j,k,4) + aux(ri,rj,rk,3) !+ mag(3)

                 dmags=(0.d0,0.d0)
                 do is=1,3
                    dmags(is)= aux(ri,rj,rk,is)
                 enddo
                 do kpol = 1, 3
                    mag(kpol)=bg(1,kpol)*dmags(1) + &
                              bg(2,kpol)*dmags(2) + &
                              bg(3,kpol)*dmags(3)
                 enddo
! rotate the magnetic moment
                 do kpol = 1, 3
                    magrot(kpol) = s(1,kpol,irot)*mag(1) + &
                                   s(2,kpol,irot)*mag(2) + &
                                   s(3,kpol,irot)*mag(3)
                 enddo
                 magrot=-magrot
                 if (sname(irot)(1:3)=='inv') magrot=-magrot
! go back to cartesian coordinates
                 do kpol = 1, 3
                    mag(kpol)=at(kpol,1)*magrot(1) + &
                              at(kpol,2)*magrot(2) + &
                              at(kpol,3)*magrot(3)
                 enddo

                 dvsym(i,j,k,2) = dvsym(i,j,k,2) + conjg(mag(1))
                 dvsym(i,j,k,3) = dvsym(i,j,k,3) + conjg(mag(2))
                 dvsym(i,j,k,4) = dvsym(i,j,k,4) + conjg(mag(3))

! go back to cartesian coordinates
!                 do kpol = 1, 3
!                    mag(kpol)=at(kpol,1)*magrot(1) + &
!                              at(kpol,2)*magrot(2) + &
!                              at(kpol,3)*magrot(3)
!                 enddo
              ELSE
                 dvsym(i,j,k,2) = dvsym(i,j,k,2) + aux(ri,rj,rk,1) !+ mag(1)
                 dvsym(i,j,k,3) = dvsym(i,j,k,3) + aux(ri,rj,rk,2) !+ mag(2)
                 dvsym(i,j,k,4) = dvsym(i,j,k,4) + aux(ri,rj,rk,3) !+ mag(3)
              END IF
           enddo
        enddo
     enddo
  enddo
  do is=2,4
!     do ipol = 1, 3
        dvsym(:,:,:,is) = dvsym(:,:,:,is) / DBLE(nsymq)
!     enddo
  enddo
  deallocate (aux)
  return
end subroutine sym_dmagb
