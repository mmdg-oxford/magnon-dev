!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE psyme (dvtosym)
  !-----------------------------------------------------------------------
  !
  ! ...  p-symmetrize the charge density.
  !
  USE kinds,     ONLY : DP
  USE fft_base, ONLY : dfftp
  USE noncollin_module, ONLY : nspin_mag
  USE mp_global, ONLY : me_pool
  USE fft_base,  ONLY : dfftp, cgather_sym
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: dvtosym (dfftp%nnr, nspin_mag, 3)
    ! the potential to symmetrize
    !-local variable
  !
#if defined (__MPI)
  !
  INTEGER :: i, is, iper, npp0
  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:,:)
    ! the potential to symmet
  !
  !
  ALLOCATE (ddvtosym ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin_mag, 3))
  npp0 = 0
  DO i = 1, me_pool
     npp0 = npp0 + dfftp%npp (i)
  ENDDO

  npp0 = npp0 * dfftp%nnp+1
  DO iper = 1, 3
     DO is = 1, nspin_mag
        CALL cgather_sym (dvtosym (:, is, iper), ddvtosym (:, is, iper) )
     ENDDO

  ENDDO

  CALL syme (ddvtosym)
  DO iper = 1, 3
     DO is = 1, nspin_mag
        CALL zcopy (dfftp%npp (me_pool+1) * dfftp%nnp, ddvtosym (npp0, is, iper), &
             1, dvtosym (1, is, iper), 1)
     ENDDO

  ENDDO

  DEALLOCATE (ddvtosym)

#endif

  RETURN

END SUBROUTINE psyme
