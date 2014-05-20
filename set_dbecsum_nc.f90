!
! Copyright (C) 2007-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!
!----------------------------------------------------------------------------
SUBROUTINE set_dbecsum_nc(dbecsum_nc, dbecsum)
!----------------------------------------------------------------------------
USE kinds, ONLY : DP
USE ions_base, ONLY : nat, ntyp => nsp, ityp
USE uspp_param, only: upf, nhm
USE noncollin_module, ONLY : nspin_mag
USE lsda_mod, ONLY : nspin
IMPLICIT NONE
INTEGER :: npe
INTEGER :: np, na
COMPLEX(DP), INTENT(IN) :: dbecsum_nc( nhm, nhm, nat, nspin)
COMPLEX(DP), INTENT(OUT) :: dbecsum( nhm*(nhm+1)/2, nat, nspin_mag)

DO np = 1, ntyp
   IF ( upf(np)%tvanp ) THEN
      DO na = 1, nat
         IF (ityp(na)==np) THEN
            IF (upf(np)%has_so) THEN
               CALL transform_dbecsum_so(dbecsum_nc,dbecsum,na)
            ELSE
               CALL transform_dbecsum_nc(dbecsum_nc,dbecsum,na)
            END IF
         END IF
      END DO
   END IF
END DO

RETURN
END SUBROUTINE set_dbecsum_nc
