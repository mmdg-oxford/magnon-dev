! This file is copied and modified from QUANTUM ESPRESSO
! Kun Cao, Henry Lambert, Feliciano Giustino
 
!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE transform_dbecsum_nc(dbecsum_nc,dbecsum,na)
!----------------------------------------------------------------------------
!
! This routine multiply dbecsum_nc by the identity and the Pauli
! matrices and saves it in dbecsum to use it in the calculation of
! the charge and magnetization.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ntyp => nsp, ityp
USE uspp_param,           ONLY : nh, nhm
USE lsda_mod,             ONLY : nspin
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : domag
!
IMPLICIT NONE

INTEGER :: na, modes
COMPLEX(DP) :: dbecsum_nc( nhm, nhm, nat , nspin)
COMPLEX(DP) :: dbecsum( nhm * (nhm + 1) /2 , nat , nspin_mag)
!
! ... local variables
!
INTEGER :: ih, jh, ijh, np

np=ityp(na)

   ijh=1
   DO ih = 1, nh(np)
      dbecsum(ijh,na,1)= dbecsum(ijh,na,1)+  &
               dbecsum_nc(ih,ih,na,1)+dbecsum_nc(ih,ih,na,4)
      IF (domag) THEN
         dbecsum(ijh,na,2)= dbecsum(ijh,na,2)+  &
                  dbecsum_nc(ih,ih,na,2)+ &
                            dbecsum_nc(ih,ih,na,3)
         dbecsum(ijh,na,3)= dbecsum(ijh,na,3)+ &
                  (0.d0,-1.d0)*(dbecsum_nc(ih,ih,na,2)- &
                            dbecsum_nc(ih,ih,na,3) )
         dbecsum(ijh,na,4)= dbecsum(ijh,na,4)+  &
                  dbecsum_nc(ih,ih,na,1)-dbecsum_nc(ih,ih,na,4)
      END IF
      ijh=ijh+1
      DO jh = ih+1, nh(np)
         dbecsum(ijh,na,1)= dbecsum(ijh,na,1) +                   &
                   dbecsum_nc(ih,jh,na,1)+dbecsum_nc(ih,jh,na,4)  &
                  +dbecsum_nc(jh,ih,na,1)+dbecsum_nc(jh,ih,na,4)
         IF (domag) THEN
            dbecsum(ijh,na,2)= dbecsum(ijh,na,2) +     &
                      dbecsum_nc(ih,jh,na,2)+                 &
                             dbecsum_nc(ih,jh,na,3)     &
                   +  dbecsum_nc(jh,ih,na,2)+           &
                             dbecsum_nc(jh,ih,na,3)
            dbecsum(ijh,na,3)= dbecsum(ijh,na,3) +        &
                      (0.d0,-1.d0)*(dbecsum_nc(ih,jh,na,2)-    &
                                    dbecsum_nc(ih,jh,na,3)     &
                   +                dbecsum_nc(jh,ih,na,2)-    &
                                    dbecsum_nc(jh,ih,na,3) )
            dbecsum(ijh,na,4)= dbecsum(ijh,na,4) +     &
                      dbecsum_nc(ih,jh,na,1)-dbecsum_nc(ih,jh,na,4)+&
                      dbecsum_nc(jh,ih,na,1)-dbecsum_nc(jh,ih,na,4)
         END IF
         ijh=ijh+1
      END DO
   END DO

RETURN
END SUBROUTINE transform_dbecsum_nc
