!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE prepare_q(do_band, do_iq, setup_pw, iq, minq)
  !-----------------------------------------------------------------------
  !
  !  This routine prepares a few variables that are needed to control
  !  the GW run after the q point has been decided, but before
  !  doing the band calculation. 
  !  In particular if ldisp=true it sets:
  !  xq : the q point for the q calculation
  !  current_iq : the current q point
  !  do_iq : if .true. q point has to be calculated
  !  setup_pw : if .true. the pw_setup has to be run
  !  do_band : if .true. the bands need to be calculated before phonon
  !  fildyn : the name of the dynamical matrix
  !  lgamma : if this is a gamma point calculation
  !  epsil and zue : if epsil and zue need to be calculated
  !  In all cases it sets:
  
  USE control_flags,   ONLY : modenum
  USE io_global,       ONLY : stdout, ionode
  USE klist,           ONLY : lgauss
  USE qpoint,          ONLY : xq
  USE disp,            ONLY : x_q, done_iq, rep_iq, done_rep_iq, comp_iq,&
                              xk_kpoints
  USE control_ph,      ONLY : ldisp, lgamma, epsil, trans, zue, zeu, &
                              start_irr, last_irr, current_iq, &
                              done_bands, tmp_dir_phq, tmp_dir_ph
  USE freq_ph,         ONLY : fpol
  USE output,          ONLY : fildyn
  USE io_files,        ONLY : prefix
  USE ph_restart,      ONLY : ph_writefile
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  LOGICAL, INTENT(INOUT) :: minq
  LOGICAL, INTENT(OUT) :: do_band, do_iq, setup_pw
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  INTEGER :: irr
  !
  do_iq=.TRUE.
  !
  ! Case 1) This q point is not calculated because not requested in this run
  !

  IF ( comp_iq(iq)==0 ) THEN
     do_iq=.FALSE.
     RETURN
  ENDIF
  !
  !  Case 2) This q point is not calculated because it has too few 
  !          representation and the starting representation is larger 
  !          than the number of available representations
  !
  current_iq = iq
  tmp_dir_phq=tmp_dir_ph
  !
     ! ... set the name for the output file
     ! ... set the q point
!MINUS Q
!       xq(1:3)  = x_q(1:3,iq)
       xq(1:3)  = xk_kpoints(1:3,iq)
       if ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 ) xq(1) = 0.1
!      xq(:) = xk_kpoints(:, 1)
!Check if it is lgamma
       lgamma = (xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0)
        ! ... for q /= 0 no calculation of the dielectric tensor,
        ! ...            Born eff. charges, electro-optic, raman or
        ! ...            frequency dependent tensor
        !
       epsil = .FALSE.
       zue   = .FALSE.
       zeu   = .FALSE.
       fpol  = .FALSE.
  !
  !  Save the current status of the run: all the flags, the list of q,
  !  and the current q, the fact that we are before the bands
  !
  CALL ph_writefile('init',0)
  !
  ! ... In the case of q != 0, we make first a non selfconsistent run
  !
  setup_pw = (.NOT.lgamma.OR.modenum /= 0).AND..NOT. done_bands
  do_band=.FALSE.
  DO irr=start_irr, MIN(ABS(last_irr),rep_iq(iq))
     IF (done_rep_iq(irr,iq) /= 1) THEN
        do_band=.TRUE.
        EXIT
     ENDIF
  ENDDO
!
!  There are two special cases. When start_irr=0 and last_irr=0 we generate only
!  the displacement patterns, and do not calculate the bands. If this q
!  has been already calculated we only diagonalize the dynamical matrix
!
  WRITE( stdout, '(/,5X,"Calculation of q = ",3F12.7)') xq
  RETURN
END SUBROUTINE prepare_q
