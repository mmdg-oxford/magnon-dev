!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE phescf()
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver for the calculation of the
  ! ... response to an electric field and related quantities.
  !
  USE io_global,       ONLY : stdout
  USE paw_variables,   ONLY : okpaw
  USE uspp,            ONLY : okvan
  USE uspp_param,      ONLY : nhm
  USE ions_base,       ONLY : nat
  USE noncollin_module,ONLY : noncolin, nspin_mag
  USE lsda_mod,        ONLY : nspin
  USE control_ph,      ONLY : convt, zeu, rec_code, rec_code_read, lnoloc, &
                              lrpa, where_rec, done_epsil, done_zeu, epsil
  USE output,          ONLY : fildrho
  USE ph_restart,      ONLY : ph_writefile
  USE phus,            ONLY : int3, int3_nc, int3_paw
  USE freq_ph
  USE ramanm,          ONLY : ramtns, lraman, elop, done_lraman, done_elop
  !
  IMPLICIT NONE
  !
  INTEGER :: iu
  !
  !
  IF ( rec_code_read >  1 ) THEN
     IF (done_epsil) call summarize_epsilon()
     IF (done_zeu) call summarize_zeu()
     IF (done_elop) call summarize_elopt()
     IF (done_lraman) call write_ramtns(6,ramtns)
     RETURN
  ENDIF
  !
  IF (okvan) THEN
     ALLOCATE (int3 ( nhm, nhm, 3, nat, nspin_mag))
     IF (okpaw) ALLOCATE (int3_paw ( nhm, nhm, 3, nat, nspin_mag))
     IF (noncolin) ALLOCATE(int3_nc( nhm, nhm, 3, nat, nspin))
  ENDIF
  !
  IF (fpol) THEN    ! calculate freq. dependent polarizability
     !
     WRITE( stdout, '(/,5X,"Frequency Dependent Polarizability Calculation",/)' )
     !
     iu = nfs
     !
     freq_loop : DO WHILE ( iu .gt. 0)
        !
        CALL solve_e_fpol( fiu(iu) )
        IF ( convt ) CALL polariz ( fiu(iu) )
        iu = iu - 1
        !
     END DO freq_loop
     !
     WRITE( stdout, '(/,5X,"End of Frequency Dependent Polarizability Calculation")' )
     !
  ENDIF
  !
  IF ((epsil.AND..NOT.done_epsil).OR.(zeu.AND..NOT.done_zeu).OR.  &
      (lraman.AND..NOT.done_lraman).OR.(elop.AND..NOT.done_elop)) THEN

     WRITE( stdout, '(/,5X,"Electric Fields Calculation")' )
     !

     CALL solve_e()
     !
     WRITE( stdout, '(/,5X,"End of electric fields calculation")' )
     !
     IF ( convt ) THEN
        !
        ! ... calculate the dielectric tensor epsilon
        !
        IF (.NOT. done_epsil) THEN
           CALL dielec()
        ELSE
           CALL summarize_epsilon()
        ENDIF
        !
        ! ... calculate the effective charges Z(E,Us) (E=scf,Us=bare)
        !
        IF (.NOT.(lrpa.OR.lnoloc).AND.(zeu.AND..NOT.done_zeu)) THEN
           CALL zstar_eu()
        ELSEIF (done_zeu) THEN
           CALL summarize_zeu()
        ENDIF
        !
        IF ( fildrho /= ' ' ) CALL punch_plot_e()
        !
     ELSE
        !
        CALL stop_ph( .FALSE. )
        !
     END IF
     !
     IF ( (lraman.AND..NOT.done_lraman) .OR. (elop.AND..NOT.done_elop) &
                  .AND..NOT.noncolin) CALL raman()
     !
     where_rec='after_diel'
     rec_code=2
     CALL ph_writefile('data',0)
  ELSE
     IF (done_epsil) call summarize_epsilon()
     IF (done_zeu) call summarize_zeu()
     IF (done_elop) call summarize_elopt()
     IF (done_lraman) call write_ramtns(6,ramtns)
  ENDIF
  !
  IF (okvan) THEN
     DEALLOCATE (int3)
     IF (okpaw) DEALLOCATE (int3_paw)
     IF (noncolin) DEALLOCATE(int3_nc)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE phescf
