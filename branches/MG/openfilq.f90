!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE openfilq()
  !----------------------------------------------------------------------------
  !
  ! ... This subroutine opens all the files necessary for the phononq
  ! ... calculation.
  !
  USE kinds,           ONLY : DP
  USE control_flags,   ONLY : modenum
  USE units_ph,        ONLY : iuwfc, iudwf, iubar, iucom, iudvkb3, &
                              iudrhous, iuebar, iudrho, iudyn, iudvscf, &
                              lrwfc, lrdwf, lrbar, lrcom, lrdvkb3, &
                              lrdrhous, lrebar, lrdrho, lint3paw, iuint3paw
  USE io_files,        ONLY : tmp_dir, diropn, seqopn
  USE control_ph,      ONLY : epsil, zue, ext_recover, trans, lgamma, &
                              tmp_dir_phq, start_irr, last_irr, xmldyn, &
                              all_done
  USE save_ph,         ONLY : tmp_dir_save
  USE ions_base,       ONLY : nat
  USE cell_base,       ONLY : at
  USE qpoint,          ONLY : xq, nksq
  USE output,          ONLY : fildyn, fildvscf
  USE wvfct,           ONLY : nbnd, npwx
  USE fft_base,        ONLY : dfftp, dffts
  USE lsda_mod,        ONLY : nspin
  USE uspp,            ONLY : nkb, okvan
  USE uspp_param,      ONLY : nhm
  USE io_files,        ONLY : prefix, iunigk
  USE noncollin_module,ONLY : npol, nspin_mag
  USE paw_variables,   ONLY : okpaw
  USE control_flags,   ONLY : twfcollect
  USE mp_global,       ONLY : me_pool
  USE io_global,       ONLY : ionode,stdout
  USE input_parameters,ONLY : nk1, nk2, nk3
  USE dfile_star,      ONLY : dvscf_star
  USE dfile_autoname,  ONLY : dfile_name
  !
  IMPLICIT NONE
  !
  INTEGER :: ios
  ! integer variable for I/O control
  CHARACTER (len=256) :: filint, fildvscf_rot
  ! the name of the file
  LOGICAL :: exst
  ! logical variable to check file existe
  !
  REAL(DP) :: edum(1,1), wdum(1,1)
  INTEGER :: ndr, ierr, iq_dummy
  INTEGER, EXTERNAL :: find_free_unit
  !
  !
  IF (LEN_TRIM(prefix) == 0) CALL errore ('openfilq', 'wrong prefix', 1)
  !
  !     There are six direct access files to be opened in the tmp area
  !
  !     The file with the wavefunctions. In the lgamma case reads those
  !     written by pw.x. In the other cases those calculated by ph.x
  !
  tmp_dir=tmp_dir_phq
  IF (lgamma.AND.modenum==0.AND.nk1.eq.0.AND.nk2.eq.0.AND.nk3.eq.0) tmp_dir=tmp_dir_save

  iuwfc = 20
  lrwfc = 2 * nbnd * npwx * npol
  CALL diropn (iuwfc, 'wfc', lrwfc, exst)
  IF (.NOT.exst.and..not.all_done) THEN
     CALL errore ('openfilq', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF
  !
  ! From now on all files are written with the _ph prefix
  !
  tmp_dir=tmp_dir_phq
  !
  !    The file with deltaV_{bare} * psi
  !
  iubar = 21
  lrbar = 2 * nbnd * npwx * npol
  CALL diropn (iubar, 'bar', lrbar, exst)
  IF (ext_recover.AND..NOT.exst) &
     CALL errore ('openfilq','file '//trim(prefix)//'.bar not found', 1)
  !
  !    The file with the solution delta psi
  !
  iudwf = 22
  lrdwf = 2 * nbnd * npwx * npol
  CALL diropn (iudwf, 'dwf', lrdwf, exst)
  IF (ext_recover.AND..NOT.exst) &
     CALL errore ('openfilq','file '//trim(prefix)//'.dwf not found', 1)
  !
  !   open a file with the static change of the charge
  !
  IF (okvan) THEN
     iudrhous = 25
     lrdrhous = 2 * dfftp%nnr * nspin_mag
     CALL diropn (iudrhous, 'prd', lrdrhous, exst)
     IF (ext_recover.AND..NOT.exst) &
        CALL errore ('openfilq','file '//trim(prefix)//'.prd not found', 1)
  ENDIF
  !
  !  Optional file(s) containing Delta\rho (opened and written in solve_e
  !  and solve_linter). Used for third-order calculations.
  !
  iudrho = 23
  lrdrho = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin_mag
  !
  !   Here the sequential files
  !   The igk at a given k (and k+q if q!=0)
  !
  iunigk = 24
  IF (nksq > 1) CALL seqopn (iunigk, 'igk', 'unformatted', exst)
  !
  !   a formatted file which contains the dynamical matrix in cartesian
  !   coordinates is opened in the current directory
  !   ... by the first node only, other nodes write on unit 6 (i.e./dev/null
  !   exception: electron-phonon calculation from saved data
  !  (iudyn is read, not written, by all nodes)
  RETURN
  !
END SUBROUTINE openfilq
