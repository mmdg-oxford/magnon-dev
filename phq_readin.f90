!
! Copyright (C) 2001-2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE phq_readin()
  !----------------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program phononq.
  !    from standard input (unit 5).
  !    A second routine readfile reads the variables saved on a file
  !    by the self-consistent program.
  !
  !
  USE kinds,         ONLY : DP
  USE parameters,    ONLY : nsx, npk
  USE constants,     ONLY : RYTOEV
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_bcast
  USE ions_base,     ONLY : amass, atm
  USE input_parameters, ONLY : max_seconds, nk1, nk2, nk3, k1, k2, k3
  USE start_k,       ONLY : reset_grid
  USE klist,         ONLY : xk, nks,wk, nkstot, lgauss, two_fermi_energies, lgauss
  USE ktetra,        ONLY : ltetra
  USE control_flags, ONLY : gamma_only, tqr, restart, lkpoint_dir
  USE uspp,          ONLY : okvan
  USE fixed_occ,     ONLY : tfixed_occ
  USE lsda_mod,      ONLY : lsda, nspin
  USE spin_orb,      ONLY : domag
  USE cellmd,        ONLY : lmovecell
  USE run_info, ONLY : title
  USE control_ph,    ONLY : maxter, alpha_mix, lgamma, lgamma_gamma, epsil, &
                            zue, zeu, xmldyn, newgrid,                      &
                            trans, reduce_io, tr2_ph, niter_ph,       &
                            nmix_ph, ldisp, recover, lrpa, lnoloc, start_irr, &
                            last_irr, start_q, last_q, current_iq, tmp_dir_ph, &
                            ext_recover, ext_restart, u_from_file, ldiag, &
                            search_sym, lqdir, electron_phonon, do_elec, dbext, qplot, dvext, &
                            transverse, symoff, man_kpoints, qpoints
  USE save_ph,       ONLY : tmp_dir_save
  USE gamma_gamma,   ONLY : asr
  USE qpoint,        ONLY : nksq, xq
  USE partial,       ONLY : atomo, nat_todo, nat_todo_input
  USE output,        ONLY : fildyn, fildvscf, fildrho
  USE disp,          ONLY : nq1, nq2, nq3, num_k_pts, xk_kpoints, kpoints,nqs
  USE io_files,      ONLY : tmp_dir, prefix
  USE noncollin_module, ONLY : i_cons, noncolin
  USE ldaU,          ONLY : lda_plus_u
  USE control_flags, ONLY : iverbosity, modenum, twfcollect
  USE io_global,     ONLY : ionode, stdout
  USE mp_global,     ONLY : nproc_pool, nproc_pool_file, &
                            nimage, my_image_id,    &
                            nproc_image_file, nproc_image, npool, &
                            get_ntask_groups, ntask_groups_file,  &
                            nbgrp
  USE paw_variables, ONLY : okpaw
  USE freq_ph,       ONLY : fpol, fiu, nfs, nfsmax
  USE ph_restart,    ONLY : ph_readfile
  USE xml_io_base,   ONLY : create_directory
  USE dfile_star,    ONLY : drho_star, dvscf_star
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios, ipol, iter, na, it, ierr
    ! integer variable for I/O control
    ! counter on polarizations
    ! counter on iterations
    ! counter on atoms
    ! counter on types
  REAL(DP) :: amass_input(nsx)
    ! save masses read from input here
  CHARACTER (LEN=256) :: outdir
  !
  CHARACTER(LEN=80)          :: card
  CHARACTER(LEN=1), EXTERNAL :: capital
  CHARACTER(LEN=6) :: int_to_char
  INTEGER                    :: i
  LOGICAL                    :: nogg
  INTEGER, EXTERNAL  :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight
  REAL(DP)           :: ar, ai
  LOGICAL, EXTERNAL  :: imatches
  LOGICAL, EXTERNAL  :: has_xml
  !
  NAMELIST / INPUTPH / tr2_ph, amass, alpha_mix, niter_ph, nmix_ph,  &
                       nat_todo, iverbosity, outdir, epsil,  &
                       trans,  zue, zeu, max_seconds, reduce_io, &
                       modenum, prefix, fildyn, fildvscf, fildrho, &
                       ldisp, nq1, nq2, nq3, recover,  &
                       fpol, asr, lrpa, lnoloc, start_irr, last_irr, &
                       start_q, last_q, nogg, ldiag, search_sym, lqdir, &
                       nk1, nk2, nk3, k1, k2, k3, &
                       drho_star, dvscf_star, qpoints, dbext, do_elec, &
                       dvext, transverse, symoff, man_kpoints

  ! tr2_ph       : convergence threshold
  ! amass        : atomic masses
  ! alpha_mix    : the mixing parameter
  ! niter_ph     : maximum number of iterations
  ! nmix_ph      : number of previous iterations used in mixing
  ! nat_todo     : number of atom to be displaced
  ! iverbosity   : verbosity control
  ! outdir       : directory where input, output, temporary files reside
  ! epsil        : if true calculate dielectric constant
  ! trans        : if true calculate phonon
  ! electron-phonon : select the kind of electron-phonon calculation
  ! zue          : if .true. calculate effective charges ( d force / dE )
  ! zeu          : if .true. calculate effective charges ( d P / du )
  ! lraman       : if true calculate raman tensor
  ! max_seconds  : maximum cputime for this run
  ! reduce_io    : reduce I/O to the strict minimum
  ! modenum      : single mode calculation
  ! prefix       : the prefix of files produced by pwscf
  ! fildyn       : output file for the dynamical matrix
  ! fildvscf     : output file containing deltavsc
  ! fildrho      : output file containing deltarho
  ! fildrho_dir  : directory where fildrho files will be stored (default: outdir or ESPRESSO_FILDRHO_DIR variable)
  ! eth_rps      : threshold for calculation of  Pc R |psi> (Raman)
  ! eth_ns       : threshold for non-scf wavefunction calculation (Raman)
  ! dek          : delta_xk used for wavefunctions derivation (Raman)
  ! recover      : recover=.true. to restart from an interrupted run
  ! asr          : in the gamma_gamma case apply acoustic sum rule
  ! start_q      : in q list does the q points from start_q to last_q
  ! last_q       :
  ! start_irr    : does the irred. representation from start_irr to last_irr
  ! last_irr     :
  ! nogg         : if .true. lgamma_gamma tricks are not used
  ! ldiag        : if .true. force diagonalization of the dyn mat
  ! lqdir        : if .true. each q writes in its own directory
  ! search_sym   : if .true. analyze symmetry if possible
  ! nk1,nk2,nk3,
  ! ik1, ik2, ik3: when specified in input it uses for the phonon run
  !                a different mesh than that used for the charge density.
  !
  ! dvscf_star%open : if .true. write in dvscf_star%dir the dvscf_q' for all q' in the
  !                   star of q with suffix dvscf_star%ext. The dvscf_q' is written in the basis dvscf_star%basis;
  !                   if dvscf_star%pat is .true. also save a pattern file.
  ! dvscf_star%dir, dvscf_star%ext, dvscf_star%basis : see dvscf_star%open
  ! drho_star%open  : like dvscf_star%open but for drho_q
  ! drho_star%dir, drho_star%ext, drho_star%basis : see drho_star%open
  
  !KC:
  ! dbext : for transverse perturbation with ground state magnetization along z
  ! direction, set dbext(1)=(1, 0), dbext(2)=(0, 1), dbext(3)=(0, 0)
  IF (ionode) THEN
  !
  ! ... Input from file ?
  !
     CALL input_from_file ( )
  !
  ! ... Read the first line of the input file
  !
     READ( 5, '(A)', IOSTAT = ios ) title
  !
  ENDIF
  !
  CALL mp_bcast(ios, ionode_id )
  CALL errore( 'phq_readin', 'reading title ', ABS( ios ) )
  CALL mp_bcast(title, ionode_id )
  !
  ! Rewind the input if the title is actually the beginning of inputph namelist
  !
  IF( imatches("&inputph", title) ) THEN
    WRITE(*, '(6x,a)') "Title line not specified: using 'default'."
    title='default'
    IF (ionode) REWIND(5, iostat=ios)
    CALL mp_bcast(ios, ionode_id )
    CALL errore('phq_readin', 'Title line missing from input.', abs(ios))
  ENDIF
  !
  ! ... set default values for variables in namelist
  !
  tr2_ph       = 1.D-4
  amass(:)     = 0.D0
  alpha_mix(:) = 0.D0
  alpha_mix(1) = 0.7D0
  niter_ph     = maxter
  nmix_ph      = 4
  nat_todo     = 0
  modenum      = 0
  iverbosity   = 0
  trans        = .TRUE.
  lrpa         = .FALSE.
  lnoloc       = .FALSE.
  epsil        = .FALSE.
  zeu          = .TRUE.
  zue          = .FALSE.
  fpol         = .FALSE.
  electron_phonon=' '
  max_seconds  =  1.E+7_DP
  reduce_io    = .FALSE.
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  prefix       = 'pwscf'
  fildyn       = 'matdyn'
  fildrho      = ' '
  fildvscf     = ' '
  ldisp        = .FALSE.
  nq1          = 0
  nq2          = 0
  nq3          = 0
  nogg         = .FALSE.
  recover      = .FALSE.
  asr          = .FALSE.
  start_irr    = 0
  last_irr     = -1000
  start_q      = 1
  last_q       =-1000
  ldiag        =.FALSE.
  lqdir        =.FALSE.
  search_sym   =.true.
  nk1       = 0
  nk2       = 0
  nk3       = 0
  k1       = 0
  k2       = 0
  k3       = 0
  !dbext(1) = (1.0, 0.0)
  !dbext(2) = (1.0, 0.0)
  !dbext(3) = (0.0, 0.0)
  dbext(1) = (0.d0, 0.d0)
  dbext(2) = (0.d0, 0.d0)
  dbext(3) = (0.d0, 0.d0)
  dvext = (0.d0, 0.d0)
!
  do_elec = .FALSE.
  transverse = .False.
  symoff = .false.
  qplot = .false.
 ! okvan=.TRUE.
  qpoints = .TRUE.
  man_kpoints = .FALSE.
  !
  drho_star%open = .FALSE.
  drho_star%basis = 'modes'
  drho_star%pat  = .TRUE.
  drho_star%ext = 'drho'
  CALL get_env( 'ESPRESSO_FILDRHO_DIR', drho_star%dir)
  IF ( TRIM( drho_star%dir ) == ' ' ) &
      drho_star%dir = TRIM(outdir)//"/Rotated_DRHO/"
  !
  dvscf_star%open = .FALSE.
  dvscf_star%basis = 'modes'
  dvscf_star%pat  = .FALSE.
  dvscf_star%ext = 'dvscf'
  CALL get_env( 'ESPRESSO_FILDVSCF_DIR', dvscf_star%dir)
  IF ( TRIM( dvscf_star%dir ) == ' ' ) &
      dvscf_star%dir = TRIM(outdir)//"/Rotated_DVSCF/"
  !
  ! ...  reading the namelist inputph
  !
  IF (ionode) READ( 5, INPUTPH, IOSTAT = ios )
  !
  CALL mp_bcast(ios, ionode_id)
  CALL errore( 'phq_readin', 'reading inputph namelist', ABS( ios ) )
  !
  IF (ionode) tmp_dir = trimcheck (outdir)

!  drho_star%dir=trimcheck(drho_star%dir)
!  dvscf_star%dir=trimcheck(dvscf_star%dir)
  ! filename for the star must always be automatically generated:
!  IF(drho_star%ext(1:5)/='auto:')  drho_star%ext  = 'auto:'//drho_star%ext
!  IF(dvscf_star%ext(1:5)/='auto:') dvscf_star%ext = 'auto:'//dvscf_star%ext

  CALL bcast_ph_input ( )
  CALL mp_bcast(nogg, ionode_id )
  !
  ! ... Check all namelist variables
  !
  IF (tr2_ph <= 0.D0) CALL errore (' phq_readin', ' Wrong tr2_ph ', 1)

  DO iter = 1, maxter
     IF (alpha_mix (iter) .LT.0.D0.OR.alpha_mix (iter) .GT.1.D0) CALL &
          errore ('phq_readin', ' Wrong alpha_mix ', iter)
  ENDDO
  IF (niter_ph.LT.1.OR.niter_ph.GT.maxter) CALL errore ('phq_readin', &
       ' Wrong niter_ph ', 1)
  IF (nmix_ph.LT.1.OR.nmix_ph.GT.5) CALL errore ('phq_readin', ' Wrong &
       &nmix_ph ', 1)
  IF (iverbosity.NE.0.AND.iverbosity.NE.1) CALL errore ('phq_readin', &
       &' Wrong  iverbosity ', 1)
!  IF (fildyn.EQ.' ') CALL errore ('phq_readin', ' Wrong fildyn ', 1)
  IF (max_seconds.LT.0.1D0) CALL errore ('phq_readin', ' Wrong max_seconds', 1)

!  IF (modenum < 0) CALL errore ('phq_readin', ' Wrong modenum ', 1)
!  IF (modenum /= 0) search_sym=.FALSE.
  
!  trans = trans .OR. ldisp
   
  !
  !  We can calculate  dielectric, raman or elop tensors and no Born effective
  !  charges dF/dE, but we cannot calculate Born effective charges dF/dE
  !  without epsil.
  !
  ! IF (zeu) zeu = epsil
  !
  !    reads the q point (just if ldisp = .false.)
  !
  !IF (ionode) THEN
  !   IF (.NOT. ldisp) &
  !      READ (5, *, iostat = ios) (xq (ipol), ipol = 1, 3)
  !END IF

 IF (qpoints) then
     num_k_pts = 0
     IF (ionode) THEN
        READ (5, *, iostat = ios) card
 !       print*, card
        IF ( TRIM(card)=='QPOINTS'.OR. &
             TRIM(card)=='qpoints'.OR. &
             TRIM(card)=='Qpoints') THEN
           READ (5, *, iostat = ios) num_k_pts
        ENDIF
        print*, card
     ENDIF
     
     CALL mp_bcast(ios, ionode_id )
     CALL errore ('phq_readin', 'reading number of qpoints', ABS(ios) )
     CALL mp_bcast(num_k_pts, ionode_id )
     nqs=num_k_pts
     if (num_k_pts > 50) call errore('phq_readin','Too many qpoints',1) 
     if (num_k_pts < 1) call errore('phq_readin','Too few qpoints',1) 
     IF (ionode) THEN
        !IF ( TRIM(card)=='K_POINTS'.OR. &
        !     TRIM(card)=='k_points'.OR. &
        !     TRIM(card)=='K_points') THEN
           DO i = 1, num_k_pts
             !should be in units of 2pi/a0 cartesian co-ordinates
              READ (5, *, iostat = ios) xk_kpoints(1,i), xk_kpoints(2,i), xk_kpoints(3,i)
           END DO
        !END IF
     END IF
     CALL mp_bcast(ios, ionode_id)
     CALL errore ('phq_readin', 'reading QPOINTS card', ABS(ios) )
     CALL mp_bcast(xk_kpoints, ionode_id )
 ELSE
     num_k_pts = 1
 ENDIF

! KC 
!  CALL mp_bcast(ios, ionode_id)
!  CALL errore ('phq_readin', 'reading xq', ABS (ios) )
!  CALL mp_bcast(xq, ionode_id )
!  IF (.NOT.ldisp) THEN
!     lgamma = xq (1) .EQ.0.D0.AND.xq (2) .EQ.0.D0.AND.xq (3) .EQ.0.D0
!     IF ( (epsil.OR.zue) .AND..NOT.lgamma) CALL errore ('phq_readin', &
!          'gamma is needed for elec.field', 1)
!  ENDIF

!  IF (zue.AND..NOT.trans) CALL errore ('phq_readin', 'trans must be &
!       &.t. for Zue calc.', 1)

  !IF (trans.AND.(lrpa.OR.lnoloc)) CALL errore('phq_readin', &
  !                  'only dielectric constant with lrpa or lnoloc',1)
  !
  ! reads the frequencies ( just if fpol = .true. )
  !
  IF ( fpol ) THEN
     nfs=0
     IF (ionode) THEN
        READ (5, *, iostat = ios) card
        IF ( TRIM(card)=='FREQUENCIES'.OR. &
             TRIM(card)=='frequencies'.OR. &
             TRIM(card)=='Frequencies') THEN
           READ (5, *, iostat = ios) nfs
        ENDIF
     ENDIF
     CALL mp_bcast(ios, ionode_id )
     CALL errore ('phq_readin', 'reading number of FREQUENCIES', ABS(ios) )
     CALL mp_bcast(nfs, ionode_id )
     if (nfs > nfsmax) call errore('phq_readin','Too many frequencies',1)
     if (nfs < 1) call errore('phq_readin','Too few frequencies',1)
     IF (ionode) THEN
        IF ( TRIM(card) == 'FREQUENCIES' .OR. &
             TRIM(card) == 'frequencies' .OR. &
             TRIM(card) == 'Frequencies' ) THEN
           DO i = 1, nfs
              READ (5, *, iostat = ios) ar, ai
              fiu(i)=dcmplx(ar, ai)/RYTOEV/1000.d0
           END DO
        END IF
     END IF
     CALL mp_bcast(ios, ionode_id)
     CALL errore ('phq_readin', 'reading FREQUENCIES card', ABS(ios) )
     CALL mp_bcast(fiu, ionode_id )
  ELSE
     nfs=0
     fiu=0.0_DP
  END IF

!Using image to separate (+-)q and (+-)fiu

if(my_image_id/=0)then
xk_kpoints=-xk_kpoints
fiu=-CONJG(fiu)
dbext=CONJG(dbext)
dvext=CONJG(dvext)
!end if
!fiu=-fiu
!write(*,*)'my_image_id', my_image_id, xk_kpoints(1,1), fiu(1),fiu(2)
end if
  !
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  !   amass will also be read from file:
  !   save its content in auxiliary variables
  !
  ! amass_input(:)= amass(:)
  !
  tmp_dir_save=tmp_dir
  tmp_dir_ph= TRIM (tmp_dir) // '_ph' // TRIM(int_to_char(my_image_id)) //'/'
  CALL create_directory(tmp_dir_ph)

  ext_restart=.FALSE.
  ext_recover=.FALSE.

1001 CONTINUE


  CALL read_file_mag ()

  !
  ! init_start_grid returns .true. if a new k-point grid is set from values
  ! read from input (this happens if nk1*nk2*nk3, else it returns .false.,
  ! leaves the current values, as read in read_file, unchanged)
  !
  !IF(.not. man_kpoints)newgrid = reset_grid (nk1, nk2, nk3, k1, k2, k3)
   newgrid = .true.


  !KC: input for the new kpoints grid, which includes a coarse grid for ground state
  !scf calculation and a dense selected grid for response only.
  !The grid may be highly nonuniform
  IF(ionode .and. man_kpoints)OPEN(1111,file='KPOINTS',status='old')
  
  !KC: kpoints must be input in cartisian coordinates
  IF (man_kpoints) then
     nkstot = 0
     IF (ionode) THEN
        READ (1111, *, iostat = ios) card
 !       print*, card
        IF ( TRIM(card)=='KPOINTS'.OR. &
             TRIM(card)=='kpoints'.OR. &
             TRIM(card)=='Kpoints') THEN
           READ (1111, *, iostat = ios) nkstot
        ENDIF
        print*, card
     ENDIF

     CALL mp_bcast(ios, ionode_id )
     CALL errore ('phq_readin', 'reading number of kpoints', ABS(ios) )
     CALL mp_bcast(nkstot, ionode_id )
     if (num_k_pts > npk) call errore('phq_readin','Too many kpoints',1)
     if (num_k_pts < 1) call errore('phq_readin','Too few kpoints',1)
     IF (ionode) THEN
           DO i = 1, nkstot
             !should be in units of 2pi/a0 cartesian co-ordinates
              READ (1111, *, iostat = ios) xk(1,i), xk(2,i), xk(3,i), wk(i)
           END DO
     END IF
     CALL mp_bcast(ios, ionode_id)
     CALL errore ('phq_readin', 'reading KPOINTS file', ABS(ios) )
     CALL mp_bcast(xk, ionode_id)
     CALL mp_bcast(wk, ionode_id)
 
 ENDIF  
  
   
  !
  tmp_dir=tmp_dir_save
  !

  IF (gamma_only) CALL errore('phq_readin',&
     'cannot start from pw.x data file using Gamma-point tricks',1)

  IF (okpaw.and.noncolin.and.domag) CALL errore('phq_readin',&
     'The phonon code with paw and domag is not available yet',1)

  IF (lmovecell) CALL errore('phq_readin', &
      'The phonon code is not working after vc-relax',1)

  IF (nproc_image /= nproc_image_file .and. .not. twfcollect)  &
     CALL errore('phq_readin',&
     'pw.x run with a different number of processors. Use wf_collect=.true.',1)

  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
     CALL errore('phq_readin',&
     'pw.x run with a different number of pools. Use wf_collect=.true.',1)

  IF (get_ntask_groups() > 1) &
     CALL errore('phq_readin','task_groups not available in phonon',1)

  IF (nbgrp > 1) &
     CALL errore('phq_readin','band parallelization not available in phonon',1)


  IF (.NOT.ldisp) lqdir=.FALSE.

  IF (i_cons /= 0) &
     CALL errore('phq_readin',&
     'The phonon code with constrained magnetization is not yet available',1)

  IF (two_fermi_energies .AND. (ltetra .OR. lgauss)) &
     CALL errore('phq_readin',&
     'The phonon code with two fermi energies is not available for metals',1)

  IF (tqr) CALL errore('phq_readin',&
     'The phonon code with Q in real space not available',1)

  IF (start_irr < 0 ) CALL errore('phq_readin', 'wrong start_irr',1)
  !
  IF (start_q <= 0 ) CALL errore('phq_readin', 'wrong start_q',1)
  !
  !  the dynamical matrix is written in xml format if fildyn ends in
  !  .xml or in the noncollinear case.
  !
  
  !
  ! If a band structure calculation needs to be done do not open a file
  ! for k point
  !
  lkpoint_dir=.FALSE.
  restart = recover
  !
  !  set masses to values read from input, if available;
  !  leave values read from file otherwise
  !
  lgamma_gamma=.FALSE.
!  IF (.NOT.ldisp) THEN
!     !
!     IF (lgamma) THEN
!        nksq = nks
!     ELSE
!        nksq = nks / 2
!     ENDIF
!  ENDIF
  !
  !
  !
  CALL allocate_part ( nat )
  !
  IF ( nat_todo < 0 .OR. nat_todo > nat ) &
     CALL errore ('phq_readin', 'nat_todo is wrong', 1)
  IF (nat_todo.NE.0) THEN
     IF (ionode) &
     READ (5, *, iostat = ios) (atomo (na), na = 1, &
          nat_todo)
     CALL mp_bcast(ios, ionode_id )
     CALL errore ('phq_readin', 'reading atoms', ABS (ios) )
     CALL mp_bcast(atomo, ionode_id )
  ENDIF
  nat_todo_input=nat_todo

  IF (epsil.AND.lgauss) &
        CALL errore ('phq_readin', 'no elec. field with metals', 1)
  IF (modenum > 0) THEN
     IF ( ldisp ) &
          CALL errore('phq_readin','Dispersion calculation and &
          & single mode calculation not possibile !',1)
     nat_todo = 0
  ENDIF

  !
  IF (ldisp .AND. (nq1 .LE. 0 .OR. nq2 .LE. 0 .OR. nq3 .LE. 0)) &
       CALL errore('phq_readin','nq1, nq2, and nq3 must be greater than 0',1)
  !
  RETURN
  !
END SUBROUTINE phq_readin
