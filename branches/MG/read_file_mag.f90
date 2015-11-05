!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE read_file_mag()
  !----------------------------------------------------------------------------
  !
  ! Wrapper routine, for compatibility
  !
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE io_global,            ONLY : stdout
  USE buffers,              ONLY : open_buffer, close_buffer
  USE wvfct,                ONLY : nbnd, npwx
  USE noncollin_module,     ONLY : npol
  USE klist,                ONLY : nks
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_onecenter,        ONLY : paw_potential
  USE uspp,                 ONLY : becsum
  USE scf,                  ONLY : rho
  USE realus,               ONLY : qpointlist, betapointlist, &
                                   init_realspace_vars,real_space
  USE dfunct,               ONLY : newd
  USE ldaU,                 ONLY : lda_plus_u, U_projection
  USE pw_restart,           ONLY : pw_readfile
  !
  IMPLICIT NONE 
  INTEGER :: ierr
  LOGICAL :: exst
  !
  ! ... Read the contents of the xml data file
  !
  CALL read_xml_file ( )
  !
  ! ... Open unit iunwfc, for Kohn-Sham orbitals
  !
  nwordwfc = nbnd*npwx*npol
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, nks, exst )
  !
  ! ... Read orbitals, write them in 'distributed' form to iunwfc
  !
  CALL pw_readfile( 'wave', ierr )
  !
  ! ... Assorted initialization: pseudopotentials, PAW
  ! ... Not sure which ones (if any) should be done here
  !
  CALL init_us_1()
  IF (lda_plus_U .AND. (U_projection == 'pseudo')) CALL init_q_aeps()
  !
  IF (okpaw) THEN
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_PAW)
  ENDIF 
  !
  IF ( real_space ) THEN
    CALL betapointlist()
    CALL init_realspace_vars()
    WRITE(stdout,'(5X,"Real space initialisation completed")')
  ENDIF
  CALL newd()
  !
  CALL close_buffer  ( iunwfc, 'KEEP' )
  !
END SUBROUTINE read_file_mag
!
!----------------------------------------------------------------------------
SUBROUTINE read_xml_file()
  !----------------------------------------------------------------------------
  !
  ! ... This routine allocates space for all quantities already computed
  ! ... in the pwscf program and reads them from the data file.
  ! ... All quantities that are initialized in subroutine "setup" when
  ! ... starting from scratch should be initialized here when restarting
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, nsp, ityp, tau, if_pos, extfor
  USE basis,                ONLY : natomwfc
  USE cell_base,            ONLY : tpiba2, alat,omega, at, bg, ibrav
  USE force_mod,            ONLY : force
  USE klist,                ONLY : nkstot, nks, xk, wk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, nbndx, et, wg, ecutwfc
  USE symm_base,            ONLY : irt, d1, d2, d3, checkallsym
  USE ktetra,               ONLY : tetra, ntetra 
  USE extfield,             ONLY : forcefield, tefield
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE grid_subroutines,     ONLY : realspace_grids_init
  USE recvec_subs,          ONLY : ggen
  USE gvect,                ONLY : gg, ngm, g, gcutm, &
                                   eigts1, eigts2, eigts3, nl, gstart
  USE fft_base,             ONLY : dfftp, dffts
  USE gvecs,                ONLY : ngms, nls, gcutms 
  USE spin_orb,             ONLY : lspinorb, domag
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE wavefunctions_module, ONLY : psic
  USE vlocal,               ONLY : strf
  USE io_files,             ONLY : tmp_dir, prefix, iunpun, nwordwfc, iunwfc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda, nspin_mag, nspin_gga
  USE pw_restart,           ONLY : pw_readfile
  USE read_pseudo_mod,      ONLY : readpp
  USE xml_io_base,          ONLY : pp_check_file
  USE uspp,                 ONLY : becsum
  USE uspp_param,           ONLY : upf
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_init,             ONLY : paw_init_onecenter, allocate_paw_internals
  USE ldaU,                 ONLY : lda_plus_u, eth, oatwfc
  USE control_flags,        ONLY : gamma_only
  USE funct,                ONLY : get_inlc, get_dft_name
  USE kernel_table,         ONLY : initialize_kernel_table
  USE esm,                  ONLY : do_comp_esm, esm_ggen_2d
  !
  INTEGER  :: i, is, ik, ibnd, nb, nt, ios, isym, ierr, inlc
  REAL(DP) :: rdum(1,1), ehart, etxc, vtxc, etotefield, charge
  REAL(DP) :: sr(3,3,48)
  CHARACTER(LEN=20) dft_name
  !
  ! ... first we get the version of the qexml file
  !     if not already read
  !
  CALL pw_readfile( 'header', ierr )
  CALL errore( 'read_xml_file ', 'unable to determine qexml version', ABS(ierr) )
  !
  ! ... then we check if the file can be used for post-processing
  !
  IF ( .NOT. pp_check_file() ) CALL infomsg( 'read_xml_file', &
               & 'file ' // TRIM( tmp_dir ) // TRIM( prefix ) &
               & // '.save not guaranteed to be safe for post-processing' )
  !
  ! ... a reset of the internal flags is necessary because some codes call
  ! ... read_xml_file() more than once
  !
  CALL pw_readfile( 'reset', ierr )
  !
  ! ... here we read the variables that dimension the system
  ! ... in parallel execution, only root proc reads the file
  ! ... and then broadcasts the values to all other procs
  !
  CALL pw_readfile( 'dim',   ierr )
  CALL errore( 'read_xml_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )
  !
  ! ... allocate space for atomic positions, symmetries, forces, tetrahedra
  !
  IF ( nat < 0 ) &
     CALL errore( 'read_xml_file', 'wrong number of atoms', 1 )
  !
  ! ... allocation
  !
  ALLOCATE( ityp( nat ) )
  ALLOCATE( tau(    3, nat ) )
  ALLOCATE( if_pos( 3, nat ) )
  ALLOCATE( force(  3, nat ) )
  ALLOCATE( extfor(  3, nat ) )
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat ) )
  !
  ALLOCATE( irt( 48, nat ) )
  ALLOCATE( tetra( 4, MAX( ntetra, 1 ) ) )
  !
  CALL set_dimensions()
  CALL realspace_grids_init ( dfftp, dffts, at, bg, gcutm, gcutms )
  !
  ! ... check whether LSDA
  !
  IF ( lsda ) THEN
     !
     nspin = 2
     npol  = 1
     !
  ELSE IF ( noncolin ) THEN
     !
     nspin        = 4
     npol         = 2
     current_spin = 1
     !
  ELSE
     !
     nspin        = 1
     npol         = 1
     current_spin = 1
     !
  END IF
  !
  if (cell_factor == 0.d0) cell_factor = 1.D0
  !
  ! ... allocate memory for eigenvalues and weights (read from file)
  !
  nbndx = nbnd
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ) )
  !
  ! ... here we read all the variables defining the system
  !
  CALL pw_readfile( 'nowave', ierr )
  !
  ! ... distribute across pools k-points and related variables.
  ! ... nks is defined by the following routine as the number 
  ! ... of k-points in the current pool
  !
  !CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks ) ! KC
  nks=nkstot    !KC
  !
  !CALL poolscatter( nbnd, nkstot, et, nks, et )
  !CALL poolscatter( nbnd, nkstot, wg, nks, wg )
  !
  ! ... check on symmetry
  !
  IF (nat > 0) CALL checkallsym( nat, tau, ityp, dfftp%nr1, dfftp%nr2, dfftp%nr3 )
  !
  !  Set the different spin indices
  !
  nspin_mag  = nspin
  nspin_lsda = nspin
  nspin_gga  = nspin
  IF (nspin==4) THEN
     nspin_lsda=1
     IF (domag) THEN
        nspin_gga=2
     ELSE
        nspin_gga=1
        nspin_mag=1
     ENDIF
  ENDIF
  !
  ! ... read pseudopotentials
  !
  CALL pw_readfile( 'pseudo', ierr )
  dft_name = get_dft_name () ! already set, should not be set again
  CALL readpp ( dft_name )
  !
  ! ... read the vdw kernel table if needed
  !
  inlc = get_inlc()
  if (inlc == 1 .or. inlc ==2 ) then
      call initialize_kernel_table()
  endif
  !
  okpaw = ANY ( upf(1:nsp)%tpawp )
  !
  IF ( .NOT. lspinorb ) CALL average_pp ( nsp )
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL pre_init()
  CALL allocate_fft()
  CALL ggen ( gamma_only, at, bg ) 
  IF (do_comp_esm) CALL esm_ggen_2d ()
  CALL gshells ( lmovecell ) 
  !
  ! ... allocate the potential and wavefunctions
  !
  CALL allocate_locpot()
  CALL allocate_nlpot()
  IF (okpaw) THEN
     CALL allocate_paw_internals()
     CALL paw_init_onecenter()
     CALL d_matrix(d1,d2,d3)
  ENDIF
  !
  IF ( lda_plus_u ) THEN
     ALLOCATE ( oatwfc(nat) )
     CALL offset_atom_wfc ( nat, oatwfc )
  ENDIF
  !
  CALL allocate_wfc()
  !
  ! ... read the charge density
  !
  CALL pw_readfile( 'rho', ierr )
  !
  ! ... re-calculate the local part of the pseudopotential vltot
  ! ... and the core correction charge (if any) - This is done here
  ! ... for compatibility with the previous version of read_file
  !
  CALL init_vloc()
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, dfftp%nr1, dfftp%nr2, &
                   dfftp%nr3, strf, eigts1, eigts2, eigts3 )
  CALL setlocal()
  CALL set_rhoc()
  !
  ! ... bring rho to G-space
  !
  DO is = 1, nspin
     !
     psic(:) = rho%of_r(:,is)
     CALL fwfft ('Dense', psic, dfftp)
     rho%of_g(:,is) = psic(nl(:))
     !
  END DO
  !
  ! ... read info needed for hybrid functionals
  !
  CALL pw_readfile('exx', ierr)
  !
  ! ... recalculate the potential
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_dimensions()
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : pi
      USE cell_base, ONLY : alat, tpiba, tpiba2
      USE gvect,     ONLY : ecutrho, gcutm
      USE wvfct,     ONLY : ecutwfc
      USE gvecs,     ONLY : gcutms, dual, doublegrid
      !
      !
      ! ... Set the units in real and reciprocal space
      !
      tpiba  = 2.D0 * pi / alat
      tpiba2 = tpiba**2
      !
      ! ... Compute the cut-off of the G vectors
      !
      gcutm = dual * ecutwfc / tpiba2
      ecutrho=dual * ecutwfc
      !
      doublegrid = ( dual > 4.D0 )
      IF ( doublegrid ) THEN
         gcutms = 4.D0 * ecutwfc / tpiba2
      ELSE
         gcutms = gcutm
      END IF
      !
    END SUBROUTINE set_dimensions
    !
END SUBROUTINE read_xml_file
