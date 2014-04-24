!
! Copyright (C) 2001-2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for the phonon program
!
MODULE modes
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the modes and the small group of q
  !
  SAVE
  !
  INTEGER :: irgq(48), nsymq, irotmq, nirr, nmodes
  ! selects the operations of the small group
  ! the number of symmetry of the small group
  ! selects the symmetry sending q <-> -q+G
  ! number of irreducible representations contained in the dynamical matrix
  ! number of modes
  ! number of crystal sym.ops. for q=0
  INTEGER, ALLOCATABLE, TARGET :: npert(:) !3 * nat )
  ! the number of perturbations per IR
  INTEGER :: npertx
  ! max number of perturbations per IR
  REAL (DP), ALLOCATABLE :: rtau(:,:,:) !3, 48, nat)
  ! coordinates of direct translations
  REAL (DP) :: gi(3,48), gimq(3)
  ! the possible G associated to each symmetry
  ! the G associated to the symmetry q<->-q+G
  COMPLEX (DP), POINTER :: &
       u(:,:),                     &!  3 * nat, 3 * nat),
       t(:,:,:,:),                 &! npertx, npertx, 48,3 * nat),
       tmq(:,:,:)                   ! npertx, npertx, 3 * nat)
  ! the transformation modes patterns
  ! the mode for deltarho
  ! the symmetry in the base of the pattern
  ! the symmetry q<->-q in the base of the pa
  LOGICAL :: &
       minus_q,  &    !  if .TRUE. there is the symmetry sending q<->-q
       invsymq        !  if .TRUE. the small group of q has inversion

  CHARACTER(15), ALLOCATABLE :: name_rap_mode(:) ! symmetry type of each mode
  INTEGER, ALLOCATABLE :: num_rap_mode(:)  ! number of the representation for
                                           ! each mode
  !
END MODULE modes
!
!
MODULE dynmat
  USE kinds, ONLY :  DP
  !
  ! ... The dynamical matrix
  !
  SAVE
  !
  COMPLEX (DP), ALLOCATABLE :: &
       dyn00(:,:),           &! 3 * nat, 3 * nat),
       dyn(:,:),             &! 3 * nat, 3 * nat)
       dyn_rec(:,:)           ! 3 * nat, 3 * nat)
  ! the initial dynamical matrix
  ! the dynamical matrix
  ! the contribution of each representation to the dynamical matrix
  REAL (DP), ALLOCATABLE :: &
       w2(:)                  ! 3 * nat)
  ! omega^2
  !
END MODULE dynmat
!
!
MODULE qpoint
  USE kinds, ONLY :  DP
  USE parameters, ONLY : npk
  !
  ! ... The q point
  !
  SAVE
  !
  INTEGER, POINTER :: igkq(:)     ! npwx)
  ! correspondence k+q+G <-> G
  INTEGER :: nksq, npwq
  ! the real number of k points
  ! the number of plane waves for q
  INTEGER, ALLOCATABLE :: ikks(:), ikqs(:)
  ! the index of k point in the list of k
  ! the index of k+q point in the list of k
  REAL (DP) :: xq(3)
  ! the coordinates of the q point
  COMPLEX (DP), ALLOCATABLE :: eigqts(:) ! nat)
  ! the phases associated to the q
  !
END MODULE qpoint
!
!
MODULE eqv
  USE kinds, ONLY :  DP
  !
  ! ... The wavefunctions at point k+q
  !
  SAVE
  !
  COMPLEX (DP), POINTER :: evq(:,:)
  !
  ! ... The variable describing the linear response problem
  !
  COMPLEX (DP), ALLOCATABLE :: dvpsi(:,:), dpsi(:,:), drhoscfs (:,:,:)
  ! the product of dV psi
  ! the change of the wavefunctions
  REAL (DP), ALLOCATABLE :: dmuxc(:,:,:)        ! nrxx, nspin, nspin),
  REAL (DP), ALLOCATABLE, TARGET :: vlocq(:,:)  ! ngm, ntyp)
  ! the derivative of the xc potential
  ! the local potential at q+G
  REAL (DP), ALLOCATABLE :: eprec(:,:) ! needed for preconditioning
  !
END MODULE eqv
!
!
MODULE efield_mod
  USE kinds, ONLY :  DP
  !
  ! ... the variables for the electric field perturbation
  !
  SAVE
  !
  REAL (DP) :: epsilon (3, 3)
  REAL (DP), ALLOCATABLE :: &
       zstareu(:,:,:),       &! 3, 3, nat),
       zstarue(:,:,:)         ! 3, nat, 3)
  ! the dielectric constant
  ! the effective charges Z(E,Us) (E=scf,Us=bare)
  ! the effective charges Z(Us,E) (Us=scf,E=bare)
  COMPLEX (DP), ALLOCATABLE :: &
       zstareu0(:,:),        &! 3, 3 * nat),
       zstarue0(:,:),        &! 3 * nat, 3)
       zstarue0_rec(:,:)      ! 3 * nat, 3)
  ! the effective charges
  !
END MODULE efield_mod
!
!
MODULE nlcc_ph
  USE kinds, ONLY :  DP
  !
  ! ... The variables needed for non-linear core correction
  !
  SAVE
  !
  COMPLEX (DP), ALLOCATABLE, TARGET :: drc(:,:) ! ngm, ntyp)
  ! contain the rhoc (without structure fac) for all atomic types
  LOGICAL :: nlcc_any
  ! .T. if any atom-type has nlcc
  !
END MODULE nlcc_ph
!
!
MODULE gc_ph
  USE kinds, ONLY :  DP
  !
  ! ... The variables needed for gradient corrected calculations
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       grho(:,:,:),              &! 3, nrxx, nspin),
       gmag(:,:,:),              &! 3, nrxx, nspin),
       vsgga(:),                 &! nrxx
       segni(:),                 &! nrxx
       dvxc_rr(:,:,:),           &! nrxx, nspin, nspin), &
       dvxc_sr(:,:,:),           &! nrxx, nspin, nspin),
       dvxc_ss(:,:,:),           &! nrxx, nspin, nspin), &
       dvxc_s(:,:,:)              ! nrxx, nspin, nspin)
  !
  ! in the noncollinear case gmag contains the gradient of the magnetization
  ! grho the gradient of rho+ and of rho-, the eigenvalues of the spin density
  ! vsgga= 0.5* (V_up-V_down) to be used in the calculation of the change
  ! of the exchange and correlation magnetic field.
  ! gradient of the unpert. density
  !
  ! derivatives of the E_xc functiona
  ! r=rho and s=|grad(rho)|
  !
END MODULE gc_ph
!
!
MODULE phus
  USE kinds, ONLY :  DP
  USE becmod, ONLY : bec_type
  !
  ! ... These are additional variables needed for the linear response
  ! ... program with the US pseudopotentials
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       alphasum(:,:,:,:),   &! nhm*(nhm+1)/2,3,nat,nspin)
                             ! used to compute modes
       dpqq(:,:,:,:)         ! (nhm, nhm, 3, ntyp)
  ! alphasum contains \sum_i <psi_i| d/du (|\beta_n><beta_m|) | psi_i> + (m-n)
  ! dipole moment of each Q
  COMPLEX (DP), ALLOCATABLE :: &
       int1(:,:,:,:,:),     &! nhm, nhm, 3, nat, nspin),&
       int2(:,:,:,:,:),     &! nhm, nhm, 3,nat, nat),&
       int3(:,:,:,:,:),     &! nhm, nhm, npert, nat, nspin),&
       int3_paw(:,:,:,:,:), &! nhm, nhm, npert, nat, nspin),&
       int4(:,:,:,:,:),     &! nhm*(nhm+1)/2, 3, 3, nat, nspin),&
       int5(:,:,:,:,:),     &! nhm*(nhm+1)/2, 3, 3, nat, nat),&
       int1_nc(:,:,:,:,:),     &! nhm, nhm, 3, nat, nspin),&
       int2_so(:,:,:,:,:,:),   &! nhm, nhm, 3, nat,nat,nspin),&
       int3_nc(:,:,:,:,:),     &! nhm, nhm, npert, nat, nspin),&
       int4_nc(:,:,:,:,:,:),   &! nhm, nhm, 3, 3, nat, nspin),&
       int5_so(:,:,:,:,:,:,:), &! nhm*(nhm+1)/2, 3, 3, nat, nat, nspin),&
!
!  These variables contains the five integrals defined in PRB 64, 35118 (2001)
!  int1 -> \int V_eff d/du (Q) d^3r
!  int2 -> \int d/du (V_loc) Q d^3r
!  int3 -> \int d\du (V_Hxc) Q d^3r
!  int4 -> \int V_eff d^2/dudu (Q) d^3r
!  int5 -> \int d/du (V_loc) d/du (Q) d^3r
!
!  int3_paw contains d/du (D^1-\tilde D^1)
!
!
       becsum_nc(:,:,:,:),     &! nhm*(nhm+1)/2,nat,npol,npol)
       becsumort(:,:,:,:),     &! nhm*(nhm+1)/2,nat,nspin,3*nat)
       alphasum_nc(:,:,:,:,:), &! nhm*(nhm+1)/2,3,nat,npol,npol)
       dpqq_so(:,:,:,:,:)       ! nhm, nhm, nspin, 3, ntyp
!
!  becsum contains \sum_i <\psi_i | \beta_n><\beta_m| \psi_i > + (m-n)
!  besumort contains alphasum+\sum_i <\psi_i | \beta_n><\beta_m| \delta \psi_i >
!  dpqq_so dipole moment of each Q multiplied by the fcoef factors
!
  type (bec_type),  ALLOCATABLE, TARGET :: &
       becp1(:)              ! (nksq); (nkbtot, nbnd)
  !
  ! becp1 contains < beta_n | \psi_i >
  !
  type (bec_type),  ALLOCATABLE, TARGET :: &
       alphap(:,:)           ! nkbtot, nbnd, 3, nksq)
  !
  ! alphap contains < d\du (\beta_n) | psi_i>
  !
END MODULE phus
!
!
MODULE partial
  USE kinds, ONLY :  DP
  !
  ! ... the variables needed for partial computation of dynamical matrix
  !
  SAVE
  !
  INTEGER, ALLOCATABLE :: &
       comp_irr(:),    &! (3*nat) : 1 if this irr.rep. has to be computed
       done_irr(:),    &! (3*nat) : 1 if this irr.rep. has been done
       atomo(:)         ! (nat) : list of the atoms that moves
  INTEGER :: nat_todo,    & ! number of atoms to compute
             nat_todo_input ! nat_todo given in input
  LOGICAL :: all_comp       ! if .TRUE. all representation have been computed
  !
END MODULE partial
!
MODULE gamma_gamma
  INTEGER, ALLOCATABLE :: &
           has_equivalent(:),  &  ! 0 if the atom has to be calculated
           with_symmetry(:),   &  ! calculated by symmetry
           n_equiv_atoms(:),   &  ! number of equivalent atoms
           equiv_atoms(:,:)       ! which atoms are equivalent

  INTEGER :: n_diff_sites,    &   ! Number of different sites
             nasr                 ! atom calculated with asr
                                  !
  LOGICAL :: asr                  ! if true apply the asr

END MODULE gamma_gamma
!
MODULE control_ph
  USE kinds, ONLY :  DP
  USE parameters, ONLY: npk
  !
  ! ... the variable controlling the phonon run
  !
  SAVE
  !
  INTEGER, PARAMETER :: maxter = 100 ! maximum number of iterations
  INTEGER :: niter_ph,      & ! maximum number of iterations (read from input)
             nmix_ph,       & ! mixing type
             nbnd_occ(npk), & ! occupated bands in metals
             start_irr,     & ! initial representation
             last_irr,      & ! last representation of this run
             current_iq,    & ! current q point
             start_q, last_q  ! initial q in the list, last_q in the list
  REAL(DP) :: tr2_ph  ! threshold for phonon calculation
  REAL(DP) :: alpha_mix(maxter), & ! the mixing parameter
              time_now,          & ! CPU time up to now
              alpha_pv             ! the alpha value for shifting the bands
  CHARACTER(LEN=10)  :: where_rec='no_recover'! where the ph run recovered
  CHARACTER(LEN=12) :: electron_phonon
  CHARACTER(LEN=256) :: flmixdpot, tmp_dir_ph, tmp_dir_phq
  INTEGER :: rec_code,    &! code for recover
             rec_code_read=-1000 ! code for recover. Not changed during the run
  LOGICAL :: lgamma,      &! if .TRUE. this is a q=0 computation
             lgamma_gamma,&! if .TRUE. this is a q=0 computation with k=0 only
             convt,       &! if .TRUE. the phonon has converged
             epsil,       &! if .TRUE. computes dielec. const and eff. charges
             done_epsil=.FALSE.,  &! .TRUE. when diel. constant is available
             trans,       &! if .TRUE. computes phonons
             zue,         &! if .TRUE. computes eff. charges as induced polarization
             done_zue=.FALSE., &! .TRUE. when the eff. charges are available
             zeu,         &! if .TRUE. computes eff. charges as induced forces
             done_zeu=.FALSE., &! .TRUE. when the eff. charges are available
             recover,     &! if .TRUE. the run restarts
             ext_restart, &! if .TRUE. there is a restart file
             ext_recover, &! if .TRUE. there is a recover file
             lrpa,        &! if .TRUE. calculates the RPA dielectric constant
             lnoloc,      &! if .TRUE. calculates the dielectric constant
                           ! neglecting local field effects
             search_sym=.TRUE.,  &! if .TRUE. search the mode symmetry
             lnscf,       &! if .TRUE. the run makes first a nscf calculation
             ldisp,       &! if .TRUE. the run calculates full phonon dispersion
             reduce_io,   &! if .TRUE. reduces needed I/O
             done_bands,  &! if .TRUE. the bands have been calculated
             bands_computed=.FALSE., & ! if .TRUE. the bands were computed
                                       ! in this run
             nogg,        &! if .TRUE. gamma_gamma tricks are disabled
             u_from_file=.FALSE.,  & ! if true the u are on file
             recover_read=.FALSE., & ! if true the recover data have been read
             ldiag=.FALSE.,        & ! if true force the diagonalization
             lqdir=.FALSE.,        & ! if true each q writes in its directory
             xmldyn=.FALSE.,   & ! if true the dynamical matrix is in xml form
             all_done, &      ! if .TRUE. all representations have been done
             newgrid=.FALSE.  ! if .TRUE. use new k-point grid nk1,nk2,nk3
  !
END MODULE control_ph
!
!
MODULE freq_ph
  !
  USE kinds,   ONLY : DP
  !
  SAVE
  !
  ! ... the variables for computing frequency dependent dielectric constant
  !
  LOGICAL :: fpol ! if .TRUE. dynamic dielectric constant is computed
  !
  INTEGER, PARAMETER :: nfsmax=50  ! # of maximum frequencies
  INTEGER :: nfs                   ! # of frequencies
  !
  REAL (KIND=DP) :: fiu(nfsmax)    ! values  of frequency
  !
END MODULE freq_ph
!
!
MODULE units_ph
  !
  ! ... the units of the files and the record lengths
  !
  SAVE
  !
  INTEGER :: &
       iuwfc,     & ! iunit with the wavefunctions
       lrwfc,     & ! the length of wavefunction record
       iuvkb,     & ! unit with vkb
       iubar,     & ! unit with the part DV_{bare}
       lrbar,     & ! length of the DV_{bare}
       iuebar,    & ! unit with the part DV_{bare} for the electric field
       lrebar,    & ! length of the DV_{bare} fro the electric field
       iudwf,     & ! unit with D psi
       iupsir,    & ! unit with evc in real space
       lrdwf,     & ! length of D psi record
       iudrhous, lrdrhous, &
       iudyn,     & ! the unit for the dynamical matrix
       iupdyn,    & ! the unit for the partial dynamical matrix
       iunrec,    & ! the unit with the recover data
       iudvscf,   & ! the unit where the delta Vscf is written
       iudrho,    & ! the unit where the delta rho is written
       lrdrho,    & ! the length of the deltarho files
       iucom,     & ! the unit of the bare commutator in US case
       lrcom,     & ! the length  of the bare commutator in US case
       iudvkb3, lrdvkb3, &
       iuint3paw, & ! the unit of the int3_paw coefficients
       lint3paw     ! the lenght of the int3_paw coefficients
  ! the unit with the products
  ! the length of the products

  logical, ALLOCATABLE :: this_dvkb3_is_on_file(:), &
                          this_pcxpsi_is_on_file(:,:)
  !
END MODULE units_ph
!
!
MODULE output
  !
  ! ... the name of the files
  !
  SAVE
  !
  CHARACTER (LEN=256) :: fildyn, fildvscf, fildrho
  ! output file for the dynamical matrix
  ! output file for deltavscf
  ! output file for deltarho
  !
END MODULE output
!
!
MODULE disp
  !
  USE kinds, ONLY: DP
  !
  SAVE
  !
  INTEGER, PARAMETER :: nqmax = 1000
  !
  INTEGER :: nq1, nq2, nq3  ! number of q-points in each direction
  INTEGER :: nqs            ! number of q points to be calculated
  REAL(DP), ALLOCATABLE :: x_q(:,:) ! coordinates of the q points
  INTEGER, ALLOCATABLE :: &
       done_iq(:),      &! if 1 this q point has been already calculated
       comp_iq(:),      &! if 1 this q point has to be calculated
       rep_iq(:),       &! number of irreducible representation per q point
       done_rep_iq(:,:),&! which representation have been already done in each q
       nsymq_iq(:),     &! dimension of the small group of q
       comp_irr_iq(:,:),&! for each q, comp_irr. Used for image parallelization
       npert_iq(:,:)     ! for each q, the number of perturbation of each irr
  !
END MODULE disp
!
!
MODULE phcom
  USE modes
  USE dynmat
  USE qpoint
  USE eqv
  USE efield_mod
  USE nlcc_ph
  USE gc_ph
  USE phus
  USE partial
  USE control_ph
  USE freq_ph
  USE units_ph
  USE output
  USE gamma_gamma
  USE disp
END MODULE phcom
