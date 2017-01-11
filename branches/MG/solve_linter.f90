!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE solve_linter (drhoscf, iw)
  !-----------------------------------------------------------------------
  !
  !    Driver routine for the solution of the linear system which
  !    defines the change of the wavefunction due to a lattice distorsion
  !    It performs the following tasks:
  !     a) computes the bare potential term Delta V | psi >
  !        and an additional term in the case of US pseudopotentials
  !     b) adds to it the screening term Delta V_{SCF} | psi >
  !     c) applies P_c^+ (orthogonalization to valence states)
  !     d) calls cgsolve_all to solve the linear system
  !     e) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, iunigk, diropn
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions_module, ONLY : evc, evc0
  USE constants,            ONLY : degspin
  USE cell_base,            ONLY : at, tpiba2
  USE ener,                 ONLY : ef
  USE klist,                ONLY : lgauss, degauss, ngauss, xk, wk, nkstot
  USE gvect,                ONLY : g
  USE gvecs,                ONLY : doublegrid, nls
  USE fft_base,             ONLY : dfftp, dffts
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE spin_orb,             ONLY : domag
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf, nhm, nh
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE paw_variables,        ONLY : okpaw
  USE paw_onecenter,        ONLY : paw_dpotential
  USE paw_symmetry,         ONLY : paw_dusymmetrize, paw_dumqsymmetrize
  USE control_ph,           ONLY : rec_code, niter_ph, nmix_ph, tr2_ph, &
                                   alpha_pv, lgamma, lgamma_gamma, convt, &
                                   nbnd_occ, alpha_mix, ldisp, rec_code_read, &
                                   where_rec, flmixdpot, ext_recover, do_elec, &
                                   transverse, dbext, lrpa, dvext, man_kpoints, &
                                   symoff, reduce_io, thresh_CG
  USE nlcc_ph,              ONLY : nlcc_any
  USE units_ph,             ONLY : iudrho, lrdrho, iudwfp, iudwfm, lrdwf, iubar, lrbar, iudwf, &
                                   iuwfc, lrwfc, iunrec, iudvscf, &
                                   this_pcxpsi_is_on_file
  USE output,               ONLY : fildrho, fildvscf
  USE phus,                 ONLY : int3_paw, becsumort
  USE eqv,                  ONLY : dvpsi, dpsi, evq, eprec, dvpsi0, dpsi0  !, dpsip,dpsim
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE modes,                ONLY : npertx, npert, u, t, irotmq, tmq, &
                                   minus_q, nsymq, rtau
  USE recover_mod,          ONLY : read_rec, write_rec
  ! used to write fildrho:
  USE dfile_autoname,       ONLY : dfile_name
  USE save_ph,              ONLY : tmp_dir_save
  ! used oly to write the restart file
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm, inter_image_comm, my_image_id
  USE mp,                   ONLY : mp_sum, mp_barrier
  USE freq_ph,       ONLY : fpol, fiu, nfs, nfsmax
  USE fft_interfaces, ONLY: fwfft, invfft
  USE ktetra,         ONLY: ltetra
  !
  implicit none

  integer :: irr, npe, imode0, convt_check
  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes

  !complex(DP) :: drhoscf (dfftp%nnr, nspin_mag)
   complex(DP) :: drhoscf (dffts%nnr, nspin_mag)  !KC
  ! output: the change of the scf charge

  real(DP) , allocatable :: h_diag (:,:)
  ! h_diag: diagonal part of the Hamiltonian
  real(DP) :: thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP) :: dos_ef, weight, aux_avg (2)
  ! Misc variables for metals
  ! dos_ef: density of states at Ef
  real(DP), external :: w0gauss, wgauss
  ! functions computing the delta and theta function

!############# for test purpose ###############
!  complex(DP), allocatable::dpsim(:,:),dpsip(:,:)
!#####################################################################


  complex(DP), allocatable, target :: dvscfin(:,:)
  ! change of the scf potential
  complex(DP), pointer :: dvscfins (:,:)
  ! change of the scf potential (smooth part only)
  complex(DP), allocatable :: drhoscfh (:,:), dvscfout (:,:)
  ! change of rho / scf potential (output)
  ! change of scf potential (output)
  complex(DP), allocatable :: ldos (:,:), ldoss (:,:), mixin(:), mixout(:), &
       dbecsum (:,:,:), dbecsum_nc(:,:,:,:), aux1 (:,:)
  ! Misc work space
  ! ldos : local density of states af Ef
  ! ldoss: as above, without augmentation charges
  ! dbecsum: the derivative of becsum
  REAL(DP), allocatable :: becsum1(:,:,:)

  logical :: conv_root,  & ! true if linear system is converged
             exst,       & ! used to open the recover file
             lmetq0,     &   ! true if xq=(0,0,0) in a metal
             convtm

  integer :: kter,       & ! counter on iterations
             iter0,      & ! starting iteration
             ipert,      & ! counter on perturbations
             ibnd,       & ! counter on bands
             iter,       & ! counter on iterations
             lter,       & ! counter on iterations of linear system
             ltaver,     & ! average counter
             lintercall, & ! average number of calls to cgsolve_all
             ik, ikk, jk, & ! counter on k points
             ikq,        & ! counter on k+q points
             ig,         & ! counter on G vectors
             ndim,       &
             is,         & ! counter on spin polarizations
             nt,         & ! counter on types
             na,         & ! counter on atoms
             nrec,       & ! the record number for dvpsi and dpsi
             ios,        & ! integer variable for I/O control
             mode,       &   ! mode index
             im
  integer  :: iq_dummy
  real(DP) :: tcpu, get_clock ! timing variables
  character(len=256) :: filename
  complex(DP) :: cw
  integer::iw
  complex(DP), allocatable :: etc(:,:)
  real(kind=DP) :: DZNRM2

  external ch_psi_all, cg_psi, cch_psi_all, DZNRM2
  !
  IF (rec_code_read > 20 ) RETURN

  npe    = 1
  imode0 = 1
  irr    = 1
  ipert  = 1
  lter   = 0

  call start_clock ('solve_linter')

   allocate (dvscfin ( dfftp%nnr , nspin_mag))
  if (doublegrid) then
     allocate (dvscfins (dffts%nnr , nspin_mag))
  else
     dvscfins => dvscfin
  endif

!############## for test purpose ######################
!  allocate (dpsip(dffts%nnr, nbnd))
!  allocate (dpsim(dffts%nnr, nbnd))
!################################################################
  
  allocate (drhoscfh ( dfftp%nnr, nspin_mag))
  allocate (dvscfout ( dfftp%nnr, nspin_mag ))
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat , nspin_mag))
  allocate (etc(nbnd, nkstot))



  IF (okpaw) THEN
     allocate (mixin(dfftp%nnr*nspin_mag+(nhm*(nhm+1)*nat*nspin_mag)/2) )
     allocate (mixout(dfftp%nnr*nspin_mag+(nhm*(nhm+1)*nat*nspin_mag)/2) )
     mixin=(0.0_DP,0.0_DP)
  ENDIF


  IF (noncolin) allocate (dbecsum_nc (nhm, nhm, nat, nspin))
  allocate (aux1 ( dffts%nnr, npol))
  allocate (h_diag ( npwx*npol, nbnd))
  !
  iter0 = 0
  convt =.FALSE.
  where_rec='no_recover'

  IF (ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     filename = dfile_name(xq, at, fildrho, TRIM(tmp_dir_save)//prefix, generate=.true., index_q=iq_dummy)
     CALL diropn (iudrho, filename, lrdrho, exst)
  END IF

  IF (convt) GOTO 155

  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !

  lmetq0 = (lgauss .or. ltetra ).and.lgamma

  if (lmetq0) then
        allocate ( ldos ( dfftp%nnr  , nspin_mag) )
        allocate ( ldoss( dffts%nnr , nspin_mag) )
        write(stdout,*)'test local dos'
        call localdos ( ldos , ldoss , dos_ef )
  write(stdout, *)'dos_ef', dos_ef
  endif


!################################# test the wavefunctions #######################

!IF(reduce_io == .true.)THEN


!Do ik=1, 2*nksq
!  Do jk = ik+1, 2*nksq
!  IF((xk (1,ik)+xk (1,jk))**2+ (xk (2,ik)+xk (2,jk))**2+ (xk (3,ik)+xk(3,jk))**2 < 1e-8)THEN
!  WRITE( stdout,'(I8,3f10.6)' )ik, xk (1,ik), xk (2,ik), xk (3,ik)
!  WRITE( stdout,'(I8,3f10.6)' )jk, xk (1,jk), xk (2,jk), xk (3,jk)
!  cw = (0.d0, 0.d0)
!   Do ibnd = 1, nbnd
!   DZNRM2 (2*npwx, evc0(1, ibnd, ik)-evc0(1, ibnd, jk), 1)
!   DZNRM2 (2*npwx, evc0(1, ibnd, ik)+evc0(1, ibnd, jk), 1)
!   evc(:,ibnd)= evc0(1, ibnd, ik)-evc0(1, ibnd,jk)
!   evq(:,ibnd)= evc0(1, ibnd, ik)+evc0(1, ibnd,jk)
!   WRITE( stdout, '(2f12.8)')DZNRM2 (npol*npwx, evc(1, ibnd),1), DZNRM2 (npol*npwx,  evq(1, ibnd), 1)
!   END DO
!   WRITE( stdout, '(4f12.8)')DCONJG(evc0(nls(1), ibnd, ik))*evc0(nls(1), ibnd, ik),DCONJG(evc0(nls(1), ibnd, jk))*evc0(nls(1), ibnd, jk)
!   END DO
!  WRITE( stdout, '(6f12.8)')evc0(10, nbnd/2, jk), evc0(npwx/2, nbnd/2, jk), evc0(50, 5, jk)
!  WRITE( stdout, '(a36,2I8,2f12.8)')'npwx, nbnd, evq(npwx/2, nbnd/2)',npwx, nbnd, &
!  &  evc0(npwx/2, nbnd/2, ik)-evc0(npwx/2, nbnd/2, jk)
!  END IF
!  END DO
!END DO


!END IF
! look at them in the real space

!################################################################################
  !
  !
  ! In this case it has recovered after computing the contribution
  ! to the dynamical matrix. This is a new iteration that has to
  ! start from the beginning.
  !
  IF (iter0==-1000) iter0=0     !start iteration
  !
  !   The outside loop is over the iterations
  !
  do kter = 1, niter_ph
!   do kter = 1, 1
     iter = kter + iter0

     ltaver = 0

     lintercall = 0
     drhoscf(:,:) = (0.d0, 0.d0)
!mark     drhoscfh(:,:)  = (0.d0, 0.d0) 
     dbecsum(:,:,:) = (0.d0, 0.d0)
     IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
     !
     if (nksq.gt.1) rewind (unit = iunigk)

     do ik = 1, nksq

        if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
100        call errore ('solve_linter', 'reading igk', abs (ios) )
        endif


        if (lgamma)  npwq = npw

        ikk = ikks(ik)
        ikq = ikqs(ik)

!############################## for test purpose ##########################

!        WRITE(stdout, *) 'ikk, igk(10)', ikk, igk(10)

!#####################################################################     

        if (lsda) current_spin = isk (ikk)

        if (.not.lgamma.and.nksq.gt.1) then
           read (iunigk, err = 200, iostat = ios) npwq, igkq
200        call errore ('solve_linter', 'reading igkq', abs (ios) )

        endif

        call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
        !
        ! reads unperturbed wavefuctions psi(k) and psi(k+q)
        !
IF(reduce_io)THEN 
        if (nksq.gt.1) then
           if (lgamma) then
              evc(:,:)=evc0(:,:,ikk)
           else
              evc(:,:)=evc0(:,:,ikk)
              evq(:,:)=evc0(:,:,ikq)
           endif

        endif
ELSE
        if (nksq.gt.1) then
           if (lgamma) then
              call davcio (evc, lrwfc, iuwfc, ikk, - 1)
           else
              call davcio (evc, lrwfc, iuwfc, ikk, - 1)
              call davcio (evq, lrwfc, iuwfc, ikq, - 1)
           endif

        endif
END IF


        !
        ! compute the kinetic energy
        !
        do ig = 1, npwq
           g2kin (ig) = ( (xk (1,ikq) + g (1, igkq(ig)) ) **2 + &
                          (xk (2,ikq) + g (2, igkq(ig)) ) **2 + &
                          (xk (3,ikq) + g (3, igkq(ig)) ) **2 ) * tpiba2
        enddo

        h_diag=0.d0
        do ibnd = 1, nbnd_occ (ikk)
           do ig = 1, npwq
              h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
           enddo
           IF (noncolin) THEN
              do ig = 1, npwq
                 h_diag(ig+npwx,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
              enddo
           END IF
        enddo
        !
        ! diagonal elements of the unperturbed hamiltonian
        !
!           mode =1
           nrec = ik         ! record position
           !
           !  and now adds the contribution of the self consistent term
           !
           if (where_rec =='solve_lint'.or.iter>1) then
              !
              ! After the first iteration dvbare_q*psi_kpoint is read from file
              !
              IF(reduce_io)THEN
              dvpsi(:,:)=dvpsi0(:,:,ik)
              ELSE
              call davcio (dvpsi, lrbar, iubar, nrec, - 1)
              ENDIF
              !
              ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
              ! dvscf_q from previous iteration (mix_potential)
              !
              call start_clock ('vpsifft')
              do ibnd = 1, nbnd_occ (ikk)
                 call cft_wave (evc (1, ibnd), aux1, +1)
                 call apply_dpot(dffts%nnr, aux1, dvscfins(1,1), current_spin)
                 !call apply_dmag(dffts%nnr, aux1, dvscfins(1,1), current_spin)
                 call cft_wave (dvpsi (1, ibnd), aux1, -1)
              enddo
              call stop_clock ('vpsifft')
              !
              !  In the case of US pseudopotentials there is an additional
              !  selfconsist term which comes from the dependence of D on
              !  V_{eff} on the bare change of the potential
              !
              !Need to check this for ultrasoft              
              !HL THIS TERM PROBABLY NEEDS TO BE INCLUDED.
              !KC: This term needs to be included for USPP.
              !KC: add the augmentation charge term for dvscf
              call adddvscf (1, ik)
           else
             ! At the first iteration dvbare_q*psi_kpoint is calculated
             ! and written to file
             ! call dvqpsi_us (ik, u (1, mode),.false. )
               call dvqpsi_mag_us (ik, .false.)
             ! add the augmentation charge term for dvext and dbext
               call adddvscf (1, ik)

            IF(niter_ph>1)then
               IF(reduce_io) THEN
               dvpsi0(:,:,ik)=dvpsi(:,:)
               ELSE
               call davcio (dvpsi, lrbar, iubar, nrec, 1)
               END IF
            END IF
               
           endif
           !
           ! Ortogonalize dvpsi to valence states: ps = <evq|dvpsi>
           ! Apply -P_c^+.
           !
           cw = fiu(iw)   
           CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, cw)
           !
           if (where_rec=='solve_lint'.or.iter > 1) then
              !
              ! starting value for delta_psi is read from iudwf
              !
              IF(reduce_io)THEN
              dpsi(:,:)=dpsi0(:,:,ik)
              ELSE
              call davcio ( dpsi, lrdwf, iudwf, nrec, -1)
              END IF 
!For frequency dependent case we will require two more wave functions
 !             call davcio ( dpsip, lrdwf, iudwfp, nrec, -1)
 !             call davcio ( dpsim, lrdwf, iudwfm, nrec, -1)
 !              dpsi(:,:) = (0.d0, 0.d0)
              !
              ! threshold for iterative solution of the linear system
              !
              !thresh = min (1.d-1 * sqrt (dr2), 1.d-2)
              thresh = thresh_CG
           else
              !
              !  At the first iteration dpsi and dvscfin are set to zero
              !
 
              dpsi(:,:) = (10.d0, 1.d0)
!              dpsim(:,:)     = (0.d0, 0.d0)
!              dpsip(:,:)     = (0.d0, 0.d0)
              dvscfin (:, :) = (0.d0, 0.d0)
              dvscfout(:, :) = (0.d0, 0.d0)
              !
              ! starting threshold for iterative solution of the linear system
              !
              thresh = thresh_CG
              if(niter_ph==1)thresh = thresh_CG ! * sqrt (tr2_ph)
           endif

           !
           ! iterative solution of the linear system (H-eS)*dpsi=dvpsi,
           ! dvpsi=-P_c^+ (dvbare+dvscf)*psi , dvscf fixed.
           !
           conv_root = .true.
           etc(:,:) = CMPLX( et(:,:), 0.0d0 , kind=DP)
!           cw = fiu(iw)

           if(real(cw).eq.0.d0.and.aimag(cw).eq.0.d0)then

           call cgsolve_all (ch_psi_all, cg_psi, et(1,ikk), dvpsi, dpsi, &
                             h_diag, npwx, npwq, thresh, ik, lter, conv_root, &
                             anorm, nbnd_occ(ikk), npol )

           else
!           call cbcg_solve(cch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsi, h_diag, &
!                     npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol, cw, .true.,0)
           call cbicgstabl(cch_psi_all, cg_psi, etc(1,ikk), dvpsi, dpsi, h_diag, &
                     npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol, cw, 2 ,.false.)
           endif
                   

!           write(stdout,*)'cbcg_solve_end'
           ltaver = ltaver + lter
           lintercall = lintercall + 1
!           if (.not.conv_root) WRITE( stdout, '(5x,"kpoint",i4," ibnd",i4,  &
!                &              " solve_linter: root not converged ",e10.3)') &
!                &              ik , ibnd, anorm

           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
           !               if (nksq.gt.1 .or. npert(irr).gt.1)
!          call davcio (dpsip, lrdwf, iudwfp, nrec, + 1)
!          call davcio (epsim, lrdwf, iudwfm, nrec, + 1)
!       if(qpol ==1) then
IF(niter_ph>1)then          
           IF(reduce_io)THEN
           dpsi0(:,:,ik)=dpsi(:,:)
           ELSE
           call davcio (dpsi, lrdwf, iudwf, nrec, + 1)
           END IF
END IF
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           weight = wk (ikk)
           IF (noncolin) THEN
!              call incdrhoscf_nc(drhoscf(1,1),weight,ik, &
!                                 dbecsum_nc(1,1,1,1), dpsip,dpsim)

               call incdrhoscf_nc(drhoscf(1,1),weight,ik, &
                                 dbecsum_nc(1,1,1,1), dpsi)
           ELSE
              call incdrhoscf(drhoscf(1,current_spin), weight, ik, &
                               dbecsum(1,1,current_spin), dpsi)
           END IF
     enddo! on k-points
         
          
          
#ifdef __MPI
     
     !
     !  The calculation of dbecsum is distributed across processors (see addusdbec)
     !  Sum over processors the contributions coming from each slice of bands
     !
     IF (noncolin) THEN
        call mp_sum ( dbecsum_nc, intra_pool_comm )
     ELSE
        call mp_sum ( dbecsum, intra_pool_comm )
     ENDIF
#endif

     if (doublegrid) then
        do is = 1, nspin_mag
              call cinterpolate (drhoscfh(1,is), drhoscf(1,is), 1)
        enddo
     else
        call zcopy (nspin_mag*dfftp%nnr, drhoscf, 1, drhoscfh, 1)
     endif
     !
     !In the noncolinear, spin-orbit case rotate dbecsum
     ! write(*,*) 'okvan, doublegrid', okvan, doublegrid
     !
     IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum)
     !
     !  Now we compute for all perturbations the total charge and potential
     !
     call addusddens (drhoscfh, dbecsum, 0)

#ifdef __MPI
     !
     !   Reduce the delta rho across pools
     !
     call mp_sum ( drhoscfh, inter_pool_comm )

     if(my_image_id/=0)then
     drhoscfh =  CONJG(drhoscfh)
     end if


     call mp_sum(drhoscfh, inter_image_comm)

     if(my_image_id/=0)then
     drhoscfh =  CONJG(drhoscfh)
     end if

     

#endif



     !
     ! q=0 in metallic case deserve special care (e_Fermi can shift)
     !
!     IF (lmetq0) call ef_shift(drhoscfh,ldos,ldoss,dos_ef,irr,npe,.false.)
     !
     !   After the loop over the perturbations we have the linear change
     !   in the charge density for each mode of this representation.
     !   Here we symmetrize them ...
     !
  IF((.not. man_kpoints) .and. (.not. symoff)) then 
     IF (.not.lgamma_gamma) THEN
!#ifdef __MPI
!        call psymdvscf (1, 1, drhoscfh)
!        IF ( noncolin.and.domag ) &
!        CALL psym_dmag( 1, 1, drhoscfh)
!#else
!        call symdvscf (1, 1, drhoscfh)
!        IF ( noncolin.and.domag ) CALL sym_dmag( 1, 1, drhoscfh)
!        call psym_dmage(dvtosym)
!        call sym_dmage(dvtosym)
!#endif
!        IF (okpaw) THEN
!           IF (minus_q) CALL PAW_dumqsymmetrize(dbecsum,npe,irr, &
!                             npertx,irotmq,rtau,xq,tmq)
!           CALL  &
!              PAW_dusymmetrize(dbecsum,npe,irr,npertx,nsymq,rtau,xq,t)
!        END IF
!#ifdef __MPI
!        call psymb (drhoscfh)
!        IF ( noncolin.and.domag ) CALL psym_dmagb(drhoscfh)
!#else
        call symb (drhoscfh)
        IF ( noncolin.and.domag ) CALL sym_dmagb(drhoscfh)
!#endif
     ENDIF
  END IF
     !
     !   ... save them on disk and
     !   compute the corresponding change in scf potential
     !

     call zcopy (dfftp%nnr*nspin_mag, drhoscfh(1,1), 1, dvscfout(1,1), 1)
     call dv_of_drho (dvscfout(1,1), .false.)
     !
     !   And we mix with the old potential
     ! KC: test dvscfout
!write(stdout, '("dvscfin before mix", 2f10.6)') dvscfin(100,2)
!write(stdout, '("dvscfout before mix", 4f10.6)') dvscfout(100,1), dvscfout(100,2)

     !
!      if(real(cw).eq.0.d0.and.aimag(cw).eq.0.d0)then
!      call mix_potential (2*dfftp%nnr*nspin_mag, dvscfout, dvscfin, &
!                         alpha_mix(kter), dr2, tr2_ph/npol, iter, &
!                         nmix_ph, flmixdpot, convt)
!      else
!     if(my_image_id==0)then
      IF(niter_ph==1)alpha_mix=1.d0
      
      if(transverse .and. .not. do_elec)then
      call mix_potential_c(dfftp%nnr*3, dvscfout(1,2), dvscfin(1,2), &
                             alpha_mix(kter), dr2, tr2_ph/npol, iter, &
                             nmix_ph, convt)
      else

!      convt = .true.
!      convtm= .false.

!      do im=1,4,3
      
!      call mix_potential_c(dfftp%nnr, dvscfout(1,im), dvscfin(1,im), &
!                             alpha_mix(kter), dr2, tr2_ph/npol, iter, &
!                             nmix_ph, convtm)
!      convt= (convt .and. convtm)
!      end do
!      write(stdout, *)im

       if(real(cw).eq.0.d0.and.aimag(cw).eq.0.d0)then
       call mix_potential (2*dfftp%nnr*nspin_mag, dvscfout, dvscfin, &
                         alpha_mix(kter), dr2, tr2_ph/npol, iter, &
                         nmix_ph, flmixdpot, convt)
       else
       call mix_potential_c(dfftp%nnr*nspin_mag, dvscfout, dvscfin, &
                             alpha_mix(kter), dr2, tr2_ph/npol, iter, &
                             nmix_ph, convt)
       end if
!       dvscfout(:,2)=dvscfin(:,2)
!       dvscfout(:,3)=dvscfin(:,3)
      end if
!        if(convt)then 
!        convt_check=1
!        else
!        convt_check=0
!        end if
!     end if
!      if(my_image_id/=0)then
!         dvscfout(:,:)=CONJG(dvscfout(:,:))
!         dvscfin(:,:)=CONJG(dvscfin(:,:))
     if(my_image_id/=0)then
!         call mix_potential_c(dfftp%nnr*nspin_mag, dvscfout, dvscfin, &
!                             alpha_mix(kter), dr2, tr2_ph/npol, iter, &
!                             nmix_ph, convt) 
         dvscfin(:,:)=CONJG(dvscfin(:,:))
!         dvscfin(:,:)=(0.d0, 0.d0)
!         convt_check=0
     end if

      call mp_sum(dvscfin, inter_image_comm)
!      call mp_sum(convt_check,inter_image_comm)
   
     if(my_image_id/=0)then
         dvscfin(:,:)=CONJG(dvscfin(:,:)/2.d0)
     else
         dvscfin(:,:)=dvscfin(:,:)/2.d0
     end if
!         if(convt_check==1)then
!          convt=.true.
!         else
!          convt=.false.
!         end if         
!     end if 

!   write(stdout, '("dvscfin after mix", 4f10.6)') dvscfin(100,1), dvscfin(100,2)
!      CALL check_all_convt(convt, iter)

     if (lmetq0.and.convt) call ef_shift (drhoscfh, ldos, ldoss, dos_ef, irr, npe, .true.)


     ! check that convergent have been reached on ALL processors in this image
      CALL check_all_convt(convt)

     if (doublegrid) then
         do is = 1, nspin_mag
             call cinterpolate (dvscfin(1,is), dvscfins(1,is), -1)
         enddo
     endif

    call newdq (dvscfin, 1)  !KC: calculate int3 for dvscf

#ifdef __MPI
     aux_avg (1) = DBLE (ltaver)
     aux_avg (2) = DBLE (lintercall)
     call mp_sum ( aux_avg, inter_pool_comm )
     averlt = aux_avg (1)  / aux_avg (2)
#else
     averlt = DBLE (ltaver)  / lintercall
#endif
     tcpu = get_clock ('PHONON')

     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / npe
     WRITE( stdout, '(5x," thresh=",e10.3, " alpha_mix = ",f6.3, &
          &      " |ddv_scf|^2 = ",e10.3 )') thresh, alpha_mix (kter) , dr2
     !
     !    Here we save the information for recovering the run from this poin
     !
     CALL flush_unit( stdout )
     !
     rec_code=10
     CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, &
                                            dvscfin, drhoscfh)

     if (check_stop_now()) call stop_smoothly_ph (.false.)
!Well spotted Kun!
!   if (convt) then
     if (doublegrid) then
        do is = 1, nspin_mag
              call cinterpolate (drhoscfh(1,is), drhoscf(1,is), -1)
        enddo
     else
        call zcopy (nspin_mag*dfftp%nnr, drhoscfh, 1, drhoscf, 1)
     endif


     if(noncolin)then 

!#################################  for test purpose ########################

write(stdout,'(5x," q,freq,drho in real space, "5f10.5,"  ",8f10.4)') xq(:), fiu(iw)*13.6057, (real(drhoscf (dffts%nnr/2,ig)), ig = 1,4), &
                                                          (aimag(drhoscf(dffts%nnr/2,ig)),ig = 1,4)

!####################################################################


         do ig=1, nspin_mag
           CALL fwfft ('Smooth', drhoscf(:,ig), dffts)
         enddo
          
!         write(stdout,'(5x, "g(1), g(2), g(10) ",9f12.5)') g(1, 1), g(2, 1), g(3,1), g(1, 2), g(2, 2), g(3,2), g(1, 10), g(2, 10), g(3,10)
         write(stdout,'(5x," q,freq,drho, "5f10.5,"  ",8f10.4)') xq(:), fiu(iw)*13.6057, (real(drhoscf (nls(1),ig)), ig = 1,4), &
                                                          (aimag(drhoscf(nls(1),ig)), ig = 1,4)
         

!################################# for test purpose ##################################
!write(stdout,'(5x," q,freq,drho, "5f10.5,"  ",8f10.4)') xq(:),fiu(iw)*13.6057, (real(drhoscf (nls(2),ig)), ig = 1,4), &
!                                                          (aimag(drhoscf(nls(2),ig)), ig = 1,4)
!write(stdout,'(5x," q,freq,drho, "5f10.5,"  ",8f10.4)') xq(:),fiu(iw)*13.6057, (real(drhoscf (nls(10),ig)), ig = 1,4), &
!                                                          (aimag(drhoscf(nls(10),ig)),ig = 1,4)
!#######################################################################################


        If(convt) write(stdout,'(5x," q,freq,convtdrho, "5f10.5,"  ",8f10.4)') xq(:), fiu(iw)*13.6057, (real(drhoscf (nls(1),ig)), ig = 1,4), &
                                                          (aimag(drhoscf(nls(1),ig)),ig = 1,4)

        if(transverse) then
        WRITE(stdout, '(5x,"transverse magnetic response" )')
        if(convt)then
          write(stdout,'(5x,"convtchiq+-", 3f10.4, f12.5,"  ",4f14.7)') xq(:), real(fiu(iw))*13605.7, &
          real((drhoscf(1,2)+(0.d0,1.d0)*drhoscf(1,3))/2.d0), &
          aimag((drhoscf(1,2)+(0.d0,1.d0)*drhoscf(1,3))/2.d0),&
          real((drhoscf(1,2)-(0.d0,1.d0)*drhoscf(1,3))/2.d0), &
          aimag((drhoscf(1,2)-(0.d0,1.d0)*drhoscf(1,3))/2.d0)
        else
          write(stdout,'(5x,"freq, chiq+-, "f12.5,"  ",4f14.7)') real(fiu(iw))*13605.7, &
          real((drhoscf(1,2)+(0.d0,1.d0)*drhoscf(1,3))/2.d0), &
          aimag((drhoscf(1,2)+(0.d0,1.d0)*drhoscf(1,3))/2.d0),&
          real((drhoscf(1,2)-(0.d0,1.d0)*drhoscf(1,3))/2.d0), &
          aimag((drhoscf(1,2)-(0.d0,1.d0)*drhoscf(1,3))/2.d0)
        end if
        end if
     end if


!######################## for test purpose  ################################

!      IF(do_elec) then
!          write(stdout,'(5x," real space q,freq,drho, "5f10.5,"  ",4f10.6)') xq(:), fiu(iw)*13.6057, drhoscf (dffts%nnr/2,1), drhoscf (dffts%nnr/4,1)

!          CALL fwfft ('Smooth', drhoscf(:,1), dffts)

!          write(stdout,'(5x," PW space q,freq,drho, "5f10.5,"  ",6f10.6)') xq(:), fiu(iw)*13.6057, drhoscf (nls(1),1), drhoscf (nls(2),1), drhoscf (nls(3),1)

!      END IF

!#################################################################
!         do ig=1, nspin_mag

     if(do_elec) then
           IF(nspin_mag==2)then
           drhoscf(:,1) = (dvscfins(:,1)+dvscfins(:,2))/2
           ELSE
           drhoscf(:,1) = dvscfins(:,1)
           END IF

           CALL fwfft ('Smooth', drhoscf(:,1), dffts)
!         enddo
          if(convt)then
          write(stdout,'(5x,"q,freq,real(convteps),im(convteps) "5f12.5," ",4f14.7)') xq(:), fiu(iw), &
          real((1.0/(1.0d0 + drhoscf (nls(1),1)))), aimag((1.0/(1.0d0+ drhoscf (nls(1),1)))), &
          (1.0d0 + real(drhoscf (nls(1),1))), aimag(drhoscf(nls(1),1))
          else
          write(stdout,'(5x,"q,freq,real(eps),im(eps) "5f12.5," ",4f14.7)') xq(:), fiu(iw), &
          real((1.0/(1.0d0 + drhoscf (nls(1),1)))), aimag((1.0/(1.0d0+ drhoscf(nls(1),1)))), &
          (1.0d0 + real(drhoscf (nls(1),1))), aimag(drhoscf(nls(1),1))
          end if
     end if
     
!   end if
     if (convt) goto 155


  enddo    ! end of iter
155 iter0=0

!    if(convt)
  
  !
  !    A part of the dynamical matrix requires the integral of
  !    the self consistent change of the potential and the variation of
  !    the charge due to the displacement of the atoms.
  !    We compute it here.
  !


  if (allocated(ldoss)) deallocate (ldoss)
  if (allocated(ldos)) deallocate (ldos)
  deallocate (h_diag)
  deallocate (aux1)
  deallocate (dbecsum)
  IF (okpaw) THEN
     if (lmetq0.and.allocated(becsum1)) deallocate (becsum1)
     deallocate (mixin)
     deallocate (mixout)
  ENDIF
  IF (noncolin) deallocate (dbecsum_nc)
  deallocate (dvscfout)
  deallocate (drhoscfh)

  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)

!  deallocate (dpsip)
!  deallocate (dpsim)

  call stop_clock ('solve_linter')
END SUBROUTINE solve_linter


SUBROUTINE check_all_convt(convt, iter)
  USE mp,        ONLY : mp_sum, mp_barrier
  USE mp_global, ONLY : nproc_image, me_image, intra_image_comm, inter_image_comm
  USE control_ph,           ONLY : niter_ph
  USE io_global,            ONLY : stdout
  IMPLICIT NONE
  LOGICAL,INTENT(inout) :: convt
  INTEGER:: convt_check, iter
  !
  ! IF(nproc_image==1) RETURN
  !
  ! ALLOCATE(convt_check(nproc_image+1))
  !
  convt_check = 1
  IF(convt) convt_check = 0
  !
  CALL mp_sum(convt_check, intra_image_comm)
  IF(convt .and. convt_check/=0) THEN
    CALL errore('check_all_convt', 'Only some of the processors within each image converged: '&
               &' something is wrong with solve_linter', 1)
  ENDIF
!  CALL mp_barrier(inter_image_comm)
  CALL mp_sum(convt_check, inter_image_comm)
  !CALL mp_sum(ios, inter_pool_comm)
  !CALL mp_sum(ios, intra_pool_com)
  !
 ! convt = (convt_check==0)
!  IF(ANY(convt_check==0).and..not.ALL(convt_check==0) ) THEN
  IF(convt .and. convt_check/=0) THEN
!  CALL errore('check_all_convt', '+q and -q did not converge with the same footing: '&
!               &' something is wrong with solve_linter', 1)
  IF(iter/=niter_ph)then
  write(stdout,'(5x, "+q and -q did not converge with the same footing. Continue &
               & until they all converges")')
  convt = .false.
  else
  write(stdout,'(5x, "+q or -q did not converge.")')
!  CALL errore('check_all_convt', '+q or -q did not converge: '&
!               &' something is wrong with solve_linter', 1)
  ENDIF
  ENDIF
  !
  RETURN
  !
END SUBROUTINE
