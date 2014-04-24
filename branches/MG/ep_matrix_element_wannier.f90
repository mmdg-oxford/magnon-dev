!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE ep_matrix_element_wannier()
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation from data saved in fildvscf
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : celldm, omega, ibrav
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, amass
  USE gvecs, ONLY: doublegrid
  USE fft_base, ONLY : dfftp, dffts
  USE noncollin_module, ONLY : nspin_mag
  USE dynmat, ONLY : dyn, w2
  USE qpoint, ONLY : xq, nksq, ikks
  USE modes,  ONLY : npert, nirr
  USE control_ph, ONLY : trans
  USE units_ph, ONLY : iudyn, lrdrho, iudvscf
  USE io_global, ONLY : stdout
  USE mp_global, ONLY : me_pool, root_pool
  USE modes, ONLY : u
  USE klist, ONLY : xk
  USE wvfct, ONLY : npwx
  USE el_phon, ONLY: elph_mat, kpq, g_kpq, igqg, xk_gamma
  !
  IMPLICIT NONE
  !
  LOGICAL :: read_dvscf_cart, ascii_dvscf
  INTEGER :: irr, imode0, ipert, is, ik
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  COMPLEX(DP), POINTER :: dvscfin(:,:,:), dvscfins (:,:,:)


  CALL start_clock ('elphon')


  ascii_dvscf=.false.
  if(elph_mat) read_dvscf_cart=.true.
  if(read_dvscf_cart) then
     write(stdout,*)
     write(stdout,*) 'Reading dvscf in cartesian coordinates !'
     write(stdout,*)
     
     u=CMPLX(0.d0,0.d0)
     do irr=1,3*nat
        u(irr,irr)=CMPLX(1.d0,0.d0)
     enddo
     

!     if(ascii_dvscf) then
!        ALLOCATE (dvrot ( nrxx , nspin , 3*nat) )
!        fildvscf_asc=trim(tmp_dir)//trim(prefix)//"."//trim(fildvscf)//'1'
!        open(unit=7899,file=fildvscf_asc,status='unknown')
!        dvrot=CMPLX(0.0,0.0)
!        do na=1,nat
!           do ipol=1,3
!              irr=(na-1)*3+ipol
!              do  k = 1, dfftp%nr3
!                 do j = 1, dfftp%nr2
!                    do i = 1, dfftp%nr1
!                       read(7899,*)   n, rep,imp
!                       dvrot(n,1,irr)=CMPLX(rep,imp)
!                    enddo
!                 enddo
!              enddo
!           enddo
!        enddo
!        close(7899)
!     endif
     
  endif

  

  !
  ! read Delta Vscf and calculate electron-phonon coefficients
  !
  imode0 = 0
  DO irr = 1, nirr
     ALLOCATE (dvscfin (dfftp%nnr, nspin_mag , npert(irr)) )
!     if(ascii_dvscf) then
!        DO ipert = 1, npert(irr)
!           dvscfin(1:dfftp%nnr,:,ipert)=dvrot(1:dfftp%nnr,:,imode0+ipert)
!        enddo
!     else
     DO ipert = 1, npert (irr)
        CALL davcio_drho ( dvscfin(1,1,ipert),  lrdrho, iudvscf, &
             imode0 + ipert,  -1 )
     END DO
!     endif
     IF (doublegrid) THEN
        ALLOCATE (dvscfins (dffts%nnr, nspin_mag , npert(irr)) )
        DO is = 1, nspin_mag
           DO ipert = 1, npert(irr)
              CALL cinterpolate (dvscfin(1,is,ipert),dvscfins(1,is,ipert),-1)
           ENDDO
        ENDDO
     ELSE
        dvscfins => dvscfin
     ENDIF
     CALL newdq (dvscfin, npert(irr))
     CALL elphel_refolded (npert (irr), imode0, dvscfins)
     !
     imode0 = imode0 + npert (irr)
     IF (doublegrid) DEALLOCATE (dvscfins)
     DEALLOCATE (dvscfin)
  ENDDO
  !
  ! now read the eigenvalues and eigenvectors of the dynamical matrix
  ! calculated in a previous run
  !
  IF (.NOT.trans) CALL readmat (iudyn, ibrav, celldm, nat, ntyp, &
       ityp, omega, amass, tau, xq, w2, dyn)
  !
  

  deallocate(xk_gamma)
  deallocate(kpq,g_kpq,igqg)

! CALL stop_clock ('elphon')
  RETURN
END SUBROUTINE ep_matrix_element_wannier

!-----------------------------------------------------------------------
SUBROUTINE elphsum_wannier(q_index)
  !-----------------------------------------------------------------------
  !
  !      Sum over BZ of the electron-phonon matrix elements el_ph_mat
  !      Original routine written by Francesco Mauri
  !      Adapted to wannier functions by Matteo Calandra
  !            Dev. Comment: missing calc_sigma_yet
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat, ityp, tau,amass,tau, ntyp => nsp, atm
  USE cell_base, ONLY : at, bg, ibrav, celldm 
  USE symm_base, ONLY : s, sr, irt, nsym, time_reversal, invs, ftau, copy_sym, inverse_s
  USE klist, ONLY : xk, nelec
  USE wvfct, ONLY : nbnd, et
  USE el_phon
  USE mp_global, ONLY : me_pool, root_pool, inter_pool_comm, npool, intra_pool_comm
  USE io_global, ONLY : stdout,ionode
  USE io_files,  ONLY : prefix
  USE qpoint, ONLY : xq, nksq
  USE dynmat, ONLY : dyn, w2
  USE modes, ONLY : u, gi, gimq, irgq, irotmq
  USE control_ph, only : lgamma
  USE lsda_mod, only : isk,nspin, current_spin,lsda
  USE mp,        ONLY: mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER :: q_index
  !
  !
  logical :: minus_qloc,sym (48)
  integer :: nq, imq, isq(48)
  INTEGER :: ik, ikk, ikq, ibnd, jbnd, ipert, jpert, nu, mu, &
         ios, iuelphmat,i,j,nt,k
  INTEGER :: iu_sym,nmodes,nsymq
  INTEGER :: io_file_unit
  !   for star_q
  REAL(DP) :: rtauloc(3,48,nat)
  !   end of star_q definitions
  real(DP) :: sxq (3, 48)
  REAL(DP) xk_dummy(3)
  character(len=80) :: filelph
  CHARACTER(len=256) ::  file_elphmat
  !
  COMPLEX(DP) :: el_ph_sum (3*nat,3*nat), dyn_corr(3*nat,3*nat)

  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER (LEN=6), EXTERNAL :: int_to_char

  nmodes=3*nat

  write(filelph,'(A5,f9.6,A1,f9.6,A1,f9.6)') 'elph.',xq(1),'.',xq(2),'.',xq(3)
  file_elphmat=trim(adjustl(prefix))//'_elph.mat.q_'// TRIM( int_to_char( q_index ) )

  
  ! parallel case: only first node writes
  IF ( .not.ionode ) THEN
     iuelphmat = 0
  ELSE
     !
     ! First I dump information for the electron-phonon interaction
     !

     !
     iuelphmat = find_free_unit()
     OPEN (unit = iuelphmat, file = file_elphmat, status = 'unknown', err = &
          111, iostat = ios, form='unformatted')
111  CALL errore ('elphon', 'opening file'//file_elphmat, ABS (ios) )
     REWIND (iuelphmat)
     xk_dummy(:)=xq(:)
     call cryst_to_cart(1,xk_dummy,at,-1)
     WRITE (iuelphmat) xk_dummy
     WRITE (iuelphmat) nelec
     WRITE (iuelphmat) elph_nbnd_min,elph_nbnd_max,nbnd
     WRITE (iuelphmat) nmodes, nksq, nat, ntyp
     WRITE (iuelphmat) ibrav,(celldm(j),j=1,6)
     WRITE(iuelphmat)  (atm(j),j=1,ntyp),(amass(j),j=1,ntyp), &
             (ityp(j),j=1,nat),((tau(j,i),j=1,3),i=1,nat)
     WRITE (iuelphmat) (w2 (nu) , nu = 1, nmodes)
     WRITE (iuelphmat) ((u(ipert,jpert),ipert=1,nmodes),jpert=1,nmodes)
     WRITE (iuelphmat) ((dyn(ipert,jpert),ipert=1,3*nat),jpert=1,3*nat)
     do ik=1,nksq
        IF (lgamma) THEN
           ikk = ik
           ikq = ik
        ELSE
           ikk = 2 * ik - 1
           ikq = ikk + 1
        ENDIF
        xk_dummy(:)=xk(:,ikk)
        call cryst_to_cart(1,xk_dummy,at,-1)
        WRITE (iuelphmat) (xk_dummy(ipert),ipert=1,3)
        WRITE (iuelphmat) (et(ibnd,ikk),ibnd=elph_nbnd_min,elph_nbnd_max)

        do nu=1,nmodes
           WRITE (iuelphmat) &
                ((el_ph_mat (jbnd, ibnd, ik, nu),jbnd=elph_nbnd_min,elph_nbnd_max),&
                ibnd=elph_nbnd_min,elph_nbnd_max)
        enddo
     enddo
  


     !
     ! Then I dump symmetry operations
     !
     minus_qloc = .true.
     sym = .false.
     sym(1:nsym) = .true.
     
     call smallg_q (xq, 0, at, bg, nsym, s, ftau, sym, minus_qloc)
     nsymq = copy_sym(nsym, sym)
     ! recompute the inverses as the order of sym.ops. has changed
     CALL inverse_s ( )
     
     ! part 2: this redoes most of the above, plus it computes irgq, gi, gimq
     CALL smallgq (xq, at, bg, s, nsym, irgq, nsymq, irotmq, &
          minus_qloc, gi, gimq)
     
     sym(1:nsym)=.true.
     call sgam_ph (at, bg, nsym, s, irt, tau, rtauloc, nat, sym)
     call star_q(xq, at, bg, nsym , s , invs , nq, sxq, &
          isq, imq, .FALSE. )
     
     
     do j=1,3
        write(iuelphmat) (at(i,j),i=1,3)
     enddo
     do j=1,3
        write(iuelphmat) (bg(i,j),i=1,3)
     enddo
     write(iuelphmat) nsym,nq,imq
     do i=1,nsym
        write(iuelphmat) i,invs(i),isq(i)
        do j=1,3
           do k=1,3
              write(iuelphmat) k,j, s(k,j,i)
           enddo
        enddo
        do j=1,nat
           write(iuelphmat) j, irt(i,j) 
        enddo
        do j=1,3
           do k=1,nat
              write(iuelphmat) j,i, rtauloc(j,i,k)  
           enddo
        enddo
        do j=1,3
           write(iuelphmat) j, sxq(j,i)
        enddo
     enddo
     
     close(iuelphmat)
  endif
  
  !
  !
  

  RETURN
  

END SUBROUTINE ELPHSUM_WANNIER
   
!
!-----------------------------------------------------------------------
SUBROUTINE elphel_refolded (npe, imode0, dvscfins)
  !-----------------------------------------------------------------------
  !
  !      Calculation of the electron-phonon matrix elements el_ph_mat
  !         <\psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)>
  !      Original routine written by Francesco Mauri
  !
  USE kinds, ONLY : DP
  USE fft_base, ONLY : dffts
  USE wavefunctions_module,  ONLY: evc
  USE io_files, ONLY: iunigk, prefix, diropn
  USE klist, ONLY: xk
  USE lsda_mod, ONLY: lsda, current_spin, isk
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE wvfct, ONLY: nbnd, npw, npwx, igk
  USE uspp, ONLY : vkb
  USE el_phon, ONLY : el_ph_mat, iunwfcwann, igqg, kpq, g_kpq, &
           xk_gamma, npwq_refolded, lrwfcr
  USE modes, ONLY : u
  USE units_ph, ONLY : iubar, lrbar, lrwfc, iuwfc
  USE eqv,      ONLY : dvpsi!, evq
  USE qpoint,   ONLY : igkq, npwq, nksq, ikks, ikqs
  USE control_ph, ONLY : trans, lgamma
  USE mp_global, ONLY: intra_pool_comm, me_pool, root_pool
  USE mp,        ONLY: mp_sum
  USE ions_base, ONLY : nat
  USE io_global, ONLY : stdout
  USE gvecs, ONLY : nls


  IMPLICIT NONE
  !
  INTEGER :: npe, imode0
  COMPLEX(DP) :: dvscfins (dffts%nnr, nspin_mag, npe)
  COMPLEX(DP), allocatable :: evq(:,:)


  ! LOCAL variables
  logical :: exst
  INTEGER :: nrec, ik, ikk, ikq, ikqg,ipert, mode, ibnd, jbnd, ir, ig, &
       ios

  COMPLEX(DP) , ALLOCATABLE :: aux1 (:,:), elphmat (:,:,:)
  COMPLEX(DP), EXTERNAL :: zdotc
  INTEGER, EXTERNAL :: find_free_unit
  !
  allocate (evq(npol*npwx,nbnd))
  ALLOCATE (aux1    (dffts%nnr, npol))
  ALLOCATE (elphmat ( nbnd , nbnd , 3*nat))

  
! iunwfcwann=find_free_unit()
! CALL diropn (iunwfcwann, 'wfc', lrwfc, exst, dvscf_dir)
! IF (.NOT.exst) THEN
!    CALL errore ('elphel_refolded', 'file '//trim(prefix)//'.wfc not found in Rotated_DVSCF', 1)
! END IF
  !
  !  Start the loops over the k-points
  !
  IF (nksq.GT.1) REWIND (unit = iunigk)

  
  DO ik = 1, nksq
     
     IF (nksq.GT.1) THEN
        READ (iunigk, err = 100, iostat = ios) npw, igk
100     CALL errore ('elphel', 'reading igk', ABS (ios) )
     ENDIF
     !
     !  ik = counter of k-points with vector k
     !  ikk= index of k-point with vector k
     !  ikq= index of k-point with vector k+q
     !       k and k+q are alternated if q!=0, are the same if q=0
     !
     IF (lgamma) npwq = npw
     ikk = ikks(ik)
     ikq = ikqs(ik)
     ikqg = kpq(ik)


     IF (lsda) current_spin = isk (ikk)
     IF (.NOT.lgamma.AND.nksq.GT.1) THEN
        READ (iunigk, err = 200, iostat = ios) npwq, igkq
200     CALL errore ('elphel', 'reading igkq', ABS (ios) )
     ENDIF

     

     !
     CALL init_us_2 (npwq, igkq, xk (1, ikq), vkb)
     !
     ! read unperturbed wavefuctions psi(k) and psi(k+q)
     !
     evc=cmplx(0.d0,0.d0)

!    Warning error in reading wfc, this could explain.
!    We read here the wfc at the Gamma point, that is
!    that saved by Wannier.

!     CALL davcio (evc, lrwfc, iunwfcwann, ik, - 1)
!     CALL davcio (evq, lrwfc, iunwfcwann, ikqg, - 1)



!     IF (nksq.GT.1) THEN
!        IF (lgamma) THEN
!           CALL davcio (evc, lrwfc, iunwfcwann, ikk, - 1)
!        ELSE
!           CALL davcio (evc, lrwfc, iunwfcwann, ik, - 1)
!           CALL davcio (evq, lrwfc, iunwfcwann, ikqg, - 1)
!        ENDIF
!     ENDIF
     !

     call read_wfc_rspace_and_fwfft( evc , ik , lrwfcr , iunwfcwann , npw , igk )


     call calculate_and_apply_phase(ik, ikqg, igqg, npwq_refolded, g_kpq,xk_gamma, evq, .true.)
     

     DO ipert = 1, npe
        nrec = (ipert - 1) * nksq + ik
        !
        !  dvbare_q*psi_kpoint is read from file (if available) or recalculated
        !
        IF (trans) THEN
           CALL davcio (dvpsi, lrbar, iubar, nrec, - 1)
        ELSE
           mode = imode0 + ipert
           ! TODO : .false. or .true. ???
           CALL dvqpsi_us (ik, u (1, mode), .FALSE. )
        ENDIF
        !
        ! calculate dvscf_q*psi_k
        !

        DO ibnd = 1, nbnd
           CALL cft_wave (evc(1, ibnd), aux1, +1)
           CALL apply_dpot(dffts%nnr, aux1, dvscfins(1,1,ipert), current_spin)
           CALL cft_wave (dvpsi(1, ibnd), aux1, -1)
        END DO
        CALL adddvscf (ipert, ik)
        !
        ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
        !
        DO ibnd =1, nbnd
           DO jbnd = 1, nbnd
              elphmat (jbnd, ibnd, ipert) = zdotc (npwq_refolded, evq (1, jbnd), 1, &
                   dvpsi (1, ibnd), 1)
              IF (noncolin) &
                 elphmat (jbnd, ibnd, ipert) = elphmat (jbnd, ibnd, ipert)+ &
                   zdotc (npwq_refolded, evq(npwx+1,jbnd),1,dvpsi(npwx+1,ibnd), 1)
           ENDDO
        ENDDO
     ENDDO
     !
     CALL mp_sum (elphmat, intra_pool_comm)
     !
     !  save all e-ph matrix elements into el_ph_mat
     !
     DO ipert = 1, npe
        DO jbnd = 1, nbnd
           DO ibnd = 1, nbnd
              el_ph_mat (ibnd, jbnd, ik, ipert + imode0) = elphmat (ibnd, jbnd, ipert)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  

!  CLOSE( UNIT = iunwfcwann, STATUS = 'KEEP' )
  !
  DEALLOCATE (elphmat)
  DEALLOCATE (aux1)
  DEALLOCATE(evq)
  !
  RETURN
END SUBROUTINE elphel_refolded
!
subroutine get_equivalent_kpq(xk,xq,kpq,g_kpq, igqg)
    !==================================================================!
    !                                                                  !
    !  Set up the k+q shell for electron-phonon coupling               ! 
    !                                                                  !
    !  This routine finds the G vectors such that                      !
    !  k+q+G=k'  with k and k' belonging to nksq                       !
    !  for each k, the G vector is stored in g_kpq                     !
    !  k'=kpq(ik)                                                      !
    !  and finally igqg(ik) is the index that allows to find           !
    !  the g vector g_kpq in the list of all the G vectors             !
    !                                                                  !
    !       Matteo Calandra                                            !
    !===================================================================  
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE cell_base,   ONLY : at, bg
  USE qpoint, ONLY : nksq, ikks
  USE gvect, ONLY: g, gg
  USE qpoint, ONLY : nksq
  USE mp_global, ONLY : intra_pool_comm
  USE mp, ONLY : mp_sum
  ! WARNING g_kpq mesh is an integer
  implicit none
  
  ! Variables that are private
  
  integer :: iqx,iqy,iqz,i,j,k,n,nn,iq,ik, ig
  integer :: kpq(nksq),g_kpq(3,nksq),igqg(nksq)
  integer, allocatable :: ig_check(:)
  real(kind=dp) :: gg_
  real(kind=dp) :: xq(3), xk(3,*)
  real(kind=dp) :: xkpq(3),Gvec(3),xq_crys(3)
  real(kind=dp), allocatable :: xk_crys(:,:)

  !
  ! nksq = number of k point per pool withour k+q
  !

  ! 
  ! The xk_point entering here must be the k and not
  !    the k+q
  !

  xq_crys=xq

  call cryst_to_cart (1, xq_crys, at, -1)


  allocate(xk_crys(3,nksq))

  do ik=1,nksq
     xk_crys(1:3,ik)=xk(1:3,ik)
  enddo
  call cryst_to_cart (nksq, xk_crys, at, -1)

  !
  ! kpt_latt are the BZ vectors in crystalline coordinates
  ! xq is the q vector in crystalline coordinates
  !
  
  do iq=1,nksq
     xkpq(:)=xk_crys(:,iq)+xq_crys(:)
     do i=1,nksq
        do iqx=-4,4
           do iqy=-4,4
              do iqz=-4,4
                 Gvec(1)=real(iqx,dp)+xkpq(1)
                 Gvec(2)=real(iqy,dp)+xkpq(2)
                 Gvec(3)=real(iqz,dp)+xkpq(3)
                   
                 if(dabs(xk_crys(1,i)-Gvec(1)).lt.1.d-6.and. &
                      dabs(xk_crys(2,i)-Gvec(2)).lt.1.d-6.and. &
                      dabs(xk_crys(3,i)-Gvec(3)).lt.1.d-6) then
                    kpq(iq)=i
                    g_kpq(1,iq)=-iqx
                    g_kpq(2,iq)=-iqy
                    g_kpq(3,iq)=-iqz
                    goto 99
                 endif
              enddo
           enddo
        enddo
     enddo
     CALL errore ('get_equivalent_kpq', 'cannot find index k+q ', 2 )
     stop
99   continue
  enddo

  !
  ! here between all the g-vectors I find the index of that
  ! related to the translation in the Brillouin zone.
  ! Warning: if G does not belong to the processor igqg is zero.
  !

  igqg=0
  do ik=1,nksq
     Gvec(:) = REAL( g_kpq(:,ik),dp )
     call cryst_to_cart (1, Gvec, bg, 1)
     gg_ = Gvec(1)*Gvec(1) + Gvec(2)*Gvec(2) + Gvec(3)*Gvec(3)
     igqg(ik)=0
     ig=1
     do while  (gg(ig) <= gg_ + 1.d-6) 
        if ( (abs(g(1,ig)-Gvec(1)) < 1.d-6) .and.  &
             (abs(g(2,ig)-Gvec(2)) < 1.d-6) .and.  &
             (abs(g(3,ig)-Gvec(3)) < 1.d-6)  ) then
           igqg(ik) = ig
           
        endif
        ig= ig +1
     end do
  end do

  allocate(ig_check(nksq))
  ig_check=igqg
  CALL mp_sum( ig_check, intra_pool_comm )
  do ik=1,nksq
     if(ig_check(ik).eq.0) &
          CALL errore('get_equivalent_kpq', &
          ' g_kpq vector is not in the list of Gs', 100*ik )

  enddo
  deallocate(xk_crys)

  
end subroutine get_equivalent_kpq

subroutine calculate_and_apply_phase(ik, ikqg, igqg, npwq_refolded, g_kpq, xk_gamma, evq, lread)
  USE kinds, ONLY : DP
  USE fft_base, ONLY : dffts
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE wvfct, ONLY: nbnd, npw, npwx,  g2kin, ecutwfc, nbnd
  USE gvect, ONLY : ngm, g
  USE gvecs, ONLY : nls
  USE cell_base, ONLY : bg, tpiba2
  USE qpoint, ONLY : nksq, npwq
   USE wavefunctions_module, ONLY : evc
!  USE eqv,      ONLY : evq
  USE noncollin_module,     ONLY : npol
  USE el_phon, ONLY:iunwfcwann, lrwfcr

  IMPLICIT NONE

  LOGICAL :: lread
  INTEGER :: ik, ikqg,   npwq_refolded
  INTEGER :: igqg(nksq)
  INTEGER :: g_kpq(3,nksq)
  REAL (DP) :: xk_gamma(3,nksq)
  complex(dp) :: evq(npwx*npol,nbnd)
!  internal
  INTEGER :: npw_, m,i
  INTEGER, allocatable :: igk_(:), igkq_(:)
  REAL(DP) :: xkqg(3), g_(3), g_scra(3,ngm)
  COMPLEX (DP), allocatable :: psi_scratch(:)
  complex(DP), allocatable :: phase(:)

  allocate(igk_(npwx), igkq_(npwx))
  allocate (psi_scratch ( dffts%nnr) )
  allocate (phase(dffts%nnr))
  call flush_unit (6)

  g_scra=g

  g_(:)=real( g_kpq(:,ik), dp )
  call cryst_to_cart (1, g_, bg, 1)
  xkqg(:)=xk_gamma(:,ikqg)+g_(:)

  npw_=0
  npwq_refolded=0
  igk_=0
  igkq_=0


  call gk_sort (xk_gamma(1,ikqg), ngm, g_scra, ecutwfc / tpiba2, npw_, igk_, g2kin)

  if(lread) then
     call read_wfc_rspace_and_fwfft( evq , ikqg , lrwfcr , iunwfcwann , npw_ , igk_ )
  endif

  call gk_sort (xkqg, ngm, g_scra, ecutwfc / tpiba2, npwq_refolded, igkq_, g2kin)

  phase(:) = CMPLX(0.d0,0.d0)

  if ( igqg(ik)>0) then
     phase( nls(igqg(ik)) ) = (1.d0,0.d0)
  endif


  CALL invfft ('Wave', phase, dffts)
  !  call cft3s (phase, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
  phase(:)=conjg(phase(:))


  if(npwq_refolded.ne.npw_) call errore('elphel_refolded', 'Warning : npwq_refolded \= npw_',-1)
  
  do m=1,nbnd
     psi_scratch = (0.d0, 0.d0)
!     psi_scratch(nls (igk_ (1:npw_) ) ) = evq (1:npw_, m)
     psi_scratch(nls (igk_ (1:npw) ) ) = evq (1:npw, m)
     CALL invfft ('Wave', psi_scratch, dffts)
     !     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
     psi_scratch(1:dffts%nnr) = psi_scratch(1:dffts%nnr) * phase(1:dffts%nnr)
     !     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
     CALL fwfft ('Wave', psi_scratch, dffts)
     evq(1:npwq_refolded,m) = psi_scratch(nls (igkq_(1:npwq_refolded) ) )
  enddo


  deallocate(psi_scratch)
  DEALLOCATE(phase)
  deallocate(igk_, igkq_)

  return
end subroutine calculate_and_apply_phase
  
