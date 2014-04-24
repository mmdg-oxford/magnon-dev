!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
Module dynamical
  !
  ! All variables read from file that need dynamical allocation
  !
  USE kinds, ONLY: DP
  complex(DP), allocatable :: dyn(:,:,:,:)
  real(DP), allocatable :: tau(:,:),  zstar(:,:,:), dchi_dtau(:,:,:,:), &
                           m_loc(:,:)
  integer, allocatable :: ityp(:)
  !
end Module dynamical
!
!--------------------------------------------------------------------
      program dynmat
!--------------------------------------------------------------------
!
!  This program
!  - reads a dynamical matrix file produced by the phonon code
!  - adds the nonanalytical part (if Z* and epsilon are read from file),
!    applies the chosen Acoustic Sum Rule (if q=0)
!  - diagonalise the dynamical matrix
!  - calculates IR and Raman cross sections (if Z* and Raman tensors
!    are read from file, respectively)
!  - writes the results to files, both for inspection and for plotting
!
!  Input data (namelist "input")
!
!  fildyn  character input file containing the dynamical matrix
!                    (default: fildyn='matdyn')
!  q(3)      real    calculate LO modes (add nonanalytic terms) along
!                    the direction q (cartesian axis, default: q=(0,0,0) )
!  amass(nt) real    mass for atom type nt, amu
!                    (default: amass is read from file fildyn)
!  asr   character   indicates the type of Acoustic Sum Rule imposed
!                     - 'no': no Acoustic Sum Rules imposed (default)
!                     - 'simple':  previous implementation of the asr used
!                     (3 translational asr imposed by correction of
!                     the diagonal elements of the dynamical matrix)
!                     - 'crystal': 3 translational asr imposed by optimized
!                     correction of the dyn. matrix (projection).
!                     - 'one-dim': 3 translational asr + 1 rotational asr
!                     imposed by optimized correction of the dyn. mat. (the
!                     rotation axis is the direction of periodicity; it
!                     will work only if this axis considered is one of
!                     the cartesian axis).
!                     - 'zero-dim': 3 translational asr + 3 rotational asr
!                     imposed by optimized correction of the dyn. mat.
!                     Note that in certain cases, not all the rotational asr
!                     can be applied (e.g. if there are only 2 atoms in a
!                     molecule or if all the atoms are aligned, etc.).
!                     In these cases the supplementary asr are cancelled
!                     during the orthonormalization procedure (see below).
!                     Finally, in all cases except 'no' a simple correction
!                     on the effective charges is performed (same as in the
!                     previous implementation).
!  axis    integer    indicates the rotation axis for a 1D system
!                     (1=Ox, 2=Oy, 3=Oz ; default =3)
!  filout  character output file containing frequencies and modes
!                    (default: filout='dynmat.out')
!  filmol  character as above, in a format suitable for 'molden'
!                    (default: filmol='dynmat.mold')
!  filxsf  character as above, in axsf format suitable for xcrysden
!                    (default: filmol='dynmat.axsf')
!
      USE kinds, ONLY: DP
      USE mp,         ONLY : mp_bcast
      USE mp_global,  ONLY : mp_startup, mp_global_end
      USE io_global,  ONLY : ionode, ionode_id, stdout
      USE environment, ONLY : environment_start, environment_end
      USE io_dyn_mat,  ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                             read_dyn_mat, read_dyn_mat_tail
      USE constants,   ONLY : amu_ry
      use dynamical
      !
      implicit none
      integer, parameter :: ntypx = 10
      character(len=256):: fildyn, filout, filmol, filxsf
      character(len=3) :: atm(ntypx)
      character(len=10) :: asr
      logical :: lread, gamma
      complex(DP), allocatable :: z(:,:)
      real(DP) :: amass(ntypx), amass_(ntypx), eps0(3,3), a0, omega, &
           at(3,3), bg(3,3), q(3), q_(3)
      real(DP), allocatable :: w2(:)
      integer :: nat, na, nt, ntyp, iout, axis, nax, nspin_mag, ios
      real(DP) :: celldm(6)
      logical :: xmldyn, lrigid, lraman
      logical, external :: has_xml
      integer :: ibrav, nqs
      integer, allocatable :: itau(:)
      namelist /input/ amass, asr, axis, fildyn, filout, filmol, filxsf, q
      !
      ! code is parallel-compatible but not parallel
      !
      CALL mp_startup()
      CALL environment_start('DYNMAT')
      !
      IF (ionode) CALL input_from_file ( )
      !
      asr  = 'no'
      axis = 3
      fildyn='matdyn'
      filout='dynmat.out'
      filmol='dynmat.mold'
      filxsf='dynmat.axsf'
      amass(:)=0.0d0
      q(:)=0.0d0
      !
      IF (ionode) read (5,input, iostat=ios)
      CALL mp_bcast(ios, ionode_id)
      CALL errore('dynmat', 'reading input namelist', ABS(ios))
      !
      CALL mp_bcast(asr,ionode_id)
      CALL mp_bcast(axis,ionode_id)
      CALL mp_bcast(amass,ionode_id)
      CALL mp_bcast(fildyn,ionode_id)
      CALL mp_bcast(filout,ionode_id)
      CALL mp_bcast(filmol,ionode_id)
      CALL mp_bcast(filxsf,ionode_id)
      CALL mp_bcast(q,ionode_id)
      !
      IF (ionode) inquire(file=fildyn,exist=lread)
      CALL mp_bcast(lread, ionode_id)
      IF (lread) THEN
         IF (ionode) WRITE(6,'(/5x,a,a)') 'Reading Dynamical Matrix from file '&
                                         , TRIM(fildyn)
      ELSE
         CALL errore('dynmat', 'File '//TRIM(fildyn)//' not found', 1)
      END IF
      !
      ntyp = ntypx ! avoids spurious out-of-bound errors
      xmldyn=has_xml(fildyn)
      IF (xmldyn) THEN
         CALL read_dyn_mat_param(fildyn,ntyp,nat)
         ALLOCATE (m_loc(3,nat))
         ALLOCATE (tau(3,nat))
         ALLOCATE (ityp(nat))
         ALLOCATE (zstar(3,3,nat))
         ALLOCATE (dchi_dtau(3,3,3,nat) )
         CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
                 celldm, at, bg, omega, atm, amass_, tau, ityp, &
                 m_loc, nqs, lrigid, eps0, zstar, lraman, dchi_dtau)
         IF (nqs /= 1) CALL errore('dynmat','only q=0 matrix allowed',1)
         a0=celldm(1) ! define alat
         at = at / a0 !  bring at in units of alat
         ALLOCATE (dyn(3,3,nat,nat) )
         CALL read_dyn_mat(nat,1,q_,dyn(:,:,:,:))
         CALL read_dyn_mat_tail(nat)
         if(asr.ne.'no') then
             call set_asr ( asr, axis, nat, tau, dyn, zstar )
         endif
         IF (ionode) THEN
            do nt=1, ntyp
               if (amass(nt) <= 0.0d0) amass(nt)=amass_(nt)
            end do
         END IF
      ELSE
         IF (ionode) THEN
            call readmat ( fildyn, asr, axis, nat, ntyp, atm, a0, &
                            at, omega, amass_, eps0, q_ )
            do nt=1, ntyp
               if (amass(nt) <= 0.0d0) amass(nt)=amass_(nt)/amu_ry
            end do
         END IF
      ENDIF
      !
      IF (ionode) THEN
         !
         ! from now on, execute on a single processor
         !
         gamma = ( abs( q_(1)**2+q_(2)**2+q_(3)**2 ) < 1.0d-8 )
         !
         IF (gamma) THEN
            allocate (itau(nat))
            do na=1,nat
               itau(na)=na
            end do
            call nonanal ( nat, nat, itau, eps0, q, zstar,omega, dyn )
            deallocate (itau)
         END IF
         !
         nax = nat
         allocate ( z(3*nat,3*nat), w2(3*nat) )
         call dyndiag(nat,ntyp,amass,ityp,dyn,w2,z)
         !
         if (filout.eq.' ') then
            iout=6
         else
            iout=4
            OPEN (unit=iout,file=filout,status='unknown',form='formatted')
         end if
         call writemodes(nax,nat,q_,w2,z,iout)
         if(iout .ne. 6) close(unit=iout)
         call writemolden (filmol, gamma, nat, atm, a0, tau, ityp, w2, z)
         call writexsf (filxsf, gamma, nat, atm, a0, at, tau, ityp, z)
         IF (gamma) call RamanIR &
           (nat, omega, w2, z, zstar, eps0, dchi_dtau)
      ENDIF
      !
      IF (xmldyn) THEN
         DEALLOCATE (m_loc)
         DEALLOCATE (tau)
         DEALLOCATE (ityp)
         DEALLOCATE (zstar)
         DEALLOCATE (dchi_dtau)
         DEALLOCATE (dyn)
      ENDIF

      CALL environment_end('DYNMAT')
      !
      CALL mp_global_end()

      end program dynmat
!
!-----------------------------------------------------------------------
      subroutine readmat ( fildyn, asr, axis, nat, ntyp, atm, &
           a0, at, omega, amass, eps0, q )
!-----------------------------------------------------------------------
!
      USE kinds, ONLY: DP
      use dynamical
      !
      implicit none
      character(len=256), intent(in) :: fildyn
      character(len=10), intent(in) :: asr
      integer, intent(in) :: axis
      integer, intent(inout) :: nat, ntyp
      character(len=3), intent(out) ::  atm(ntyp)
      real(DP), intent(out) :: amass(ntyp), a0, at(3,3), omega, &
           eps0(3,3), q(3)
      !
      character(len=80) :: line
      real(DP) :: celldm(6), dyn0r(3,3,2)
      integer :: ibrav, nt, na, nb, naa, nbb, i, j, k, ios
      logical :: qfinito, noraman
      !
      !
      noraman=.true.
      open (unit=1,file=fildyn,status='old',form='formatted')
      read(1,'(a)') line
      read(1,'(a)') line
      read(1,*) ntyp,nat,ibrav,celldm
      !
      if (ibrav==0) then
         read(1,'(a)') line
         read(1,*) ((at(i,j),i=1,3),j=1,3)
      end if
      !
      allocate ( dyn (3,3,nat,nat) )
      allocate ( dchi_dtau (3,3,3,nat) )
      allocate (zstar(3,3,nat) )
      allocate ( tau (3,nat) )
      allocate (ityp (nat) )
      !
      call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
      a0=celldm(1) ! define alat
      at = at / a0 !  bring at in units of alat

      do nt=1,ntyp
         read(1,*) i,atm(nt),amass(nt)
      end do
      do na=1,nat
         read(1,*) i,ityp(na), (tau(j,na),j=1,3)
      end do
      read(1,'(a)') line
      read(1,'(a)') line
      read(1,'(a)') line
      read(1,'(a)') line
      read(line(11:80),*) (q(i),i=1,3)
      qfinito=q(1).ne.0.0 .or. q(2).ne.0.0 .or. q(3).ne.0.0
      if (qfinito .and. asr .ne. 'no') &
           call errore('readmat','Acoustic Sum Rule for q != 0 ?',1)
      do na = 1,nat
         do nb = 1,nat
            read (1,*) naa, nbb
            if (na.ne.naa .or. nb.ne.nbb) then
               call errore ('readmat','mismatch in reading file',1)
            end if
            read (1,*) ((dyn0r(i,j,1), dyn0r(i,j,2), j=1,3), i=1,3)
            dyn(:,:,na,nb) = CMPLX( dyn0r(:,:,1), dyn0r(:,:,2) ,kind=DP)
         end do
      end do
      write(6,'(5x,a)') '...Force constants read'
      !
      if (.not.qfinito) then
         ios=0
         read(1,*,iostat=ios)
         read(1,'(a)',iostat=ios) line
         if (ios .ne. 0 .or. line(1:23).ne.'     Dielectric Tensor:') then
            write(6,'(5x,a)') '...epsilon and Z* not read (not found on file)'
            do na=1,nat
               do j=1,3
                  do i=1,3
                     zstar(i,j,na)=0.0d0
                  end do
               end do
            end do
            do j=1,3
               do i=1,3
                  eps0(i,j)=0.0d0
               end do
               eps0(j,j)=1.0d0
            end do
         else
            read(1,*)
            read(1,*) ((eps0(i,j), j=1,3), i=1,3)
            read(1,*)
            read(1,*)
            read(1,*)
            do na = 1,nat
               read(1,*)
               read(1,*) ((zstar(i,j,na), j=1,3),i=1,3)
            end do
            write(6,'(5x,a)') '...epsilon and Z* read'
 20         read(1,'(a)',end=10,err=10) line
            if (line(1:10) == '     Raman') go to 25
            go to 20
 25         read(1,*,end=10,err=10)
            do na = 1,nat
               do i = 1, 3
                  read(1,*,end=10,err=10)
                  read(1,*,end=10,err=10) &
                       ((dchi_dtau(k,j,i,na), j=1,3), k=1,3)
               end do
            end do
            write(6,'(5x,a)') '...Raman cross sections read'
            noraman=.false.
10          continue
         end if
      end if
      if (noraman) dchi_dtau=0.d0
!
      if(asr.ne.'no') then
          call set_asr ( asr, axis, nat, tau, dyn, zstar )
      endif
!
      close(unit=1)
!
      return
      end subroutine readmat
!
!-----------------------------------------------------------------------
subroutine RamanIR (nat, omega, w2, z, zstar, eps0, dchi_dtau)
  !-----------------------------------------------------------------------
  !
  !   write IR and Raman cross sections
  !   on input: z = eigendisplacements (normalized as <z|M|z>)
  !             zstar = effective charges (units of e)
  !             dchi_dtau = derivatives of chi wrt atomic displacement
  !                         (units: A^2)
 USE kinds, ONLY: DP
 USE constants, ONLY : fpi, BOHR_RADIUS_ANGS, RY_TO_THZ, RY_TO_CMM1, amu_ry
 implicit none
 ! input
 integer, intent(in) :: nat
 real(DP) omega, w2(3*nat), zstar(3,3,nat), eps0(3,3), &
      dchi_dtau(3,3,3,nat), chi(3,3)
 complex(DP) z(3*nat,3*nat)
 ! local
 integer na, nu, ipol, jpol, lpol
 logical noraman
 real(DP), allocatable :: infrared(:), raman(:,:,:)
 real(DP):: polar(3), cm1thz, freq, irfac
 real(DP):: cmfac, alpha, beta2
 !
 !
 cm1thz = RY_TO_THZ/RY_TO_CMM1
 !
 !   conversion factor for IR cross sections from
 !   (Ry atomic units * e^2)  to  (Debye/A)^2/amu
 !   1 Ry mass unit = 2 * mass of one electron = 2 amu
 !   1 e = 4.80324x10^(-10) esu = 4.80324 Debye/A
 !     (1 Debye = 10^(-18) esu*cm = 0.2081928 e*A)
 !
 irfac = 4.80324d0**2/2.d0*amu_ry
 !
 write (6,'(/5x,"Polarizability (A^3 units)")')
 !
 !  correction to molecular polarizabilities from Clausius-Mossotti formula
 !  (for anisotropic systems epsilon is replaced by its trace)
 !
 cmfac = 3.d0 / ( 2.d0 + (eps0(1,1) + eps0(2,2) + eps0(3,3))/3.d0 )
 !
 write (6,'(5x,"multiply by",f9.6," for Clausius-Mossotti correction")') cmfac
 do jpol=1,3
    do ipol=1,3
       if (ipol == jpol) then
          chi(ipol,jpol) = (eps0(ipol,jpol)-1.d0)
       else
          chi(ipol,jpol) = eps0(ipol,jpol)
       end if
    end do
 end do
 do ipol=1,3
    write (6,'(5x,3f12.6)') (chi(ipol,jpol)*BOHR_RADIUS_ANGS**3*omega/fpi, &
         jpol=1,3)
 end do
 !
 allocate(infrared (3*nat))
 allocate(raman(3,3,3*nat))
 !
 noraman=.true.
 do nu = 1,3*nat
    do ipol=1,3
       polar(ipol)=0.0d0
    end do
    do na=1,nat
       do ipol=1,3
          do jpol=1,3
             polar(ipol) = polar(ipol) +  &
                  zstar(ipol,jpol,na)*z((na-1)*3+jpol,nu)
          end do
       end do
    end do
    !
    infrared(nu) = 2.d0*(polar(1)**2+polar(2)**2+polar(3)**2)*irfac
    !
    do ipol=1,3
       do jpol=1,3
          raman(ipol,jpol,nu)=0.0d0
          do na=1,nat
             do lpol=1,3
                raman(ipol,jpol,nu) = raman(ipol,jpol,nu) + &
                     dchi_dtau(ipol,jpol,lpol,na) * z((na-1)*3+lpol,nu)
             end do
          end do
          noraman=noraman .and. abs(raman(ipol,jpol,nu)).lt.1.d-12
       end do
    end do
    !   Raman cross sections are in units of bohr^4/(Ry mass unit)
 end do
 !
 write (6,'(/5x,"IR activities are in (D/A)^2/amu units")')
 if (noraman) then
    write (6,'(/"# mode   [cm-1]    [THz]      IR")')
 else
    write (6,'(5x,"Raman activities are in A^4/amu units")')
    write (6,'(5x,"multiply Raman by",f9.6," for Clausius-Mossotti", &
         & " correction")') cmfac**2
    write (6,'(/"# mode   [cm-1]    [THz]      IR          Raman   depol.fact")')
 end if
 !
 do nu = 1,3*nat
    !
    freq = sqrt(abs(w2(nu)))*RY_TO_CMM1
    if (w2(nu).lt.0.0) freq = -freq
    !
    ! alpha, beta2: see PRB 54, 7830 (1996) and refs quoted therein
    !
    if (noraman) then
       write (6,'(i5,f10.2,2f10.4)') &
         nu, freq, freq*cm1thz, infrared(nu)
    else
       alpha = (raman(1,1,nu) + raman(2,2,nu) + raman(3,3,nu))/3.d0
       beta2 = ( (raman(1,1,nu) - raman(2,2,nu))**2 + &
                 (raman(1,1,nu) - raman(3,3,nu))**2 + &
                 (raman(2,2,nu) - raman(3,3,nu))**2 + 6.d0 * &
          (raman(1,2,nu)**2 + raman(1,3,nu)**2 + raman(2,3,nu)**2) )/2.d0
       write (6,'(i5,f10.2,2f10.4,f15.4,f10.4)') &
            nu, freq, freq*cm1thz, infrared(nu), &
            (45.d0*alpha**2 + 7.0d0*beta2)*amu_ry, &
             3.d0*beta2/(45.d0*alpha**2 + 4.0d0*beta2)
    end if
 end do
 !
 deallocate (raman)
 deallocate (infrared)
 return
 !
end subroutine RamanIR
!
!----------------------------------------------------------------------
subroutine set_asr ( asr, axis, nat, tau, dyn, zeu )
  !-----------------------------------------------------------------------
  !
  !  Impose ASR - refined version by Nicolas Mounet
  !
  USE kinds, ONLY: DP
  implicit none
  character(len=10), intent(in) :: asr
  integer, intent(in) :: axis, nat
  real(DP), intent(in) :: tau(3,nat)
  real(DP), intent(inout) :: zeu(3,3,nat)
  complex(DP), intent(inout) :: dyn(3,3,nat,nat)
  !
  integer :: i,j,n,m,p,k,l,q,r,na, nb, na1, i1, j1
  real(DP), allocatable:: dynr_new(:,:,:,:,:), zeu_new(:,:,:)
  real(DP), allocatable :: u(:,:,:,:,:)
  ! These are the "vectors" associated with the sum rules
  !
  integer u_less(6*3*nat),n_less,i_less
  ! indices of the vectors u that are not independent to the preceding ones,
  ! n_less = number of such vectors, i_less = temporary parameter
  !
  integer, allocatable :: ind_v(:,:,:)
  real(DP), allocatable :: v(:,:)
  ! These are the "vectors" associated with symmetry conditions, coded by
  ! indicating the positions (i.e. the four indices) of the non-zero elements
  ! (there should be only 2 of them) and the value of that element.
  ! We do so in order to use limit the amount of memory used.
  !
  real(DP), allocatable :: w(:,:,:,:), x(:,:,:,:)
  real(DP) sum, scal, norm2
  ! temporary vectors and parameters
  !
  real(DP), allocatable :: zeu_u(:,:,:,:)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer zeu_less(6*3),nzeu_less,izeu_less
  ! indices of the vectors zeu_u that are not independent to the preceding
  ! ones, nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(DP), allocatable :: zeu_w(:,:,:), zeu_x(:,:,:)
  ! temporary vectors
  !
  ! Initialization
  ! n is the number of sum rules to be considered (if asr.ne.'simple')
  ! 'axis' is the rotation axis in the case of a 1D system (i.e. the rotation
  !  axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if axis='3')
  !
  if ( (asr.ne.'simple') .and. (asr.ne.'crystal') .and. (asr.ne.'one-dim') &
                         .and.(asr.ne.'zero-dim')) then
     call errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
  endif
  if(asr.eq.'crystal') n=3
  if(asr.eq.'one-dim') then
     write(6,'("asr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(asr.eq.'zero-dim') n=6
  !
  ! ASR on effective charges
  !
  if(asr.eq.'simple') then
     do i=1,3
        do j=1,3
           sum=0.0d0
           do na=1,nat
              sum = sum + zeu(i,j,na)
           end do
           do na=1,nat
              zeu(i,j,na) = zeu(i,j,na) - sum/nat
           end do
        end do
     end do
  else
     ! generating the vectors of the orthogonal of the subspace to project
     ! the effective charges matrix on
     !
     allocate ( zeu_new(3,3,nat) )
     allocate (zeu_u(6*3,3,3,nat) )
     zeu_u(:,:,:,:)=0.0d0
     do i=1,3
        do j=1,3
           do na=1,nat
              zeu_new(i,j,na)=zeu(i,j,na)
           enddo
        enddo
     enddo
     !
     p=0
     do i=1,3
        do j=1,3
           ! These are the 3*3 vectors associated with the
           ! translational acoustic sum rules
           p=p+1
           zeu_u(p,i,j,:)=1.0d0
           !
        enddo
     enddo
     !
     if (n.eq.4) then
        do i=1,3
           ! These are the 3 vectors associated with the
           ! single rotational sum rule (1D system)
           p=p+1
           do na=1,nat
              zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
              zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
           enddo
           !
        enddo
     endif
     !
     if (n.eq.6) then
        do i=1,3
           do j=1,3
              ! These are the 3*3 vectors associated with the
              ! three rotational sum rules (0D system - typ. molecule)
              p=p+1
              do na=1,nat
                 zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
                 zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
              enddo
              !
           enddo
        enddo
     endif
     !
     ! Gram-Schmidt orthonormalization of the set of vectors created.
     !
     allocate ( zeu_w(3,3,nat), zeu_x(3,3,nat) )
     nzeu_less=0
     do k=1,p
        zeu_w(:,:,:)=zeu_u(k,:,:,:)
        zeu_x(:,:,:)=zeu_u(k,:,:,:)
        do q=1,k-1
           r=1
           do izeu_less=1,nzeu_less
              if (zeu_less(izeu_less).eq.q) r=0
           enddo
           if (r.ne.0) then
              call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
              zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
           endif
        enddo
        call sp_zeu(zeu_w,zeu_w,nat,norm2)
        if (norm2.gt.1.0d-16) then
           zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
        else
           nzeu_less=nzeu_less+1
           zeu_less(nzeu_less)=k
        endif
     enddo
     !
     !
     ! Projection of the effective charge "vector" on the orthogonal of the
     ! subspace of the vectors verifying the sum rules
     !
     zeu_w(:,:,:)=0.0d0
     do k=1,p
        r=1
        do izeu_less=1,nzeu_less
           if (zeu_less(izeu_less).eq.k) r=0
        enddo
        if (r.ne.0) then
           zeu_x(:,:,:)=zeu_u(k,:,:,:)
           call sp_zeu(zeu_x,zeu_new,nat,scal)
           zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
        endif
     enddo
     !
     ! Final substraction of the former projection to the initial zeu, to get
     ! the new "projected" zeu
     !
     zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
     call sp_zeu(zeu_w,zeu_w,nat,norm2)
     write(6,'(5x,"Acoustic Sum Rule: || Z*(ASR) - Z*(orig)|| = ",E15.6)') &
          SQRT(norm2)
     !
     ! Check projection
     !
     !write(6,'("Check projection of zeu")')
     !do k=1,p
     !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
     !  call sp_zeu(zeu_x,zeu_new,nat,scal)
     !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
     !enddo
     !
     do i=1,3
        do j=1,3
           do na=1,nat
              zeu(i,j,na)=zeu_new(i,j,na)
           enddo
        enddo
     enddo
     deallocate (zeu_w, zeu_x)
     deallocate (zeu_u)
     deallocate (zeu_new)
  endif
  !
  ! ASR on dynamical matrix
  !
  if(asr.eq.'simple') then
     do i=1,3
        do j=1,3
           do na=1,nat
              sum=0.0d0
              do nb=1,nat
                 if (na.ne.nb) sum=sum + DBLE (dyn(i,j,na,nb))
              end do
              dyn(i,j,na,na) = CMPLX(-sum, 0.d0,kind=DP)
           end do
        end do
     end do
     !
  else
     ! generating the vectors of the orthogonal of the subspace to project
     ! the dyn. matrix on
     !
     allocate (u(6*3*nat,3,3,nat,nat))
     allocate (dynr_new(2,3,3,nat,nat))
     u(:,:,:,:,:)=0.0d0
     do i=1,3
        do j=1,3
           do na=1,nat
              do nb=1,nat
                 dynr_new(1,i,j,na,nb) = DBLE (dyn(i,j,na,nb) )
                 dynr_new(2,i,j,na,nb) =AIMAG (dyn(i,j,na,nb) )
              enddo
           enddo
        enddo
     enddo
     !
     p=0
     do i=1,3
        do j=1,3
           do na=1,nat
              ! These are the 3*3*nat vectors associated with the
              ! translational acoustic sum rules
              p=p+1
              do nb=1,nat
                 u(p,i,j,na,nb)=1.0d0
              enddo
              !
           enddo
        enddo
     enddo
     !
     if (n.eq.4) then
        do i=1,3
           do na=1,nat
              ! These are the 3*nat vectors associated with the
              ! single rotational sum rule (1D system)
              p=p+1
              do nb=1,nat
                 u(p,i,axis,na,nb)=0.0d0
                 u(p,i,MOD(axis,3)+1,na,nb)=-tau(MOD(axis+1,3)+1,nb)
                 u(p,i,MOD(axis+1,3)+1,na,nb)=tau(MOD(axis,3)+1,nb)
              enddo
              !
           enddo
        enddo
     endif
     !
     if (n.eq.6) then
        do i=1,3
           do j=1,3
              do na=1,nat
                 ! These are the 3*3*nat vectors associated with the
                 ! three rotational sum rules (0D system - typ. molecule)
                 p=p+1
                 do nb=1,nat
                    u(p,i,j,na,nb)=0.0d0
                    u(p,i,MOD(j,3)+1,na,nb)=-tau(MOD(j+1,3)+1,nb)
                    u(p,i,MOD(j+1,3)+1,na,nb)=tau(MOD(j,3)+1,nb)
                 enddo
                 !
              enddo
           enddo
        enddo
     endif
     !
     allocate (ind_v(9*nat*nat,2,4))
     allocate (v(9*nat*nat,2))
     m=0
     do i=1,3
        do j=1,3
           do na=1,nat
              do nb=1,nat
                 ! These are the vectors associated with the symmetry constraints
                 q=1
                 l=1
                 do while((l.le.m).and.(q.ne.0))
                    if ((ind_v(l,1,1).eq.i).and.(ind_v(l,1,2).eq.j).and. &
                         (ind_v(l,1,3).eq.na).and.(ind_v(l,1,4).eq.nb)) q=0
                    if ((ind_v(l,2,1).eq.i).and.(ind_v(l,2,2).eq.j).and. &
                         (ind_v(l,2,3).eq.na).and.(ind_v(l,2,4).eq.nb)) q=0
                    l=l+1
                 enddo
                 if ((i.eq.j).and.(na.eq.nb)) q=0
                 if (q.ne.0) then
                    m=m+1
                    ind_v(m,1,1)=i
                    ind_v(m,1,2)=j
                    ind_v(m,1,3)=na
                    ind_v(m,1,4)=nb
                    v(m,1)=1.0d0/DSQRT(2.0d0)
                    ind_v(m,2,1)=j
                    ind_v(m,2,2)=i
                    ind_v(m,2,3)=nb
                    ind_v(m,2,4)=na
                    v(m,2)=-1.0d0/DSQRT(2.0d0)
                 endif
              enddo
           enddo
        enddo
     enddo
     !
     ! Gram-Schmidt orthonormalization of the set of vectors created.
     ! Note that the vectors corresponding to symmetry constraints are already
     ! orthonormalized by construction.
     !
     allocate ( w(3,3,nat,nat), x(3,3,nat,nat))
     n_less=0
     do k=1,p
        w(:,:,:,:)=u(k,:,:,:,:)
        x(:,:,:,:)=u(k,:,:,:,:)
        do l=1,m
           !
           call sp2(x,v(l,:),ind_v(l,:,:),nat,scal)
           do r=1,2
              i=ind_v(l,r,1)
              j=ind_v(l,r,2)
              na=ind_v(l,r,3)
              nb=ind_v(l,r,4)
              w(i,j,na,nb)=w(i,j,na,nb)-scal*v(l,r)
           enddo
        enddo
        if (k.le.(9*nat)) then
           na1=MOD(k,nat)
           if (na1.eq.0) na1=nat
           j1=MOD((k-na1)/nat,3)+1
           i1=MOD((((k-na1)/nat)-j1+1)/3,3)+1
        else
           q=k-9*nat
           if (n.eq.4) then
              na1=MOD(q,nat)
              if (na1.eq.0) na1=nat
              i1=MOD((q-na1)/nat,3)+1
           else
              na1=MOD(q,nat)
              if (na1.eq.0) na1=nat
              j1=MOD((q-na1)/nat,3)+1
              i1=MOD((((q-na1)/nat)-j1+1)/3,3)+1
           endif
        endif
        do q=1,k-1
           r=1
           do i_less=1,n_less
              if (u_less(i_less).eq.q) r=0
           enddo
           if (r.ne.0) then
              call sp3(x,u(q,:,:,:,:),i1,na1,nat,scal)
              w(:,:,:,:) = w(:,:,:,:) - scal* u(q,:,:,:,:)
           endif
        enddo
        call sp1(w,w,nat,norm2)
        if (norm2.gt.1.0d-16) then
           u(k,:,:,:,:) = w(:,:,:,:) / DSQRT(norm2)
        else
           n_less=n_less+1
           u_less(n_less)=k
        endif
     enddo
     !
     ! Projection of the dyn. "vector" on the orthogonal of the
     ! subspace of the vectors verifying the sum rules and symmetry contraints
     !
     w(:,:,:,:)=0.0d0
     do l=1,m
        call sp2(dynr_new(1,:,:,:,:),v(l,:),ind_v(l,:,:),nat,scal)
        do r=1,2
           i=ind_v(l,r,1)
           j=ind_v(l,r,2)
           na=ind_v(l,r,3)
           nb=ind_v(l,r,4)
           w(i,j,na,nb)=w(i,j,na,nb)+scal*v(l,r)
        enddo
     enddo
     do k=1,p
        r=1
        do i_less=1,n_less
           if (u_less(i_less).eq.k) r=0
        enddo
        if (r.ne.0) then
           x(:,:,:,:)=u(k,:,:,:,:)
           call sp1(x,dynr_new(1,:,:,:,:),nat,scal)
           w(:,:,:,:) = w(:,:,:,:) + scal* u(k,:,:,:,:)
        endif
     enddo
     !
     ! Final substraction of the former projection to the initial dyn,
     ! to get the new "projected" dyn
     !
     dynr_new(1,:,:,:,:)=dynr_new(1,:,:,:,:) - w(:,:,:,:)
     call sp1(w,w,nat,norm2)
     write(6,'(5x,"Acoustic Sum Rule: ||dyn(ASR) - dyn(orig)||= ",E15.6)') &
          DSQRT(norm2)
     !
     ! Check projection
     !
     !write(6,'("Check projection")')
     !do l=1,m
     !  call sp2(dynr_new(1,:,:,:,:),v(l,:),ind_v(l,:,:),nat,scal)
     !  if (DABS(scal).gt.1d-10) write(6,'("l= ",I8," dyn|v(l)= ",F15.10)') l,scal
     !enddo
     !do k=1,p
     !  x(:,:,:,:)=u(k,:,:,:,:)
     !  call sp1(x,dynr_new(1,:,:,:,:),nat,scal)
     !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," dyn|u(k)= ",F15.10)') k,scal
     !enddo
     !
     deallocate ( w, x )
     deallocate ( v )
     deallocate ( ind_v )
     deallocate ( u )
     !
     do i=1,3
        do j=1,3
           do na=1,nat
              do nb=1,nat
                 dyn (i,j,na,nb) = &
                      CMPLX(dynr_new(1,i,j,na,nb), dynr_new(2,i,j,na,nb) ,kind=DP)
              enddo
           enddo
        enddo
     enddo
     deallocate ( dynr_new )
  endif
  !
  return
end subroutine set_asr
!
!
!----------------------------------------------------------------------
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  integer i,j,na,nat
  real(DP) zeu_u(3,3,nat)
  real(DP) zeu_v(3,3,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp_zeu
!
!
!----------------------------------------------------------------------
subroutine sp1(u,v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two dyn. matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  integer i,j,na,nb,nat
  real(DP) u(3,3,nat,nat)
  real(DP) v(3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        do nb=1,nat
          scal=scal+u(i,j,na,nb)*v(i,j,na,nb)
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp1
!
!----------------------------------------------------------------------
subroutine sp2(u,v,ind_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two dyn. matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat) space). u is coded in the usual way
  ! but v is coded as explained when defining the vectors corresponding to the
  ! symmetry constraints
  !
  USE kinds, ONLY: DP
  implicit none
  integer i,nat
  real(DP) u(3,3,nat,nat)
  integer ind_v(2,4)
  real(DP) v(2)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,2
    scal=scal+u(ind_v(i,1),ind_v(i,2),ind_v(i,3),ind_v(i,4))*v(i)
  enddo
  !
  return
  !
end subroutine sp2
!
!----------------------------------------------------------------------
subroutine sp3(u,v,i,na,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! like sp1, but in the particular case when u is one of the u(k)%vec
  ! defined in set_asr (before orthonormalization). In this case most of the
  ! terms are zero (the ones that are not are characterized by i and na), so
  ! that a lot of computer time can be saved (during Gram-Schmidt).
  !
  USE kinds, ONLY: DP
  implicit none
  integer i,j,na,nb,nat
  real(DP) u(3,3,nat,nat)
  real(DP) v(3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do j=1,3
    do nb=1,nat
      scal=scal+u(i,j,na,nb)*v(i,j,na,nb)
    enddo
  enddo
  !
  return
  !
end subroutine sp3
