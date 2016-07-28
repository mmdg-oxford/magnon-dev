!
! Copyright (C) 2001-2014 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE dfpt_tetra_mod
  !--------------------------------------------------------------------------
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC dfpt_tetra_beta, dfpt_tetra_dbl, dfpt_tetra_occ, dfpt_tetra_dlta, &
  &      dfpt_tetra_main
  !
  REAL(dp),ALLOCATABLE,SAVE :: &
  & dfpt_tetra_beta(:,:,:), & ! (nbnd, nbnd, nksq)
  & dfpt_tetra_dbl(:,:,:), & ! (nbnd, nbnd, nksq) Double occupancy
  & dfpt_tetra_occ(:,:), & ! (nbnd, nksq) Occupation
  & dfpt_tetra_dlta(:,:) ! (nbnd, nksq) delta(e_F - e_k)
  !
  INTEGER,SAVE :: &
  & sttk, & ! the first k in this Prosess Element
  & lstk, & ! the last k in this PE
  & nkBZ, & ! nk1 * nk2 * nk3
  & ivvec(3,20,6), & ! 
  & nsymqbz, & ! # of q-conserving symmetries of BZ; 
  & nsymbz ! # of symmetries of BZ
  !
  REAL(dp),SAVE :: &
  & wlsm(4,20)
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & grid(:,:), &   ! (3,nkBZ)
  & indx(:,:,:), & ! (nk1,nk2,nk3)
  & symq(:,:,:), & ! (3,3,nsymqbz)
  & sym(:,:,:)     ! (3,3,nsymbz)
  !
  REAL(dp),ALLOCATABLE,SAVE :: &
  & eig1(:,:), & ! (nbnd,nkBZ) eig(k)
  & eig2(:,:), & ! (nbnd,nkBZ) eig(k+q)
  & occ(:,:),  & ! (nbnd,nkBZ) occupations
  & dfde(:,:), & ! (nbnd,nkBZ) df / de
  & beta(:,:,:)  ! (nbnd,nbnd,nkBZ) 
  !
CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_main()
  !----------------------------------------------------------------------------
  !
  ! This routine computes the integration weight which is represented in Eq. (B 28)
  ! in PRB 64, 235118 (2001).
  !
  USE kinds,      ONLY : dp
  USE mp,         ONLY : mp_sum, mp_max
  USE mp_global,   ONLY : mpime, nproc, world_comm
  USE io_global,  ONLY : stdout
  USE symm_base,  ONLY : s, nsym, invsym
  USE wvfct,      ONLY : nbnd, wg
  USE klist,      ONLY : xk, wk, nks
  USE cell_base,  ONLY : at
  USE qpoint,     ONLY : nksq, xq
  USE control_ph, ONLY : lgamma, alpha_pv, nbnd_occ
  USE modes, ONLY : nsymq, minus_q
  USE start_k,    ONLY : nk1, nk2, nk3
  !
  IMPLICIT NONE
  !
  INTEGER :: nkpp, rest, ikv2(3), ik, ik2, ikk, ikq, ib, isym, nksym
  REAL(dp) :: kv1(3), kv2(3)
  !
  ! Symmetries of BZ of this crystal
  !
  nsymbz = nsym
  IF(.NOT. invsym) nsymbz = nsym * 2
  !
  ALLOCATE(sym(3,3,nsymbz))
  !
  sym(1:3,1:3,1:nsym) = s(1:3,1:3,1:nsym)
  IF(.NOT. invsym) sym(1:3,1:3,nsym+1:nsym+nsym) = - s(1:3,1:3,1:nsym)  
  !
  ! Symmetries of this q
  !
  nsymqbz = nsymq
  IF(minus_q) nsymqbz = nsymq * 2
  !
  ALLOCATE(symq(3,3,nsymqbz))
  !
  symq(1:3,1:3,1:nsymq) = s(1:3,1:3,1:nsymq)
  IF(minus_q) symq(1:3,1:3,nsymq+1:nsymq+nsymq) = - s(1:3,1:3,1:nsymq)
  !
  ! Compute eig1, eig2, indx, grid
  !
  CALL dfpt_tetra_eig()
  !
  ALLOCATE(dfde(nbnd,nkBZ), occ(nbnd,nkBZ), beta(nbnd,nbnd,nkBZ))
  !
  ! Work Sharing
  !
  nkpp = nkBZ / nproc
  rest = mod(nkBZ, nproc)
  IF(mpime < rest) THEN
     sttk = (nkpp + 1) *  mpime + 1
     lstk = (nkpp + 1) * (mpime + 1)
  ELSE 
     sttk = nkpp *  mpime + 1  + rest
     lstk = nkpp * (mpime + 1) + rest
  END IF
  !
  ! Define type of tetra
  !
  CALL dfpt_tetra_init()
  !
  !alpha_pv = - MINVAL(eig1(1:nbnd,1:nkBZ))
   alpha_pv = 0.d0
  !
  ! Calculate occupation & its derivative
  !
  CALL dfpt_tetra_dlta_occ()
  !
  ! Integration weight(1st & 2nd term)
  !
  beta(1:nbnd,1:nbnd,1:nkBZ) = 0d0
  CALL dfpt_tetra_calc_beta1()
  WRITE(stdout,*) "[dfpt_tetra]  beta1 : ", SUM(beta)
  CALL dfpt_tetra_calc_beta2()
  WRITE(stdout,*) "[dfpt_tetra]  beta2 : ", SUM(beta)
  !
  CALL mp_sum(beta(1:nbnd,1:nbnd,1:nkBZ), world_comm)
  !
  dfde(       1:nbnd,1:nkBZ) = dfde(        1:nbnd,1:nkBZ) / REAL(6 * nkBZ, dp)
  occ(        1:nbnd,1:nkBZ) =  occ(        1:nbnd,1:nkBZ) / REAL(6 * nkBZ, dp)
  beta(1:nbnd,1:nbnd,1:nkBZ) = beta(1:nbnd, 1:nbnd,1:nkBZ) / REAL(6 * nkBZ, dp)
  !
  WRITE(stdout,'(a,e15.7)') "[dfpt_tetra]  Dos(E_F)[/Ry] : ", SUM(dfde) * 2d0
  WRITE(stdout,'(a,e15.7)') "[dfpt_tetra]  # of electrons : ", SUM(occ) * 2d0
  !
  ! Map k-point 
  !
  IF(ALLOCATED(dfpt_tetra_dbl)) DEALLOCATE(dfpt_tetra_dbl) 
  IF(ALLOCATED(dfpt_tetra_occ    )) DEALLOCATE(dfpt_tetra_occ    )
  IF(ALLOCATED(dfpt_tetra_dlta   )) DEALLOCATE(dfpt_tetra_dlta   )
  !
  ALLOCATE(dfpt_tetra_dbl(nbnd,nbnd,nksq), &
  &        dfpt_tetra_occ(nbnd,nksq), dfpt_tetra_dlta(nbnd,nksq))
  !
  dfpt_tetra_dbl(1:nbnd,1:nbnd,1:nksq) = 0.0_dp
  dfpt_tetra_occ(       1:nbnd,1:nksq) = 0.0_dp
  dfpt_tetra_dlta(      1:nbnd,1:nksq) = 0.0_dp
  !  
  DO ik = 1, nksq
     !
     IF (lgamma) THEN
        ikk = ik
     ELSE
        ikk = 2 * ik - 1
     END IF
     !
     kv1(1:3) = MATMUL(xk(1:3,ikk), at(1:3,1:3))
     nksym = 0
     !
     DO isym = 1, nsymqbz
        !
        kv2(1:3) = MATMUL(REAL(symq(1:3,1:3,isym), dp), kv1(1:3)) * REAL((/nk1, nk2, nk3/), dp)
        ikv2(1:3) = NINT(kv2(1:3))
        !
        IF(ANY(ABS(kv2(1:3) - REAL(ikv2(1:3), dp)) > 1d-5)) CYCLE
        !
        ikv2(1:3) = MODULO(ikv2(1:3), (/nk1, nk2, nk3/)) + 1
        ik2 = indx(ikv2(1), ikv2(2), ikv2(3))
        !
        nksym = nksym + 1
        dfpt_tetra_dbl(1:nbnd,1:nbnd,ik) = dfpt_tetra_dbl( 1:nbnd,1:nbnd,ik) &
        &                                          + beta( 1:nbnd,1:nbnd,ik2)
        dfpt_tetra_occ(       1:nbnd,ik) = dfpt_tetra_occ( 1:nbnd,ik) &
        &                                +            occ( 1:nbnd,ik2) 
        dfpt_tetra_dlta(      1:nbnd,ik) = dfpt_tetra_dlta(1:nbnd,ik) &
        &                                           + dfde(1:nbnd,ik2) 
        !
     END DO
     !
     dfpt_tetra_dbl(1:nbnd,1:nbnd,ik) = dfpt_tetra_dbl(1:nbnd,1:nbnd,ik) / REAL(nksym, dp)
     dfpt_tetra_occ(       1:nbnd,ik) = dfpt_tetra_occ(       1:nbnd,ik) / REAL(nksym, dp)
     dfpt_tetra_dlta(      1:nbnd,ik) = dfpt_tetra_dlta(      1:nbnd,ik) / REAL(nksym, dp)
     !
  END DO
  !
  ! Compute wg & nbnd_occ
  ! 
  DO ik = 1, nks
     !
     kv1(1:3) = MATMUL(xk(1:3,ik), at(1:3,1:3))
     nksym = 0
     !
     wg(1:nbnd,ik) = 0d0
     !
     DO isym = 1, nsymbz
        !
        kv2(1:3) = MATMUL(REAL(sym(1:3,1:3,isym), dp), kv1(1:3)) * REAL((/nk1, nk2, nk3/), dp)
        ikv2(1:3) = NINT(kv2(1:3))
        !
        IF(ANY(ABS(kv2(1:3) - REAL(ikv2(1:3), dp)) > 1d-5)) CYCLE
        !
        ikv2(1:3) = MODULO(ikv2(1:3), (/nk1, nk2, nk3/)) + 1
        ik2 = indx(ikv2(1), ikv2(2), ikv2(3))
        !
        nksym = nksym + 1
        wg(1:nbnd,ik) = wg(1:nbnd,ik) + occ(1:nbnd,ik2) 
        !
     END DO
     !
     wg(1:nbnd,ik) = wg(1:nbnd,ik) / REAL(nksym, dp) * REAL(nkBZ, dp)
     !
     !nbnd_occ(ik) = 0
     !do ib = 1, nbnd
     !  IF(wg(ib,ik) > 1d-8) nbnd_occ(ik) = ib
     !END DO
     nbnd_occ(ik) = nbnd
     !
     wg(1:nbnd,ik) = wg(1:nbnd,ik) * wk(ik)
     !
  END DO ! ik = 1, nks
  !
  ! Integration weight(3rd term)
  !
  beta(1:nbnd,1:nbnd,1:nkBZ) = 0d0
  CALL dfpt_tetra_calc_beta3()
  WRITE(stdout,*) "[dfpt_tetra]  beta3 : ", SUM(beta)
  !
  beta(1:nbnd,1:nbnd,1:nkBZ) = beta(1:nbnd, 1:nbnd,1:nkBZ) / REAL(6 * nkBZ, dp)
  !
  kv1(1:3) = MATMUL(xq(1:3), at(1:3,1:3))
  kv1(1:3) = ABS(kv1(1:3) - REAL(NINT(kv1(1:3)), dp))
  !
  IF(lgamma .OR. MAXVAL(kv1(1:3)) < 1d-8) THEN
     !
     WRITE(stdout,*) "[dfpt_tetra]  Add Drude term"
     !
     DO ik = sttk, lstk
        DO ib = 1, nbnd
           beta(ib,ib,ik) = beta(ib, ib, ik) &
           &              - dfde(    ib, ik) * 0.5_dp
        END DO
     END DO
  END IF
  !
  CALL mp_sum(beta(1:nbnd,1:nbnd,1:nkBZ), world_comm)
  !
  ! Map k-point 
  !
  IF(ALLOCATED(dfpt_tetra_beta)) DEALLOCATE(dfpt_tetra_beta)
  ALLOCATE(dfpt_tetra_beta(nbnd,nbnd,nksq))
  !
  dfpt_tetra_beta(1:nbnd,1:nbnd,1:nksq) = 0d0
  !
  DO ik = 1, nksq
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     END IF
     !
     kv1(1:3) = MATMUL(xk(1:3,ikk), at(1:3,1:3))
     nksym = 0
     !
     DO isym = 1, nsymqbz
        !
        kv2(1:3) = MATMUL(REAL(symq(1:3,1:3,isym), dp), kv1(1:3)) * REAL((/nk1, nk2, nk3/), dp)
        ikv2(1:3) = NINT(kv2(1:3))
        !
        IF(ANY(ABS(kv2(1:3) - REAL(ikv2(1:3), dp)) > 1d-5)) CYCLE
        !
        ikv2(1:3) = MODULO(ikv2(1:3), (/nk1, nk2, nk3/)) + 1
        ik2 = indx(ikv2(1),ikv2(2),ikv2(3))
        !
        nksym = nksym + 1
        dfpt_tetra_beta(1:nbnd,1:nbnd,ik) = dfpt_tetra_beta(1:nbnd,1:nbnd,ik) &
        &                                            + beta(1:nbnd,1:nbnd,ik2) 
        !
     END DO ! isym = 1, nsymqbz
     !
     dfpt_tetra_beta(1:nbnd,1:nbnd,ik) = &
     &     dfpt_tetra_beta(1:nbnd,1:nbnd,ik) / REAL(nksym, dp)
     !
  END DO ! ik = 1, nksq
  !
  dfpt_tetra_beta(1:nbnd,1:nbnd,1:nksq) = dfpt_tetra_dbl(1:nbnd,1:nbnd,1:nksq) &
  &                                    + dfpt_tetra_beta(1:nbnd,1:nbnd,1:nksq) * alpha_pv
  !
  dfpt_tetra_occ(        1:nbnd,1:nksq) = dfpt_tetra_occ(        1:nbnd,1:nksq) * REAL(nkBZ, dp)
  dfpt_tetra_dlta(       1:nbnd,1:nksq) = dfpt_tetra_dlta(       1:nbnd,1:nksq) * REAL(nkBZ, dp)
  dfpt_tetra_beta(1:nbnd,1:nbnd,1:nksq) = dfpt_tetra_beta(1:nbnd,1:nbnd,1:nksq) * REAL(nkBZ, dp)
  dfpt_tetra_dbl( 1:nbnd,1:nbnd,1:nksq) = dfpt_tetra_dbl( 1:nbnd,1:nbnd,1:nksq) * REAL(nkBZ, dp)
  !
  DEALLOCATE(eig1,eig2,grid,indx,symq,sym)
  DEALLOCATE(occ,dfde,beta)
  !
END SUBROUTINE dfpt_tetra_main
!
!--------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_eig()
  !--------------------------------------------------------------------------
  !
  ! This routine collect eig_{k} & eig_{k + q} in whole of the Brillouin zone.
  !
  USE kinds,      ONLY : dp
  USE cell_base,  ONLY : at
  USE mp,         ONLY : mp_sum
  USE klist,      ONLY : xk
  USE wvfct,      ONLY : nbnd, et
  USE symm_base,  ONLY : s
  USE mp_global,   ONLY : inter_pool_comm
  USE start_k,    ONLY : nk1, nk2, nk3
  USE qpoint,     ONLY : nksq
  USE control_ph, ONLY : lgamma
  USE ener,       ONLY : ef
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ik2, ikk, ikq, isym, ikv2(3), i1, i2, i3
  REAL(dp) :: kv1(3), kv2(3)
  !
  nkBZ = nk1 * nk2 * nk3
  !
  ALLOCATE(eig1(nbnd, nkBZ), eig2(nbnd, nkBZ), &
  &        grid(3,nkBZ), indx(nk1,nk2,nk3))
  !
  ! K-grid
  !
  ik = 0
  DO i3 = 1, nk3
     DO i2 = 1, nk2
        DO i1 = 1, nk1
           !
           ik = ik + 1
           grid(1:3,ik) = (/ i1, i2, i3 /) - 1
           indx(i1, i2, i3) = ik
           !
        END DO
     END DO
  END DO
  !
  eig1(1:nbnd,1:nkBZ) = 0d0
  eig2(1:nbnd,1:nkBZ) = 0d0
  !
  DO ik = 1, nksq
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     END IF
     !
     kv1(1:3) = MATMUL(xk(1:3,ikk), at(1:3,1:3))
     !
     DO isym = 1, nsymqbz
        !
        kv2(1:3) = MATMUL(REAL(symq(1:3,1:3,isym), dp), kv1(1:3)) * REAL((/nk1, nk2, nk3/), dp)
        ikv2(1:3) = NINT(kv2(1:3))
        !
        IF(ANY(ABS(kv2(1:3) - REAL(ikv2(1:3), dp)) > 1d-5)) CYCLE
        !
        ikv2(1:3) = MODULO(ikv2(1:3), (/nk1, nk2, nk3/)) + 1
        ik2 = indx(ikv2(1),ikv2(2),ikv2(3))
        !
        eig1(1:nbnd,ik2) = et(1:nbnd,ikk) - ef
        eig2(1:nbnd,ik2) = et(1:nbnd,ikq) - ef
        !
     END DO
     !
  END DO
  !
  CALL mp_sum(eig1, inter_pool_comm )
  CALL mp_sum(eig2, inter_pool_comm )
  !
END SUBROUTINE dfpt_tetra_eig
!
!--------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_init()
  !--------------------------------------------------------------------------
  !
  ! This routine compute 20 points and weights for the optimized tetrahedron metod.
  !
  use kinds, only : dp
  USE cell_base, ONLY : bg
  USE io_global, ONLY : stdout
  USE start_k,    ONLY : nk1, nk2, nk3
  USE opt_tetra_mod, ONLY : tetra_type
  !
  IMPLICIT NONE
  !
  INTEGER :: itype, i1, i2, i3, it, &
  &          divvec(4,4), ivvec0(4)
  REAL(dp) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  bvec2(1:3,1) = bg(1:3,1) / REAL(nk1, dp)
  bvec2(1:3,2) = bg(1:3,2) / REAL(nk2, dp)
  bvec2(1:3,3) = bg(1:3,3) / REAL(nk3, dp)
  !
  bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
  bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  !
  ! length of delta bvec
  !
  DO i1 = 1, 4
     l(i1) = DOT_PRODUCT(bvec3(1:3,i1),bvec3(1:3,i1))
  END DO
  !
  itype = MINLOC(l(1:4),1)
  !
  ! start & last
  !
  ivvec0(1:4) = (/ 0, 0, 0, 0 /)
  !
  divvec(1:4,1) = (/ 1, 0, 0, 0 /)
  divvec(1:4,2) = (/ 0, 1, 0, 0 /)
  divvec(1:4,3) = (/ 0, 0, 1, 0 /)
  divvec(1:4,4) = (/ 0, 0, 0, 1 /)
  !
  ivvec0(itype) = 1
  divvec(itype, itype) = - 1
  !
  it = 0
  DO i1 = 1, 3
     DO i2 = 1, 3
        IF(i2 == i1) CYCLE
        DO i3 = 1, 3
           IF(i3 == i1 .OR. i3 == i2) CYCLE
           !
           it = it + 1
           !
           ivvec(1:3,1,it) = ivvec0(1:3)
           ivvec(1:3,2,it) = ivvec(1:3,1,it) + divvec(1:3,i1)
           ivvec(1:3,3,it) = ivvec(1:3,2,it) + divvec(1:3,i2)
           ivvec(1:3,4,it) = ivvec(1:3,3,it) + divvec(1:3,i3)
           !
        END DO
     END DO
  END DO
  !
  ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
  !
  ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,12,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
  !
  ivvec(1:3,13,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,14,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,15,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
  !
  ivvec(1:3,17,1:6) =  ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6)
  ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
  ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
  ivvec(1:3,20,1:6) =  ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6) + ivvec(1:3,1,1:6)
  !
  IF(tetra_type == 1) THEN
     !
     WRITE(stdout,*) "[dfpt_tetra]  Linear tetrahedron method is used."
     !
     wlsm(1:4,1:20) = 0.0_dp
     wlsm(1,1) = 1.0_dp
     wlsm(2,2) = 1.0_dp
     wlsm(3,3) = 1.0_dp
     wlsm(4,4) = 1.0_dp
     !
  ELSE IF(tetra_type == 2) THEN
     !
     WRITE(stdout,*) "[dfpt_tetra]  Optimized tetrahedron method is used."
     !
     wlsm(1, 1: 4) = REAL((/1440,    0,   30,    0/), dp)
     wlsm(2, 1: 4) = REAL((/   0, 1440,    0,   30/), dp)
     wlsm(3, 1: 4) = REAL((/  30,    0, 1440,    0/), dp)
     wlsm(4, 1: 4) = REAL((/   0,   30,    0, 1440/), dp)
     !
     wlsm(1, 5: 8) = REAL((/ -38,    7,   17,  -28/), dp)
     wlsm(2, 5: 8) = REAL((/ -28,  -38,    7,   17/), dp)
     wlsm(3, 5: 8) = REAL((/  17,  -28,  -38,    7/), dp)
     wlsm(4, 5: 8) = REAL((/   7,   17,  -28,  -38/), dp)
     !
     wlsm(1, 9:12) = REAL((/ -56,    9,  -46,    9/), dp)
     wlsm(2, 9:12) = REAL((/   9,  -56,    9,  -46/), dp)
     wlsm(3, 9:12) = REAL((/ -46,    9,  -56,    9/), dp)
     wlsm(4, 9:12) = REAL((/   9,  -46,    9,  -56/), dp)
     !
     wlsm(1,13:16) = REAL((/ -38,  -28,   17,    7/), dp)
     wlsm(2,13:16) = REAL((/   7,  -38,  -28,   17/), dp)
     wlsm(3,13:16) = REAL((/  17,    7,  -38,  -28/), dp)
     wlsm(4,13:16) = REAL((/ -28,   17,    7,  -38/), dp)
     !
     wlsm(1,17:20) = REAL((/ -18,  -18,   12,  -18/), dp)
     wlsm(2,17:20) = REAL((/ -18,  -18,  -18,   12/), dp)
     wlsm(3,17:20) = REAL((/  12,  -18,  -18,  -18/), dp)
     wlsm(4,17:20) = REAL((/ -18,   12,  -18,  -18/), dp)
     !
     wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260_dp
     !
  ELSE
     !
     CALL errore("dfpt_tetra_init", "tetra_type is invalid", tetra_type)
     !
  END IF
  !
END SUBROUTINE dfpt_tetra_init
!
!--------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_dlta_occ()
  !--------------------------------------------------------------------------
  !
  ! This routine compute occupation & its derivative
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE start_k, ONLY : nk1, nk2, nk3
  USE mp,         ONLY : mp_sum
  USE mp_global,   ONLY : world_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ikk, it, ib, jb, kb, ikv(3), ii
  REAL(dp) :: ei(nbnd,4), e(4), a(4,4), C(3), w0(4,4), tmp(5,4), &
  &           wdos1(4,nbnd), wocc1(4,nbnd), wdos(4), wocc(4), wdos2, wocc2
  !
  occ(1:nbnd,1:nkBZ) = 0.0_dp
  dfde(1:nbnd,1:nkBZ) = 0.0_dp
  !
  w0(1:4,1:4) = 0d0
  DO ii = 1, 4
     w0(ii,ii) = 1d0
  END DO
  !
  DO ik = sttk, lstk
     !
     DO it = 1, 6
        !
        ei(1:nbnd,1:4) = 0d0
        !
        DO ii = 1, 20
           !
           ikv(1:3) = grid(1:3,ik) + ivvec(1:3,ii,it)
           ikv(1:3) = MODULO(ikv(1:3),(/nk1, nk2, nk3/)) + 1
           ikk = indx(ikv(1),ikv(2),ikv(3))
           !
           DO ib = 1, nbnd
              ei(ib,1:4) = ei(ib,1:4) &
              &          + wlsm(1:4,ii) * eig1(ib,ikk)
           END DO
           !
        END DO
        !
        DO ib = 1, nbnd
           !
           tmp(1,1:4) = ei(ib,1:4)
           tmp(2:5,1:4) = w0(1:4, 1:4)
           !
           CALL dfpt_tetra_sort(5,4,tmp)
           !
           e(1:4) = tmp(1,1:4)
           !
           DO ii = 1, 4
              a(ii,1:4) = ( 0_dp - e(1:4) ) / (e(ii) - e(1:4))
           END DO
           !
           IF(e(1) <= 0_dp .AND. 0_dp < e(2)) THEN
              !
              ! DOS
              !
              C(1) = a(2,1) * a(3,1) * a(4,1) / (0_dp - e(1) )
              wdos(1) = a(1,2) + a(1,3) + a(1,4)
              wdos(2:4) = a(2:4,1)
              wdos(1:4) = wdos(1:4) * C(1)
              !
              ! OCCUPATION
              !
              C(1) = a(2,1) * a(3,1) *  a(4,1) / 4d0
              wocc(1) = C(1) * (1d0 + a(1,2) + a(1,3) + a(1,4))
              wocc(2) = C(1) * a(2,1)
              wocc(3) = C(1) * a(3,1)
              wocc(4) = C(1) * a(4,1)
              !
           ELSE IF(e(2) <= 0_dp .AND. 0_dp < e(3)) THEN
              !
              ! DOS
              !
              C(1) = a(2,3) * a(3,1) + a(3,2) * a(2,4)
              wdos(1) = a(1,4) * C(1) + a(1,3) * a(3,1) * a(2,3)
              wdos(2) = a(2,3) * C(1) + a(2,4) * a(2,4) * a(3,2)
              wdos(3) = a(3,2) * C(1) + a(3,1) * a(3,1) * a(2,3)
              wdos(4) = a(4,1) * C(1) + a(4,2) * a(2,4) * a(3,2)
              C(1)  = 1d0 / ( e(4) - e(1) )
              wdos(1:4) = wdos(1:4) * C(1)
              !
              ! OCCUPATION
              !
              C(1) = a(4,1) * a(3,1) / 4d0
              C(2) = a(4,1) * a(3,2) * a(1,3) / 4d0
              C(3) = a(4,2) * a(3,2) * a(1,4) / 4d0
              !
              wocc(1) = C(1) + (C(1) + C(2)) * a(1,3) + (C(1) + C(2) + C(3)) * a(1,4)
              wocc(2) = C(1) + C(2) + C(3) + (C(2) + C(3)) * a(2,3) + C(3) * a(2,4)
              wocc(3) = (C(1) + C(2)) * a(3,1) + (C(2) + C(3)) * a(3,2)
              wocc(4) = (C(1) + C(2) + C(3)) * a(4,1) + C(3) * a(4,2)
              !
           ELSE IF(e(3) <= 0_dp .AND. 0_dp < e(4)) THEN
              !
              ! DOS
              !
              C(1) = a(1,4) * a(2,4) * a(3,4) / ( e(4) - 0_dp )
              wdos(1:3)  = a(1:3,4)
              wdos(4)  = a(4,1) + a(4,2) + a(4,3)
              wdos(1:4) = wdos(1:4) * C(1)
              !
              ! OCCUPATION
              !
              C(1) = a(1,4) * a(2,4) * a(3,4)
              !
              wocc(1) = 1d0 - C(1) * a(1,4)
              wocc(2) = 1d0 - C(1) * a(2,4)
              wocc(3) = 1d0 - C(1) * a(3,4)
              wocc(4) = 1d0 - C(1) * (1d0 + a(4,1) + a(4,2) + a(4,3))
              !
              wocc(1:4) = wocc(1:4) / 4d0
              !
           ELSE IF(e(4) <= 0_dp) THEN
              !
              wdos(1:4) = 0d0
              wocc(1:4) = 1d0 / 4d0
              !
           ELSE
              !
              wdos(1:4) = 0d0
              wocc(1:4) = 0d0
              !
           END IF
           !
           wdos1(1:4, ib) = MATMUL(tmp(2:5, 1:4), wdos(1:4))
           wocc1(1:4, ib) = MATMUL(tmp(2:5, 1:4), wocc(1:4))
           !
        END DO ! ib
        !
        DO ii = 1, 20
           !
           ikv(1:3) = grid(1:3,ik) + ivvec(1:3,ii,it)
           ikv(1:3) = MODULO(ikv(1:3), (/nk1, nk2, nk3/)) + 1
           ikk = indx(ikv(1),ikv(2),ikv(3))
           !
           dfde(1:nbnd,ikk) = dfde(1:nbnd,ikk) &
           &            + MATMUL(wlsm(1:4,ii), wdos1(1:4,1:nbnd))
           occ( 1:nbnd,ikk) = occ( 1:nbnd,ikk) &
           &            + MATMUL(wlsm(1:4,ii), wocc1(1:4,1:nbnd))
           !
        END DO
        !
     END DO ! it
     !
  END DO ! ik
  !
  ! Average weights of degenerated states
  !
  DO ik = sttk, lstk
     DO ib = 1, nbnd
        !
        wdos2 = dfde(ib,ik)
        wocc2 = occ( ib,ik)
        !
        DO jb = ib + 1, nbnd
           !
           IF(ABS(eig1(ib,ik) - eig1(jb,ik)) < 1e-6_dp) THEN
              wdos2 = wdos2 + dfde(jb,ik)
              wocc2 = wocc2 +  occ(jb,ik)
           ELSE
              !
              DO kb = ib, jb - 1
                 dfde(kb,ik) = wdos2 / real(jb - ib, dp)
                 occ( kb,ik) = wocc2 / real(jb - ib, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
     END DO
  END DO
  !
  CALL mp_sum(occ(1:nbnd,1:nkBZ), world_comm)
  CALL mp_sum(dfde(1:nbnd,1:nkBZ), world_comm)
  !
END SUBROUTINE dfpt_tetra_dlta_occ
!
!--------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_calc_beta1()
  !--------------------------------------------------------------------------
  !
  ! This routine compute the first term of (B 28) in PRB 64, 235118 (2001).
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE start_k, ONLY : nk1, nk2, nk3
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, it, ib, jb, kb, ikv(3), ii, ikk
  REAL(dp) :: thr = 1d-8, V, e(4), a(4,4)
  REAL(dp) :: ei(nbnd,4), ej(nbnd,4), ei2(4), ej2(nbnd,4)
  REAL(dp) :: tmp(6,nbnd,4), tmp2(6,nbnd,4)
  REAL(dp) :: w1(4,nbnd), w2(4,nbnd,4), beta2(nbnd)
  REAL(dp),ALLOCATABLE :: w0(:,:,:)
  !
  ALLOCATE(w0(4,nbnd,4))
  !
  w0(1:4,1:nbnd,1:4) = 0d0
  DO ii = 1, 4
     w0(ii,1:nbnd,ii) = 1d0
  END DO
  tmp(1:6,1:nbnd,1:4) = 0d0
  !
  DO ik = sttk, lstk
     !
     DO it = 1, 6
        !
        ei(1:nbnd, 1:4) = 0d0
        ej(1:nbnd, 1:4) = 0d0
        !
        DO ii = 1, 20
           !
           ikv(1:3) = grid(1:3,ik) + ivvec(1:3,ii,it)
           ikv(1:3) = MODULO(ikv(1:3),(/nk1, nk2, nk3/)) + 1
           ikk = indx(ikv(1),ikv(2),ikv(3))
           !
           DO ib = 1, nbnd
              !
              ei(ib, 1:4) = ei(ib, 1:4) + wlsm(1:4,ii) * eig1(ib, ikk)
              ej(ib, 1:4) = ej(ib, 1:4) + wlsm(1:4,ii) * eig2(ib, ikk)
              !               
           END DO
           !
        END DO
        !
        DO ib = 1, nbnd
           !
           w1(1:4,1:nbnd) = 0d0
           !
           tmp(1,        1, 1:4) = ei(         ib, 1:4)
           tmp(2,   1:nbnd, 1:4) = ej(     1:nbnd, 1:4)
           tmp(3:6, 1:nbnd, 1:4) = w0(1:4, 1:nbnd, 1:4)
           !
           CALL dfpt_tetra_sort(6 * nbnd, 4, tmp)
           !
           e(1:4) = tmp(1, 1, 1:4)
           !
           DO ii = 1, 4
              a(ii,1:4) = ( 0d0 - e(1:4)) / (e(ii) - e(1:4))
           END DO
           !
           IF(e(1) <= 0d0 .AND. 0d0 < e(2)) THEN
              !
              ! A - 1
              !
              V = a(2,1) * a(3,1) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,2) &
                 &                  + tmp(1:6,1:nbnd,2) * a(2,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,1) * a(1,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,1)
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,1) * a(1,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,1)
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + v * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
           ELSE IF( e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
              !
              ! B - 1
              !
              V = a(3,1) * a(4,1) * a(2,4)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,1) * a(1,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,1) 
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! B - 2
              !
              V = a(3,2) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1:2) = tmp(1:6,1:nbnd,1:2)
                 tmp2(1:6,1:nbnd,3)   = tmp(1:6,1:nbnd,2) * a(2,3) &
                 &                    + tmp(1:6,1:nbnd,3) * a(3,2) 
                 tmp2(1:6,1:nbnd,4)   = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! B - 3
              !
              V = a(2,3) * a(3,1) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,2) * a(2,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,2) 
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
           ELSE IF( e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
              !
              ! C - 1
              !
              V = a(4,3)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1:3) = tmp(1:6,1:nbnd,1:3)
                 tmp2(1:6,1:nbnd,4)   = tmp(1:6,1:nbnd,3) * a(3,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! C - 2
              !
              V = a(3,4) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1:2) = tmp(1:6,1:nbnd,1:2)
                 tmp2(1:6,1:nbnd,3)   = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,2) 
                 tmp2(1:6,1:nbnd,4)   = tmp(1:6,1:nbnd,3) * a(3,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! C - 3
              !
              V = a(3,4) * a(2,4) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,2) 
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,3) * a(3,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
           ELSE IF( e(4) <= 0d0 ) THEN
              !
              ! D - 1
              !
              V = 1d0
              !             
              tmp2(1:6,1:nbnd,1:4) = tmp(1:6,1:nbnd,1:4)
              !
              ei2(            1:4) = tmp2(  1,      1, 1:4)
              ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
              w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
              ! 
              CALL dfpt_tetra2_theta(ei2,ej2,w2)
              !
              w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
              &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
              !
           END IF
           !
           DO ii = 1, 20
              !
              ikv(1:3) = grid(1:3,ik) + ivvec(1:3,ii,it)
              ikv(1:3) = MODULO(ikv(1:3), (/nk1, nk2, nk3/)) + 1
              ikk = indx(ikv(1),ikv(2),ikv(3))
              !
              beta(1:nbnd,ib,ikk) = beta(1:nbnd,ib,ikk) &
              &               + MATMUL(wlsm(1:4,ii), w1(1:4,1:nbnd))
              !               
           END DO
           !
        END DO ! ib
        !
     END DO ! it
     !
  END DO ! ik
  !
  ! Average weights of degenerated states
  !
  DO ik = sttk, lstk
     DO ib = 1, nbnd
        !
        beta2(1:nbnd) = beta(1:nbnd,ib,ik)
        !
        DO jb = ib + 1, nbnd
           !
           IF(ABS(eig1(ib,ik) - eig1(jb,ik)) < 1e-6_dp) THEN
              beta2(1:nbnd) = beta2(1:nbnd) + beta(1:nbnd,jb,ik)
           ELSE
              !
              DO kb = ib, jb - 1
                 beta(1:nbnd,kb,ik) = beta2(1:nbnd) / real(jb - ib, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
        beta2(1:nbnd) = beta(ib,1:nbnd,ik)
        !
        DO jb = ib + 1, nbnd
           !
           IF(ABS(eig1(ib,ik) - eig1(jb,ik)) < 1e-6_dp) THEN
              beta2(1:nbnd) = beta2(1:nbnd) + beta(jb,1:nbnd,ik)
           ELSE
              !
              DO kb = ib, jb - 1
                 beta(kb,1:nbnd,ik) = beta2(1:nbnd) / real(jb - ib, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
     END DO
  END DO
  !
  DEALLOCATE(w0)
  !
END SUBROUTINE dfpt_tetra_calc_beta1
!
!----------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_sort(n1,n2,a)
  !----------------------------------------------------------------------------
  !
  ! Simple sort
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n1, n2
  REAL(dp),INTENT(INOUT) :: a(n1,n2) 
  !
  INTEGER :: i, m
  REAL(dp) :: am, atmp(n1)
  !
  DO i = 1, n2 - 1
     am = MINVAL(a(1,i+1:n2) )
     m  = MINLOC(a(1,i+1:n2),1) + i
     IF(a(1,i) > am) THEN
        atmp(1:n1) = a(1:n1,m)
        a(1:n1,m) = a(1:n1,i)
        a(1:n1,i) = atmp(1:n1)
     END IF
  END DO
  !
END SUBROUTINE dfpt_tetra_sort
!
!--------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_calc_beta2()
  !--------------------------------------------------------------------------
  !
  ! This routine compute the second term of (B 28) in PRB 64, 235118 (2001).
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE start_k, ONLY : nk1, nk2, nk3
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ikk, it, ib, jb, kb, ikv(3), ii
  REAL(dp) :: thr = 1d-8, V, e(4), a(4,4)
  REAL(dp) :: ei(nbnd,4), ej(nbnd,4), ei2(4), ej2(nbnd,4)
  REAL(dp) :: tmp(6,nbnd,4), tmp2(6,nbnd,4)
  REAL(dp) :: w1(4,nbnd), w2(4,nbnd,4), beta2(nbnd)
  REAL(dp),ALLOCATABLE :: w0(:,:,:)
  !
  ALLOCATE(w0(4,nbnd,4))
  !
  w0(1:4,1:nbnd,1:4) = 0d0
  DO ii = 1, 4
     w0(ii,1:nbnd,ii) = 1d0
  END DO
  tmp(1:6,1:nbnd,1:4) = 0d0
  !
  DO ik = sttk, lstk
     !
     DO it = 1, 6
        !
        ei(1:nbnd, 1:4) = 0d0
        ej(1:nbnd, 1:4) = 0d0
        !
        DO ii = 1, 20
           !
           ikv(1:3) = grid(1:3,ik) + ivvec(1:3,ii,it)
           ikv(1:3) = MODULO(ikv(1:3),(/nk1, nk2, nk3/)) + 1
           ikk = indx(ikv(1),ikv(2),ikv(3))
           !
           DO ib = 1, nbnd
              !
              ei(ib, 1:4) = ei(ib, 1:4) + wlsm(1:4,ii) * eig2(ib, ikk)
              ej(ib, 1:4) = ej(ib, 1:4) + wlsm(1:4,ii) * eig1(ib, ikk)
              !
           END DO
           !
        END DO
        !
        DO ib = 1, nbnd
           !
           w1(1:4,1:nbnd) = 0d0
           !
           tmp(1,        1, 1:4) = ei(         ib, 1:4)
           tmp(2,   1:nbnd, 1:4) = ej(     1:nbnd, 1:4)
           tmp(3:6, 1:nbnd, 1:4) = w0(1:4, 1:nbnd, 1:4)
           !
           CALL dfpt_tetra_sort(6 * nbnd, 4, tmp)
           !
           e(1:4) = tmp(1, 1, 1:4)
           !
           DO ii = 1, 4
              a(ii,1:4) = ( 0d0 - e(1:4)) / (e(ii) - e(1:4))
           END DO
           !
           IF(e(1) <= 0d0 .AND. 0d0 < e(2)) THEN
              !
              ! A - 1
              !
              V = a(2,1) * a(3,1) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,2) &
                 &                  + tmp(1:6,1:nbnd,2) * a(2,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,1) * a(1,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,1)
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,1) * a(1,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,1)
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
           ELSE IF( e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
              !
              ! B - 1
              !
              V = a(3,1) * a(4,1) * a(2,4)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,1) * a(1,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,1) 
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! B - 2
              !
              V = a(3,2) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1:2) = tmp(1:6,1:nbnd,1:2)
                 tmp2(1:6,1:nbnd,3)   = tmp(1:6,1:nbnd,2) * a(2,3) &
                 &                    + tmp(1:6,1:nbnd,3) * a(3,2) 
                 tmp2(1:6,1:nbnd,4)   = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! B - 3
              !
              V = a(2,3) * a(3,1) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,2) * a(2,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,2) 
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
           ELSE IF( e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
              !
              ! C - 1
              !
              V = a(4,3)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1:3) = tmp(1:6,1:nbnd,1:3)
                 tmp2(1:6,1:nbnd,4)   = tmp(1:6,1:nbnd,3) * a(3,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! C - 2
              !
              V = a(3,4) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1:2) = tmp(1:6,1:nbnd,1:2)
                 tmp2(1:6,1:nbnd,3)   = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,2) 
                 tmp2(1:6,1:nbnd,4)   = tmp(1:6,1:nbnd,3) * a(3,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! C - 3
              !
              V = a(3,4) * a(2,4) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,2) 
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,3) * a(3,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_theta(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
           ELSE IF( e(4) <= 0d0 ) THEN
              !
              ! D - 1
              !
              V = 1d0
              !             
              tmp2(1:6,1:nbnd,1:4) = tmp(1:6,1:nbnd,1:4)
              !
              ei2(            1:4) = tmp2(  1,      1, 1:4)
              ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
              w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
              ! 
              CALL dfpt_tetra2_theta(ei2,ej2,w2)
              !
              w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
              &      + V * SUM(w2(1:4,1:nbnd,1:4), 3)
              !
           END IF
           !
           DO ii = 1, 20
              !
              ikv(1:3) = grid(1:3,ik) + ivvec(1:3,ii,it)
              ikv(1:3) = MODULO(ikv(1:3), (/nk1, nk2, nk3/)) + 1
              ikk = indx(ikv(1),ikv(2),ikv(3))
              !
              beta(ib,1:nbnd,ikk) = beta(ib,1:nbnd,ikk) &
              &                + MATMUL(wlsm(1:4,ii), w1(1:4,1:nbnd))
              !               
           END DO
           !
        END DO ! ib
        !
     END DO ! it
     !
  END DO ! ik
  !
  ! Average weights of degenerated states
  !
  DO ik = sttk, lstk
     DO ib = 1, nbnd
        !
        beta2(1:nbnd) = beta(1:nbnd,ib,ik)
        !
        DO jb = ib + 1, nbnd
           !
           IF(ABS(eig1(ib,ik) - eig1(jb,ik)) < 1e-6_dp) THEN
              beta2(1:nbnd) = beta2(1:nbnd) + beta(1:nbnd,jb,ik)
           ELSE
              !
              DO kb = ib, jb - 1
                 beta(1:nbnd,kb,ik) = beta2(1:nbnd) / real(jb - ib, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
        beta2(1:nbnd) = beta(ib,1:nbnd,ik)
        !
        DO jb = ib + 1, nbnd
           !
           IF(ABS(eig1(ib,ik) - eig1(jb,ik)) < 1e-6_dp) THEN
              beta2(1:nbnd) = beta2(1:nbnd) + beta(jb,1:nbnd,ik)
           ELSE
              !
              DO kb = ib, jb - 1
                 beta(kb,1:nbnd,ik) = beta2(1:nbnd) / real(jb - ib, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
     END DO
  END DO
  !
  DEALLOCATE(w0)
  !
END SUBROUTINE dfpt_tetra_calc_beta2
!
!----------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_calc_beta3()
  !--------------------------------------------------------------------------
  !
  ! This routine compute the third term of (B 28) in PRB 64, 235118 (2001).
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE start_k, ONLY : nk1, nk2, nk3
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, it, ib, jb, kb, ikv(3), ii, ikk
  REAL(dp) :: thr = 1d-8, V, e(4), a(4,4)
  REAL(dp) :: ei(nbnd,4), ej(nbnd,4), ei2(4), ej2(nbnd,4)
  REAL(dp) :: tmp(6,nbnd,4), tmp2(6,nbnd,4)
  REAL(dp) :: w1(4,nbnd), w2(4,nbnd,4), beta2(nbnd)
  REAL(dp),ALLOCATABLE :: w0(:,:,:)
  !
  ALLOCATE(w0(4,nbnd,4))
  !
  w0(1:4,1:nbnd,1:4) = 0d0
  DO ii = 1, 4
     w0(ii,1:nbnd,ii) = 1d0
  END DO
  tmp(1:6,1:nbnd,1:4) = 0d0
  !
  DO ik = sttk, lstk
     !
     DO it = 1, 6
        !
        ei(1:nbnd, 1:4) = 0d0
        ej(1:nbnd, 1:4) = 0d0
        !
        DO ii = 1, 20
           !
           ikv(1:3) = grid(1:3,ik) + ivvec(1:3,ii,it)
           ikv(1:3) = MODULO(ikv(1:3),(/nk1, nk2, nk3/)) + 1
           ikk = indx(ikv(1),ikv(2),ikv(3))
           !
           DO ib = 1, nbnd
              !
              ei(ib, 1:4) = ei(ib, 1:4) + wlsm(1:4,ii) * eig1(ib, ikk)
              ej(ib, 1:4) = ej(ib, 1:4) + wlsm(1:4,ii) * eig2(ib, ikk)
              !               
           END DO
           !
        END DO
        !
        DO ib = 1, nbnd
           !
           w1(1:4,1:nbnd) = 0d0
           !
           tmp(1,        1, 1:4) = ei(         ib, 1:4)
           tmp(2,   1:nbnd, 1:4) = ej(     1:nbnd, 1:4)
           tmp(3:6, 1:nbnd, 1:4) = w0(1:4, 1:nbnd, 1:4)
           !
           CALL dfpt_tetra_sort(6 * nbnd, 4, tmp)
           !
           e(1:4) = tmp(1, 1, 1:4)
           !
           DO ii = 1, 4
              a(ii,1:4) = (0d0 - e(1:4)) / (e(ii) - e(1:4))
           END DO
           !
           IF(e(1) <= 0d0 .AND. 0d0 < e(2)) THEN
              !
              ! A - 1
              !
              V = a(2,1) * a(3,1) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,2) &
                 &                  + tmp(1:6,1:nbnd,2) * a(2,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,1) * a(1,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,1)
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,1) * a(1,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,1)
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_lindhard(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      - V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
           ELSE IF( e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
              !
              ! B - 1
              !
              V = a(3,1) * a(4,1) * a(2,4)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,1) * a(1,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,1) 
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_lindhard(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      - V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! B - 2
              !
              V = a(3,2) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1:2) = tmp(1:6,1:nbnd,1:2)
                 tmp2(1:6,1:nbnd,3)   = tmp(1:6,1:nbnd,2) * a(2,3) &
                 &                    + tmp(1:6,1:nbnd,3) * a(3,2) 
                 tmp2(1:6,1:nbnd,4)   = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_lindhard(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      - V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! B - 3
              !
              V = a(2,3) * a(3,1) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,2) * a(2,3) &
                 &                  + tmp(1:6,1:nbnd,3) * a(3,2) 
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,2) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_lindhard(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      - V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
           ELSE IF( e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
              !
              ! C - 1
              !
              V = a(4,3)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1:3) = tmp(1:6,1:nbnd,1:3)
                 tmp2(1:6,1:nbnd,4)   = tmp(1:6,1:nbnd,3) * a(3,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_lindhard(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      - V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! C - 2
              !
              V = a(3,4) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1:2) = tmp(1:6,1:nbnd,1:2)
                 tmp2(1:6,1:nbnd,3)   = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,2) 
                 tmp2(1:6,1:nbnd,4)   = tmp(1:6,1:nbnd,3) * a(3,4) &
                 &                    + tmp(1:6,1:nbnd,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_lindhard(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      - V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
              ! C - 3
              !
              V = a(3,4) * a(2,4) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tmp2(1:6,1:nbnd,1) = tmp(1:6,1:nbnd,1)
                 tmp2(1:6,1:nbnd,2) = tmp(1:6,1:nbnd,1) * a(1,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,1) 
                 tmp2(1:6,1:nbnd,3) = tmp(1:6,1:nbnd,2) * a(2,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,2) 
                 tmp2(1:6,1:nbnd,4) = tmp(1:6,1:nbnd,3) * a(3,4) &
                 &                  + tmp(1:6,1:nbnd,4) * a(4,3) 
                 !
                 ei2(            1:4) = tmp2(  1,      1, 1:4)
                 ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
                 w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
                 ! 
                 CALL dfpt_tetra2_lindhard(ei2,ej2,w2)
                 !
                 w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
                 &      - V * SUM(w2(1:4,1:nbnd,1:4), 3)
                 !
              END IF
              !
           ELSE IF( e(4) <= 0d0 ) THEN
              !
              ! D - 1
              !
              V = 1d0
              !             
              tmp2(1:6,1:nbnd,1:4) = tmp(1:6,1:nbnd,1:4)
              !
              ei2(            1:4) = tmp2(  1,      1, 1:4)
              ej2(    1:nbnd, 1:4) = tmp2(  2, 1:nbnd, 1:4)
              w2(1:4, 1:nbnd, 1:4) = tmp2(3:6, 1:nbnd, 1:4)
              ! 
              CALL dfpt_tetra2_lindhard(ei2,ej2,w2)
              !
              w1(1:4,1:nbnd) = w1(1:4,1:nbnd) &
              &      - V * SUM(w2(1:4,1:nbnd,1:4), 3)
              !
           END IF
           !
           DO ii = 1, 20
              !
              ikv(1:3) = grid(1:3,ik) + ivvec(1:3,ii,it)
              ikv(1:3) = MODULO(ikv(1:3), (/nk1, nk2, nk3/)) + 1
              ikk = indx(ikv(1),ikv(2),ikv(3))
              !
              beta(1:nbnd,ib,ikk) = beta(1:nbnd,ib,ikk) &
              &                   + MATMUL(wlsm(1:4,ii), w1(1:4,1:nbnd))
              !               
           END DO
           !
        END DO ! ib
        !
     END DO ! it
     !
  END DO ! ik
  !
  ! Average weights of degenerated states
  !
  DO ik = sttk, lstk
     DO ib = 1, nbnd
        !
        beta2(1:nbnd) = beta(1:nbnd,ib,ik)
        !
        DO jb = ib + 1, nbnd
           !
           IF(ABS(eig1(ib,ik) - eig1(jb,ik)) < 1e-6_dp) THEN
              beta2(1:nbnd) = beta2(1:nbnd) + beta(1:nbnd,jb,ik)
           ELSE
              !
              DO kb = ib, jb - 1
                 beta(1:nbnd,kb,ik) = beta2(1:nbnd) / real(jb - ib, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
        beta2(1:nbnd) = beta(ib,1:nbnd,ik)
        !
        DO jb = ib + 1, nbnd
           !
           IF(ABS(eig1(ib,ik) - eig1(jb,ik)) < 1e-6_dp) THEN
              beta2(1:nbnd) = beta2(1:nbnd) + beta(jb,1:nbnd,ik)
           ELSE
              !
              DO kb = ib, jb - 1
                 beta(kb,1:nbnd,ik) = beta2(1:nbnd) / real(jb - ib, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
     END DO
  END DO
  !
  DEALLOCATE(w0)
  !
END SUBROUTINE dfpt_tetra_calc_beta3
!
!-----------------------------------------------------------
SUBROUTINE dfpt_tetra2_theta(ei,ej,w)
  !
  ! This routine compute theta(ei - ej)
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: ei(4), ej(nbnd,4)
  REAL(dp),INTENT(INOUT) :: w(4,nbnd,4)
  !
  INTEGER :: ii, ib
  REAL(dp) :: V, w2(4,4), thr = 1d-8
  REAL(dp) :: tmp(5,4), e(4), a(4,4)
  !
  DO ib = 1, nbnd
     !
     tmp(  1,1:4) = ej(ib,1:4) - ei(1:4)
     tmp(2:5,1:4) = w(1:4,ib,1:4)
     CALL dfpt_tetra_sort(5, 4, tmp)
     e(       1:4) = tmp(1,1:4)
     w(1:4,ib,1:4) = 0d0
     !
     DO ii = 1, 4
        a(ii,1:4) = (0d0 - e(1:4)) / (e(ii) - e(1:4))
     END DO
     !
     IF(ABS(e(1)) < thr .AND. ABS(e(4)) < thr) THEN
        !
        ! Theta(0) = 0.5
        !
        V = 0.5_dp * 0.25_dp
        !
        w2(1:4,1:4) = tmp(2:5, 1:4)
        !
        w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
        !
     ELSE IF((e(1) <= 0d0 .AND. 0d0 < e(2)) .OR. (e(1) < 0d0 .AND. 0d0 <= e(2))) THEN
        !
        ! A - 1
        !
        V = 0.25d0 * a(2,1) * a(3,1) * a(4,1)
        !
        IF(V > thr) THEN
           !
           w2(1:4,1) = tmp(2:5,1)
           w2(1:4,2) = tmp(2:5,1) * a(1,2) + tmp(2:5,2) * a(2,1) 
           w2(1:4,3) = tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1)
           w2(1:4,4) = tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1)
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        END IF
        !
     ELSE IF((e(2) <= 0d0 .AND. 0d0 < e(3)) .OR. (e(2) < 0d0 .AND. 0d0 <= e(3))) THEN
        !
        ! B - 1
        !
        V = 0.25d0 * a(3,1) * a(4,1) * a(2,4)
        !
        IF(V > thr) THEN
           !
           w2(1:4,1) = tmp(2:5,1)
           w2(1:4,2) = tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) 
           w2(1:4,3) = tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) 
           w2(1:4,4) = tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        END IF
        !
        ! B - 2
        !
        V = 0.25d0 * a(3,2) * a(4,2)
        !
        IF(V > thr) THEN
           !
           w2(1:4,1:2) = tmp(2:5,1:2)
           w2(1:4,3)   = tmp(2:5,2) * a(2,3) + tmp(2:5,3) * a(3,2) 
           w2(1:4,4)   = tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        END IF
        !
        ! B - 3
        !
        V = 0.25d0 * a(2,3) * a(3,1) * a(4,2)
        !
        IF(V > thr) THEN
           !
           w2(1:4,1) = tmp(2:5,1)
           w2(1:4,2) = tmp(2:5,1) * a(1,3) + tmp(2:5,3) * a(3,1) 
           w2(1:4,3) = tmp(2:5,2) * a(2,3) + tmp(2:5,3) * a(3,2) 
           w2(1:4,4) = tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        END IF
        !
     ELSE IF((e(3) <= 0d0 .AND. 0d0 < e(4)) .OR. (e(3) < 0d0 .AND. 0d0 <= e(4))) THEN
        !
        ! C - 1
        !
        V = 0.25d0 * a(4,3)
        !
        IF(V > thr) THEN
           !
           w2(1:4,1:3) = tmp(2:5,1:3)
           w2(1:4,4)   = tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        END IF
        !
        ! C - 2
        !
        V = 0.25d0 * a(3,4) * a(4,2)
        !
        IF(V > thr) THEN
           !
           w2(1:4,1:2) = tmp(2:5,1:2)
           w2(1:4,3)   = tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) 
           w2(1:4,4)   = tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        END IF
        !
        ! C - 3
        !
        V = 0.25d0 * a(3,4) * a(2,4) * a(4,1)
        !
        IF(V > thr) THEN
           !
           w2(1:4,1) = tmp(2:5,1)
           w2(1:4,2) = tmp(2:5,1) * a(1,4) + tmp(2:5,4) * a(4,1) 
           w2(1:4,3) = tmp(2:5,2) * a(2,4) + tmp(2:5,4) * a(4,2) 
           w2(1:4,4) = tmp(2:5,3) * a(3,4) + tmp(2:5,4) * a(4,3) 
           !
           w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
           !
        END IF
        !
     ELSE IF(e(4) <= 0d0) THEN
        !
        ! D - 1
        !
        V = 0.25d0
        !             
        w2(1:4,1:4) = tmp(2:5,1:4)
        !
        w(1:4, ib,1:4) = w(1:4, ib, 1:4) + w2(1:4, 1:4) * V
        !
     END IF
     !
  END DO
  !
END SUBROUTINE dfpt_tetra2_theta
!
!-----------------------------------------------------------------------
SUBROUTINE dfpt_tetra2_lindhard(ei,ej,w)
  !---------------------------------------------------------------------
  !
  ! This routine take the unoccupied region.
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: ei(4), ej(nbnd,4)
  REAL(dp),INTENT(INOUT) :: w(4,nbnd,4)
  !
  INTEGER :: ii, ib
  REAL(dp) :: V, ei2(4), ej2(4), w2(4,4), thr = 1d-8
  REAL(dp) :: tmp(6,4), tmp2(6,4), e(4), a(4,4)
  !
  DO ib = 1, nbnd
     !
     tmp(  1, 1:4) = ej(ib,   1:4)
     tmp(  2, 1:4) = ei(      1:4)
     tmp(3:6, 1:4) = w(1:4,ib,1:4)
     CALL dfpt_tetra_sort(6, 4, tmp)
     e(       1:4) = tmp(1,1:4)
     w(1:4,ib,1:4) = 0d0
     !
     DO ii = 1, 4
        a(ii,1:4) = ( 0d0 - e(1:4) ) / (e(ii) - e(1:4))
     END DO
     !
     IF(0d0 <= e(1)) THEN
        !
        ! A - 1
        !
        V = 1d0
        !
        tmp2(1:6,1:4) = tmp(1:6,1:4)
        !
        ej2(   1:4) = tmp2(  1,1:4)
        ei2(   1:4) = tmp2(  2,1:4)
        w2(1:4,1:4) = tmp2(3:6,1:4)
        !
        CALL dfpt_tetra_lindhard(ei2,ej2,w2)
        w(1:4,ib,1:4) = w(1:4,ib,1:4) + w2(1:4,1:4) * V
        !
     ELSE IF((e(1) < 0d0 .AND. 0d0 <= e(2)) .OR. (e(1) <= 0d0 .AND. 0d0 < e(2))) THEN
        !
        ! B - 1
        !
        V = a(1,2)
        !
        IF(V > thr) THEN
           !
           tmp2(1:6,1)   = tmp(1:6,1) * a(1,2) + tmp(1:6,2) * a(2,1)
           tmp2(1:6,2:4) = tmp(1:6,2:4)
           !
           ej2(   1:4) = tmp2(1,  1:4)
           ei2(   1:4) = tmp2(2,  1:4)
           w2(1:4,1:4) = tmp2(3:6,1:4)
           !
           CALL dfpt_tetra_lindhard(ei2,ej2,w2)
           w(1:4,ib,1:4) = w(1:4,ib,1:4) + w2(1:4,1:4) * V
           !       
        END IF
        !
        ! B - 2
        !
        V = a(1,3) * a(2,1)
        !
        IF(V > thr) THEN
           !
           tmp2(1:6,1) = tmp(1:6,1) * a(1,2) + tmp(1:6,2) * a(2,1)
           tmp2(1:6,2) = tmp(1:6,1) * a(1,3) + tmp(1:6,3) * a(3,1)
           tmp2(1:6,3:4) = tmp(1:6,3:4)
           !
           ej2(   1:4) = tmp2(  1,1:4)
           ei2(   1:4) = tmp2(  2,1:4)
           w2(1:4,1:4) = tmp2(3:6,1:4)
           !
           CALL dfpt_tetra_lindhard(ei2,ej2,w2)
           w(1:4,ib,1:4) = w(1:4,ib,1:4) + w2(1:4,1:4) * V
           !         
        END IF
        !
        ! B - 3
        !
        V = a(1,4) * a(2,1) * a(3,1)
        !
        IF(V > thr) THEN
           !
           tmp2(1:6,1) = tmp(1:6,1) * a(1,2) + tmp(1:6,2) * a(2,1)
           tmp2(1:6,2) = tmp(1:6,1) * a(1,3) + tmp(1:6,3) * a(3,1)
           tmp2(1:6,3) = tmp(1:6,1) * a(1,4) + tmp(1:6,4) * a(4,1)
           tmp2(1:6,4) = tmp(1:6,4)
           !
           ej2(   1:4) = tmp2(  1,1:4)
           ei2(   1:4) = tmp2(  2,1:4)
           w2(1:4,1:4) = tmp2(3:6,1:4)
           !
           CALL dfpt_tetra_lindhard(ei2,ej2,w2)
           w(1:4,ib,1:4) = w(1:4,ib,1:4) + w2(1:4,1:4) * V
           !       
        END IF
        !          
     ELSE IF((e(2) < 0d0 .AND. 0d0 <= e(3)) .OR. (e(2) <= 0d0 .AND. 0d0 < e(3))) THEN
        !          
        ! C - 1
        !
        V = a(2,4) * a(1,4) * a(3,1)
        !
        IF(V > thr) THEN
           !
           tmp2(1:6,1) = tmp(1:6,1) * a(1,3) + tmp(1:6,3) * a(3,1)
           tmp2(1:6,2) = tmp(1:6,1) * a(1,4) + tmp(1:6,4) * a(4,1)
           tmp2(1:6,3) = tmp(1:6,2) * a(2,4) + tmp(1:6,4) * a(4,2)
           tmp2(1:6,4) = tmp(1:6,4)
           !
           ej2(   1:4) = tmp2(  1,1:4)
           ei2(   1:4) = tmp2(  2,1:4)
           w2(1:4,1:4) = tmp2(3:6,1:4)
           !
           CALL dfpt_tetra_lindhard(ei2,ej2,w2)
           w(1:4,ib,1:4) = w(1:4,ib,1:4) + w2(1:4,1:4) * V
           !      
        END IF
        !
        ! C - 2
        !
        V = a(1,3) * a(2,3)
        !
        IF(V > thr) THEN
           !
           tmp2(1:6,1) = tmp(1:6,1) * a(1,3) + tmp(1:6,3) * a(3,1)
           tmp2(1:6,2) = tmp(1:6,2) * a(2,3) + tmp(1:6,3) * a(3,2)
           tmp2(1:6,3:4) = tmp(1:6,3:4)
           !
           ej2(   1:4) = tmp2(  1,1:4)
           ei2(   1:4) = tmp2(  2,1:4)
           w2(1:4,1:4) = tmp2(3:6,1:4)
           !
           CALL dfpt_tetra_lindhard(ei2,ej2,w2)
           w(1:4,ib,1:4) = w(1:4,ib,1:4) + w2(1:4,1:4) * V
           !
        END IF
        !
        ! C - 3
        ! 
        V = a(1,3) * a(2,4) * a(3,2)
        !
        IF(V > thr) THEN
           !
           tmp2(1:6,1) = tmp(1:6,1) * a(1,3) + tmp(1:6,3) * a(3,1)
           tmp2(1:6,2) = tmp(1:6,2) * a(2,3) + tmp(1:6,3) * a(3,2)
           tmp2(1:6,3) = tmp(1:6,2) * a(2,4) + tmp(1:6,4) * a(4,2)
           tmp2(1:6,4) = tmp(1:6,4)
           !
           ej2(   1:4) = tmp2(  1,1:4)
           ei2(   1:4) = tmp2(  2,1:4)
           w2(1:4,1:4) = tmp2(3:6,1:4)
           !
           CALL dfpt_tetra_lindhard(ei2,ej2,w2)
           w(1:4,ib,1:4) = w(1:4,ib,1:4) + w2(1:4,1:4) * V
           !
        END IF
        !          
     ELSE IF((e(3) < 0d0 .AND. 0d0 <= e(4)) .OR. (e(3) <= 0d0 .AND. 0d0 < e(4))) THEN
        !
        ! D - 1
        !
        V = a(3,4) * a(2,4) * a(1,4) 
        !          
        IF(V > thr) THEN
           !
           tmp2(1:6,1) = tmp(1:6,1) * a(1,4) + tmp(1:6,4) * a(4,1)
           tmp2(1:6,2) = tmp(1:6,2) * a(2,4) + tmp(1:6,4) * a(4,2)
           tmp2(1:6,3) = tmp(1:6,3) * a(3,4) + tmp(1:6,4) * a(4,3)
           tmp2(1:6,4) = tmp(1:6,4)
           !          
           ej2(   1:4) = tmp2(  1,1:4)
           ei2(   1:4) = tmp2(  2,1:4)
           w2(1:4,1:4) = tmp2(3:6,1:4)
           !
           CALL dfpt_tetra_lindhard(ei2,ej2,w2)
           w(1:4,ib,1:4) = w(1:4,ib,1:4) + w2(1:4,1:4) * V
           !        
        END IF
        !
     END IF
     !
  END DO
  !
END SUBROUTINE dfpt_tetra2_lindhard
!
!-----------------------------------------------------------------
SUBROUTINE dfpt_tetra_lindhard(ei,ej,w)
  !---------------------------------------------------------------
  !
  ! This routine compute 1 / (e_{k+q} - e_{k})
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: ei(4), ej(4)
  REAL(dp),INTENT(INOUT) :: w(4,4)
  !
  INTEGER :: ii
  REAL(dp) :: tmp(5,4), w2(4), de(4), lnd(4), thr, thr2
  !
  tmp( 1, 1:4) = ej(1:4) - ei(1:4)
  tmp(2:5,1:4) = w(1:4,1:4)
  CALL dfpt_tetra_sort(5, 4, tmp)
  de(   1:4) = tmp(  1,1:4)
  w(1:4,1:4) = tmp(2:5,1:4)
  !
  thr = MAXVAL(de(1:4)) * 1d-3
  thr2 = 1d-8
  !
  DO ii = 1, 4
     IF(de(ii) < thr2) THEN
        IF(ii == 3) THEN
           CALL errore("dfpt_tetra_lindhard", "Nesting occurs.", 0)
        END IF
        lnd(ii) = 0d0
        de(ii) = 0d0
     ELSE
        lnd(ii) = LOG(de(ii))
     END IF
  END DO
  !
  IF(ABS(de(4) - de(3)) < thr ) THEN
     IF(ABS(de(4) - de(2)) < thr ) THEN
        IF(ABS(de(4) - de(1)) < thr ) THEN
           !
           ! de(4) = de(3) = de(2) = de(1)
           !
           w2(4) = 0.25d0 / de(4)
           w2(3) = w2(4)
           w2(2) = w2(4)
           w2(1) = w2(4)
           !
        ELSE
           !
           ! de(4) = de(3) = de(2)
           !
           w2(4) = dfpt_tetra_lindhard_1211(de(4),de(1),lnd(4),lnd(1))
           w2(3) = w2(4)
           w2(2) = w2(4)
           w2(1) = dfpt_tetra_lindhard_1222(de(1),de(4),lnd(1),lnd(4))
           !
           IF(ANY(w2(1:4) < 0d0)) THEN
              WRITE(*,'(100e15.5)') de(1:4)
              WRITE(*,'(100e15.5)') w2(1:4)
              CALL errore("dfpt_tetra_lindhard", "4=3=2", 0)
           END IF
           !
        END IF
     ELSE IF(ABS(de(2) - de(1)) < thr ) THEN
        !
        ! de(4) = de(3), de(2) = de(1)
        !
        w2(4) = dfpt_tetra_lindhard_1221(de(4),de(2), lnd(4),lnd(2))
        w2(3) = w2(4)
        w2(2) = dfpt_tetra_lindhard_1221(de(2),de(4), lnd(2),lnd(4))
        w2(1) = w2(2)
        !
        IF(ANY(w2(1:4) < 0d0)) THEN
           WRITE(*,'(100e15.5)') de(1:4)
           WRITE(*,'(100e15.5)') w2(1:4)
           CALL errore("dfpt_tetra_lindhard", "4=3 2=1", 0)
        END IF
        !
     ELSE
        !
        ! de(4) = de(3)
        !
        w2(4) = dfpt_tetra_lindhard_1231(de(4),de(1),de(2),lnd(4),lnd(1),lnd(2))
        w2(3) = w2(4)
        w2(2) = dfpt_tetra_lindhard_1233(de(2),de(1),de(4),lnd(2),lnd(1),lnd(4))
        w2(1) = dfpt_tetra_lindhard_1233(de(1),de(2),de(4),lnd(1),lnd(2),lnd(4))
        !
        IF(ANY(w2(1:4) < 0d0)) THEN
           WRITE(*,'(100e15.5)') de(1:4)
           WRITE(*,'(100e15.5)') w2(1:4)
           CALL errore("dfpt_tetra_lindhard", "4=3", 0)
        END IF
        !
     END IF
  ELSE IF(ABS(de(3) - de(2)) < thr) THEN
     IF(ABS(de(3) - de(1)) < thr) THEN
        !
        ! de(3) = de(2) = de(1)
        !
        w2(4) = dfpt_tetra_lindhard_1222(de(4),de(3), lnd(4),lnd(3))
        w2(3) = dfpt_tetra_lindhard_1211(de(3),de(4), lnd(3),lnd(4))
        w2(2) = w2(3)
        w2(1) = w2(3)
        !
        IF(ANY(w2(1:4) < 0d0)) THEN
           WRITE(*,'(100e15.5)') de(1:4)
           WRITE(*,'(100e15.5)') w2(1:4)
           CALL errore("dfpt_tetra_lindhard", "3=2=1", 0)
        END IF
        !
     ELSE
        !
        ! de(3) = de(2)
        !
        w2(4) = dfpt_tetra_lindhard_1233(de(4),de(1),de(3),lnd(4),lnd(1),lnd(3))
        w2(3) = dfpt_tetra_lindhard_1231(de(3),de(1),de(4),lnd(3),lnd(1),lnd(4))
        w2(2) = w2(3)
        w2(1) = dfpt_tetra_lindhard_1233(de(1),de(4),de(3),lnd(1),lnd(4),lnd(3))
        !
        IF(ANY(w2(1:4) < 0d0)) THEN
           WRITE(*,'(100e15.5)') de(1:4)
           WRITE(*,'(100e15.5)') w2(1:4)
           CALL errore("dfpt_tetra_lindhard", "3=2", 0)
        END IF
        !
     END IF
  ELSE IF(ABS(de(2) - de(1)) < thr) THEN
     !
     ! de(2) = de(1)
     !
     w2(4) = dfpt_tetra_lindhard_1233(de(4),de(3),de(2),lnd(4),lnd(3),lnd(2))
     w2(3) = dfpt_tetra_lindhard_1233(de(3),de(4),de(2),lnd(3),lnd(4),lnd(2))
     w2(2) = dfpt_tetra_lindhard_1231(de(2),de(3),de(4),lnd(2),lnd(3),lnd(4))
     w2(1) = w2(2)
     !
     IF(ANY(w2(1:4) < 0d0)) THEN
        WRITE(*,'(100e15.5)') de(1:4)
        WRITE(*,'(100e15.5)') w2(1:4)
        CALL errore("dfpt_tetra_lindhard", "2=1", 0)
     END IF
     !
  ELSE
     !
     ! DIFferent each other.
     !
     w2(4) = dfpt_tetra_lindhard_1234(de(4),de(1),de(2),de(3),lnd(4),lnd(1),lnd(2),lnd(3))
     w2(3) = dfpt_tetra_lindhard_1234(de(3),de(1),de(2),de(4),lnd(3),lnd(1),lnd(2),lnd(4))
     w2(2) = dfpt_tetra_lindhard_1234(de(2),de(1),de(3),de(4),lnd(2),lnd(1),lnd(3),lnd(4))
     w2(1) = dfpt_tetra_lindhard_1234(de(1),de(2),de(3),de(4),lnd(1),lnd(2),lnd(3),lnd(4))
     !      
     IF(ANY(w2(1:4) < 0d0)) THEN
        WRITE(*,'(100e15.5)') de(1:4)
        WRITE(*,'(100e15.5)') w2(1:4)
        CALL errore("dfpt_tetra_lindhard", "Something wrong.", 0)
     END IF
     !
  END IF
  !
  DO ii = 1, 4
     w(1:4,ii) = w2(ii) * w(1:4,ii)
  END DO
  !
END SUBROUTINE dfpt_tetra_lindhard
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1234(g1,g2,g3,g4,lng1,lng2,lng3,lng4) RESULT(w)
  !--------------------------------------------------------------------------
  !
  ! g1, g2, g3, g4 are different each other
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1,g2,g3,g4,lng1,lng2,lng3,lng4
  REAL(dp) :: w
  !
  REAL(dp) :: w2, w3, w4
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3/(g3 - g1)
  w4 = ((lng4 - lng1)/(g4 - g1)*g4 - 1d0)*g4/(g4 - g1)
  w2 = ((w2 - w3)*g2)/(g2 - g3)
  w4 = ((w4 - w3)*g4)/(g4 - g3)
  w = (w4 - w2)/(g4 - g2)
  !
END FUNCTION dfpt_tetra_lindhard_1234
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1231(g1,g2,g3,lng1,lng2,lng3) RESULT(w)
  !----------------------------------------------------------------------------
  !
  ! g4 = g1
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1,g2,g3,lng1,lng2,lng3
  REAL(dp) :: w
  !
  REAL(dp) :: w2, w3
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2**2/(g2 - g1) - g1/( &
  &   2d0)
  w2 = w2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3**2/(g3 - g1) - g1/( &
  &   2d0)
  w3 = w3/(g3 - g1)
  w = (w3 - w2)/(g3 - g2)
  !
END FUNCTION dfpt_tetra_lindhard_1231
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1233(g1,g2,g3,lng1,lng2,lng3) RESULT(w)
  !--------------------------------------------------------------------------
  !
  ! g4 = g3
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1,g2,g3,lng1,lng2,lng3
  REAL(dp) :: w
  !
  REAL(dp) :: w2, w3
  !
  w2 = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w2 = (g2*w2)/(g2 - g1)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = (g3*w3)/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = 1d0 - (2d0*w3*g1)/(g3 - g1)
  w3 = w3/(g3 - g1)
  w = (g3*w3 - g2*w2)/(g3 - g2)
  !
END FUNCTION dfpt_tetra_lindhard_1233
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1221(g1,g2,lng1,lng2) RESULT(w)
  !----------------------------------------------------------------------------
  !
  ! g4 = g1 and g3 = g2
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1, g2, lng1, lng2
  REAL(dp) :: w
  !
  w = 1d0 - (lng2 - lng1)/(g2 - g1)*g1
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(g2 - g1)
  w = w/(2d0*(g2 - g1))
  !
END FUNCTION dfpt_tetra_lindhard_1221
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1222(g1,g2,lng1,lng2) RESULT(w)
  !----------------------------------------------------------------------------
  !
  ! g4 = g3 = g2
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1, g2, lng1, lng2
  REAL(dp) :: w
  !
  w = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w = (2d0*g1*w)/(g2 - g1) - 1d0
  w = (3d0*g1*w)/(g2 - g1) + 1d0
  w = w/(2d0*(g2 - g1))
  !
END FUNCTION dfpt_tetra_lindhard_1222
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1211(g1,g2,lng1,lng2) RESULT(w)
  !----------------------------------------------------------------------------
  !
  ! g4 = g3 = g1
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1,g2,lng1,lng2
  REAL(dp) :: w
  !
  w = -1d0 + (lng2 - lng1)/(g2 - g1)*g2
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(2d0*(g2 - g1))
  w = w/(3d0*(g2 - g1))
  !
END FUNCTION dfpt_tetra_lindhard_1211
!
END MODULE dfpt_tetra_mod
