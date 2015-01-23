!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM magnon
  !-----------------------------------------------------------------------
  ! ... This is the main driver of the phonon code.
  ! ... It reads all the quantities calculated by pwscf, it
  ! ... checks if some recover file is present and determines
  ! ... which calculation needs to be done. Finally, it makes
  ! ... a loop over the q points. At a generic q, if necessary it
  ! ... recalculates the band structure calling pwscf again.
  !     At q=0 it can calculate the linear response
  ! ... to an electric field perturbation and hence the dielectric
  ! ... constant, the Born effective charges and the polarizability
  ! ... at imaginary frequencies.
  ! ... At q=0, from the second order response to an electric field,
  ! ... it can calculate also the electro-optic and the raman tensors.
  ! ... Presently implemented:
  ! ... dynamical matrix (q/=0)   NC [4], US [4], PAW [4]
  ! ... dynamical matrix (q=0)    NC [5], US [5], PAW [4]
  ! ... dielectric constant       NC [5], US [5], PAW [3]
  ! ... born effective charges    NC [5], US [5], PAW [3]
  ! ... polarizability (iu)       NC [2], US [2]
  ! ... electron-phonon           NC [3], US [3]
  ! ... electro-optic             NC [1]
  ! ... raman tensor              NC [1]
  !
  ! NC = norm conserving pseudopotentials
  ! US = ultrasoft pseudopotentials
  ! PAW = projector augmented-wave
  ! [1] LDA, 
  ! [2] [1] + GGA, 
  ! [3] [2] + LSDA/sGGA, 
  ! [4] [3] + Spin-orbit/nonmagnetic,
  ! [5] [4] + Spin-orbit/magnetic (experimental when available)
  !
  ! Not implemented in ph.x:
  ! [6] [5] + constraints on the magnetization
  ! [7] [6] + Hubbard U
  ! [8] [7] + Hybrid functionals
  ! [9]  ?  + External Electric field
  ! [10] ?  + nonperiodic boundary conditions.

  USE io_global,       ONLY : stdout
  USE disp,            ONLY : nqs, num_k_pts, xk_kpoints, comp_iq
  USE control_ph,      ONLY : epsil, trans, bands_computed,lgamma
  USE output,          ONLY : fildrho
  USE check_stop,      ONLY : check_stop_init
  USE ph_restart,      ONLY : ph_writefile, destroy_status_run
  USE save_ph,         ONLY : clean_input_variables
  USE mp_global,       ONLY: mp_startup, nimage
  USE image_io_routines, ONLY : io_image_start
  USE environment,     ONLY: environment_start
  USE qpoint, ONLY:xq
  USE freq_ph,       ONLY : fpol, fiu, nfs, nfsmax

  !
  IMPLICIT NONE
  !
  INTEGER :: iq,iq1,i
  LOGICAL :: do_band, do_iq, setup_pw
  CHARACTER (LEN=9)   :: code = 'MAGNON'
  CHARACTER (LEN=256) :: auxdyn
  !
  ! Initialize MPI, clocks, print initial messages
  !
#ifdef __MPI
  CALL mp_startup ( start_images=.true. )
  IF (nimage>1) CALL io_image_start()
#endif
  CALL environment_start ( code )
  !
  WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
  !
  ! ... and begin with the initialization part
  !
  WRITE(stdout, '(/5x, "Reading variables")') 
  CALL phq_readin()
  
  !do iq=1, num_k_pts
  !WRITE(stdout,*)xk_kpoints(1,iq), xk_kpoints(2,iq), xk_kpoints(3,iq)
  !end do 
  WRITE(stdout, '(/5x, "Finished reading variables")') 
  WRITE(stdout, '(7x, "Imag. Frequencies: ")')
  DO i = 1, nfs
       WRITE(stdout,'(8x, i4, 4x, 2f9.4)')i, fiu(i)*13.605
  ENDDO

  !
  CALL check_stop_init()
  !
  ! ... Checking the status of the calculation and if necessary initialize
  ! ... the q mesh
  !
  CALL check_initial_status(auxdyn)
! WRITE(stdout,*)'xq'
! WRITE(stdout,*)xq
! WRITE(stdout,*)'lgamma',lgamma

  WRITE(stdout,*)'comp_iq'
  do iq=1, num_k_pts
  WRITE(stdout,*)comp_iq(iq)
  end do
  !
  !DO iq = 1, nqs
  
  DO iq = 1, num_k_pts
     !
     print*, "Number of qpoints to calc.", num_k_pts

     CALL prepare_q(do_band, do_iq, setup_pw, iq)
    ! WRITE(stdout,)
     !
     !  If necessary the bands are recalculated

     !
     IF (setup_pw) CALL run_pwscf(do_band)
     !
     !  Initialize the quantities which do not depend on
     !  the linear response of the system
     !
     CALL initialize_ph()

     !
     !  magnon perturbation
     !
     CALL magscf()
    !calculates dielectric tensor:
    !CALL phescf()
     CALL clean_pw_ph(iq)
  END DO

  CALL ph_writefile('init',0)
  CALL clean_input_variables()
  CALL collect_grid_files()
  CALL destroy_status_run()
  !
  IF (bands_computed) CALL print_clock_pw()
  !
  CALL stop_ph( .TRUE. )
  !
  STOP
  !
END PROGRAM magnon
