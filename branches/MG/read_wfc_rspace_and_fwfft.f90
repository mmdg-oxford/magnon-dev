!
! This routine reads a wavefunction in real space and transform it in Fourier space
!
! Not tested for the non-collinear case.
!
!     Matteo Calandra
!
subroutine read_wfc_rspace_and_fwfft( evc , ik , lrec ,  iunit , n_plane_waves , igmap )
  use kinds,           ONLY : DP
  use wvfct,       ONLY : npwx, nbnd
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE fft_base,            ONLY : dffts, cscatter_smooth
  USE fft_interfaces,      ONLY : fwfft
  USE gvecs,                ONLY : nls
  USE io_global,             ONLY : ionode_id, ionode
  USE mp_global,            ONLY : inter_pool_comm, mpime
  USE mp,                   ONLY : mp_bcast
  
!
  INTEGER, INTENT (IN)      :: ik                    ! k-point to read
  INTEGER, INTENT (IN)      :: lrec                  ! length of the record
  INTEGER, INTENT (IN)      :: n_plane_waves         ! number of plane waves
  INTEGER, INTENT (IN)      :: iunit                 ! input iunit from where to read
  INTEGER, INTENT (IN)      :: igmap(npwx)           ! index for the mapping of the g
  COMPLEX(DP), INTENT (OUT) :: evc(npol*npwx,nbnd)   ! wavefunction in g space
! internal
  INTEGER                   :: ibnd, ig
  COMPLEX(DP), ALLOCATABLE  :: evc_r(:,:), dist_evc_r(:,:)

  allocate( evc_r( dffts%nnr, npol ) )
  allocate( dist_evc_r( dffts%nr1x*dffts%nr2x*dffts%nr3x , nspin_mag) )
  
  !
  ! Fourier transform it in reciprocal space
  !

  do ibnd=1,nbnd
     !
     ! read wfc in real space
     !
#ifdef __MPI
     !
     ! ... First task reads and broadcasts ddrho to all pools
     !

     IF ( ionode ) &
          CALL davcio (dist_evc_r, lrec, iunit, (ik-1)*nbnd+ibnd, - 1)
     
     CALL mp_bcast( dist_evc_r, ionode_id, inter_pool_comm )
          !
     ! ... distributes ddrho between between the tasks of the pool
     !
     DO is = 1, nspin_mag
        !
        CALL cscatter_smooth ( dist_evc_r(:,is), evc_r(:,is) )
        !
     END DO
     !
     !     call mp_bcast( evc_r, ionode_id, inter_pool_comm )
#else    
     CALL davcio (evc_r, lrec, iunit, (ik-1)*nbnd+ibnd, - 1)
#endif
     
     call fwfft('Wave',evc_r(:,1),dffts)
     do ig = 1, n_plane_waves
        evc (ig,ibnd) = evc_r (nls (igmap (ig) ), 1 )
     enddo
     
     IF (noncolin) THEN
        CALL fwfft ('Wave', evc_r(:,2), dffts)
        DO ig = 1, n_plane_waves
           evc (ig+npwx,ibnd) = evc_r (nls(igmap(ig)),2)
        ENDDO
     ENDIF
     
  enddo

  deallocate( evc_r )
  deallocate( dist_evc_r )

end subroutine read_wfc_rspace_and_fwfft
