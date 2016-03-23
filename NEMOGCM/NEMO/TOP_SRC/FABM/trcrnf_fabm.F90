MODULE trcrnf_fabm
   !!======================================================================
   !!                    ***  MODULE  trcrnf_fabm  ***
   !! TOP:  river runoff for FABM
   !!       This is currently hard-coding the run-off variables and should
   !!       be replaced by a more flexible structure asap.
   !!=====================================================================
#if defined key_fabm 

   USE oce
   USE dom_oce
   USE sbc_oce
   USE sbcrnf
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE fldread
   USE in_out_manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_rnf_fabm       ! routine called in trcsms_fabm module
   PUBLIC   trc_rnf_fabm_alloc ! routine call in trcsms_fabm module
   PUBLIC   trc_rnf_fabm_init

   INTEGER, PARAMETER         ::   jprnf_fabm = 2  !: number of runoff variables for FABM states
   LOGICAL                    ::   ln_rnf_no3      !: oxidised nitrogen river runoffs attribute specified in a file
   LOGICAL                    ::   ln_rnf_po4      !: phosphate river runoffs attribute specified in a file
   TYPE(FLD_N)                ::   sn_no3_rnf        !: information about the salinities of runoff file to be read
   TYPE(FLD_N)                ::   sn_po4_rnf        !: information about the temperatures of runoff file to be read
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rnf_fabm_b, rnf_fabm  !: before and now T & S runoff contents   [K.m/s & PSU.m/s]   
   TYPE(FLD),        ALLOCATABLE, DIMENSION(:) ::   sf_no3_rnf     ! structure: river runoff salinity (file information, fields read)  
   TYPE(FLD),        ALLOCATABLE, DIMENSION(:) ::   sf_po4_rnf     ! structure: river runoff temperature (file information, fields read)  
   !! * Substitutions  
#  include "domzgr_substitute.h90"  

   CONTAINS
  
      INTEGER FUNCTION trc_rnf_fabm_alloc()
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE sbc_rnf_alloc  ***
      !!----------------------------------------------------------------------
         ALLOCATE( rnf_fabm_b(jpi,jpj,jpts) , rnf_fabm (jpi,jpj,jpts) , STAT=trc_rnf_fabm_alloc )
         !
         IF( lk_mpp            )   CALL mpp_sum ( trc_rnf_fabm_alloc )
         IF( trc_rnf_fabm_alloc > 0 )   CALL ctl_warn('trc_rnf_fabm_alloc: allocation of arrays failed')
      END FUNCTION trc_rnf_fabm_alloc

      SUBROUTINE trc_rnf_fabm_init

         CONTINUE

      END SUBROUTINE trc_rnf_fabm_init

      SUBROUTINE trc_rnf_fabm( kt )

         INTEGER, INTENT(in) ::   kt          ! ocean time step
         !
         INTEGER  ::   ji, jj    ! dummy loop indices
         INTEGER  ::   z_err = 0 ! dummy integer for error handling
         !!----------------------------------------------------------------------
         IF( kt /= nit000 ) THEN
            rnf_fabm_b(:,:,:) = rnf_fabm(:,:,:)               ! Swap
         ENDIF

         !                                            !   Update runoff   !
         !                                            !-------------------!
         !
         IF(   ln_rnf_no3   )   CALL fld_read ( kt, nn_fsbc, sf_no3_rnf )    ! read runoffs of oxidised nitrogen if required
         IF(   ln_rnf_po4   )   CALL fld_read ( kt, nn_fsbc, sf_po4_rnf )    ! idem runoffs of phosphate if required
 

      END SUBROUTINE trc_rnf_fabm

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rnf_fabm (kt)              ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_rnf_fabm: You should not have seen this print! error?', kt
   END SUBROUTINE trc_rnf_fabm
#endif
   
   !!======================================================================
END MODULE trcrnf_fabm
