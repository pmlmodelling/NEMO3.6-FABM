MODULE trcwri_fabm
   !!======================================================================
   !!                       *** MODULE trcwri_fabm ***
   !!    fabm :   Output of FABM tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && key_fabm && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                           FABM model
   !!----------------------------------------------------------------------
   !! trc_wri_fabm   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE trcsms_fabm, only: trc_sms_fabm_check_mass
   USE par_fabm
   USE st2d_fabm
   USE fabm, only: fabm_get_bulk_diagnostic_data, fabm_get_horizontal_diagnostic_data

   IMPLICIT NONE
   PRIVATE

#if defined key_tracer_budget
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: trb_temp ! slwa
#endif

   INTERFACE trc_wri_fabm
       MODULE PROCEDURE wri_fabm,wri_fabm_fl
   END INTERFACE trc_wri_fabm


   PUBLIC trc_wri_fabm 

#  include "top_substitute.h90"
CONTAINS

   SUBROUTINE wri_fabm_fl(kt,fl)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )               :: fl
      INTEGER, INTENT( in )               :: kt

#if defined key_tracer_budget
      INTEGER              :: jn
      CHARACTER (len=20)   :: cltra
      REAL(wp), DIMENSION(jpi,jpj,jpk)    :: trpool,st2dpool !temporary storage pools
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
! depth integrated
! for strict budgetting write this out at end of timestep as an average between 'now' and 'after' at kt
      DO jn = 1, jp_fabm1
         trpool(:,:,:) = 0.5 * (trn(:,:,:,jp_fabm0+jn-1)*fse3t_a(:,:,:) +  &
                             tr_temp(:,:,:jp_fabm0+jn-1jn)*fse3t(:,:,:) )
         cltra = TRIM( model%state_variables(jn)%name )//"e3t"     ! depth integrated output
         IF( kt == nittrc000 ) write(6,*)'output pool ',cltra
         CALL iom_put( cltra, trpool)
      END DO
      DO jn = 1, jp_fabm1_surface
         st2dpool(:,:) = 0.5 * (fabm_st2dn(:,:,jn) +  &
                             fabm_st2d_temp(:,:,jn))
         cltra = TRIM( model%surface_state_variables(jn)%name )//"e3t"     ! depth integrated output
         IF( kt == nittrc000 ) write(6,*)'output pool ',cltra
         CALL iom_put( cltra, st2dpool)
      END DO
      DO jn = 1, jp_fabm1_bottom
         st2dpool(:,:) = 0.5 * (fabm_st2dn(:,:,jp_fabm_surface+jn) +  &
                             fabm_st2d_temp(:,:,jp_fabm_surface+jn))
         cltra = TRIM( model%surface_state_variables(jn)%name )//"e3t"     ! depth integrated output
         IF( kt == nittrc000 ) write(6,*)'output pool ',cltra
         CALL iom_put( cltra, st2dpool)
      END DO
#else
      CONTINUE
#endif

   END SUBROUTINE wri_fabm_fl

   SUBROUTINE wri_fabm(kt)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )               :: kt
      INTEGER              :: jn

#if defined key_tracer_budget
      IF( kt == nittrc000 ) THEN
         ALLOCATE(tr_temp(jpi,jpj,jpk,jp_my_trc),fabm_st2d_temp(jpi,jpj,jp_fabm1_surface+jp_fabm_bottom)  ! slwa
      ENDIF
      tr_temp(:,:,:,:)=trn(:,:,:,:) ! slwa save for tracer budget (unfiltered trn)
      fabm_st2d_temp(:,:,:,:)=fabm_st2dn(:,:,:,:)
#else
      DO jn = 1, jp_fabm
         CALL iom_put( model%state_variables(jn)%name, trn(:,:,:,jp_fabm0+jn-1) )
      END DO
      DO jn = 1, jp_fabm_surface
         CALL iom_put( model%surface_state_variables(jn)%name, fabm_st2dn(:,:,jn) )
      END DO
      DO jn = 1, jp_fabm_bottom
         CALL iom_put( model%bottom_state_variables(jn)%name, fabm_st2dn(:,:,jp_fabm_surface+jn) )
      END DO
#endif
      ! write 3D diagnostics in the file
      ! ---------------------------------------
      DO jn = 1, size(model%diagnostic_variables)
         IF (model%diagnostic_variables(jn)%save) &
             CALL iom_put( model%diagnostic_variables(jn)%name, fabm_get_bulk_diagnostic_data(model,jn))
      END DO

      ! write 2D diagnostics in the file
      ! ---------------------------------------
      DO jn = 1, size(model%horizontal_diagnostic_variables)
         IF (model%horizontal_diagnostic_variables(jn)%save) &
             CALL iom_put( model%horizontal_diagnostic_variables(jn)%name, fabm_get_horizontal_diagnostic_data(model,jn))
      END DO
      !
      CALL trc_sms_fabm_check_mass

   END SUBROUTINE wri_fabm

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_fabm
CONTAINS
   SUBROUTINE trc_wri_fabm                     ! Empty routine  
   END SUBROUTINE trc_wri_fabm
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri_fabm.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_fabm
