MODULE trcini_fabm
   !!======================================================================
   !!                         ***  MODULE trcini_fabm  ***
   !! TOP :   initialisation of the FABM tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                               FABM tracers
   !!----------------------------------------------------------------------
   !! trc_ini_fabm   : FABM model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE trcsms_fabm
   USE fabm_config,ONLY: fabm_create_model_from_yaml_file
   USE fabm,ONLY: type_external_variable

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_fabm   ! called by trcini.F90 module
   PUBLIC   nemo_fabm_init

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE nemo_fabm_init()
      INTEGER :: jn
      INTEGER, PARAMETER :: xml_unit = 1979

      ! Allow FABM to parse fabm.yaml. This ensures numbers of variables are known.
      call fabm_create_model_from_yaml_file(model)

      jp_fabm = size(model%state_variables)
      jp_fabm_bottom = size(model%bottom_state_variables)
      jp_fabm_surface = size(model%surface_state_variables)
      jp_fabm0 = jptra + 1
      jp_fabm1 = jptra + jp_fabm
      jp_fabm_m1=jptra
      jptra = jptra + jp_fabm
      jpdia2d = jpdia2d + size(model%horizontal_diagnostic_variables)
      jpdia3d = jpdia3d + size(model%diagnostic_variables)
      jpdiabio = jpdiabio + jp_fabm
      
      IF (lwp) THEN
         ! write field_def_fabm.xml on lead process
         OPEN(UNIT=xml_unit,FILE='field_def_fabm.xml',ACTION='WRITE',STATUS='REPLACE')

         WRITE (xml_unit,1000) '<field_definition level="1" prec="4" operation="average" enabled=".TRUE." default_value="1.e20" >'

         WRITE (xml_unit,1000) ' <field_group id="ptrc_T" grid_ref="grid_T_3D">'
         DO jn=1,jp_fabm
            CALL write_variable_xml(xml_unit,model%state_variables(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'
      
         WRITE (xml_unit,1000) ' <field_group id="sf1D_T" grid_ref="grid_T_2D">'
         DO jn=1,jp_fabm_surface
            CALL write_variable_xml(xml_unit,model%surface_state_variables(jn))
         END DO
         DO jn=1,jp_fabm_bottom
            CALL write_variable_xml(xml_unit,model%bottom_state_variables(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) ' <field_group id="diad_T" grid_ref="grid_T_2D">'
         DO jn=1,size(model%diagnostic_variables)
            CALL write_variable_xml(xml_unit,model%diagnostic_variables(jn),.TRUE.)
         END DO
         DO jn=1,size(model%horizontal_diagnostic_variables)
            CALL write_variable_xml(xml_unit,model%horizontal_diagnostic_variables(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) '</field_definition>'

         CLOSE(xml_unit)
      END IF
      IF( lk_mpp )   CALL mppsync !Ensure field_def_fabm is ready.

1000 FORMAT (A)

   END SUBROUTINE nemo_fabm_init

   SUBROUTINE write_variable_xml(xml_unit,variable,flag3D)
      INTEGER,INTENT(IN) :: xml_unit
      LOGICAL,INTENT(IN),OPTIONAL :: flag3D
      CLASS (type_external_variable),INTENT(IN) :: variable

      CHARACTER(LEN=20) :: missing_value
      LOGICAL :: is3D

      ! Check variable dimension, default is 2D:
      IF (present(flag3D)) THEN
          is3D=flag3D
      ELSE
          is3D=.FALSE.
      ENDIF

      WRITE (missing_value,'(E9.3)') variable%missing_value
      IF (is3D) THEN
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
      ELSE
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
      ENDIF
          
   END SUBROUTINE write_variable_xml

   SUBROUTINE trc_ini_fabm
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_fabm  ***  
      !!
      !! ** Purpose :   initialization for FABM model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------

      !                       ! Allocate FABM arrays
      IF( trc_sms_fabm_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_fabm: unable to allocate FABM arrays' )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_fabm: initialisation of FABM model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      
   END SUBROUTINE trc_ini_fabm

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No FABM model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE nemo_fabm_init
   END SUBROUTINE nemo_fabm_init

   SUBROUTINE trc_ini_fabm            ! Empty routine
   END SUBROUTINE trc_ini_fabm
#endif

   !!======================================================================
END MODULE trcini_fabm
