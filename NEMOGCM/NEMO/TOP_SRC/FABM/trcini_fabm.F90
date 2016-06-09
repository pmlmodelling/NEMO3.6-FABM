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
   USE par_fabm
   USE trcsms_fabm
   USE fabm_config,ONLY: fabm_create_model_from_yaml_file
   USE fabm,ONLY: type_external_variable, fabm_initialize_library
   USE fabm_version,ONLY: fabm_commit_id=>git_commit_id, &
                          fabm_branch_name=>git_branch_name
   USE fabm_types,ONLY: type_version,first_module_version


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
      
         WRITE (xml_unit,1000) ' <field_group id="sf_T" grid_ref="grid_T_2D">'
         DO jn=1,jp_fabm_surface
            CALL write_variable_xml(xml_unit,model%surface_state_variables(jn))
         END DO
         DO jn=1,jp_fabm_bottom
            CALL write_variable_xml(xml_unit,model%bottom_state_variables(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) ' <field_group id="diad_T" grid_ref="grid_T_2D">'
         DO jn=1,size(model%diagnostic_variables)
            CALL write_variable_xml(xml_unit,model%diagnostic_variables(jn),3)
         END DO
         DO jn=1,size(model%horizontal_diagnostic_variables)
            CALL write_variable_xml(xml_unit,model%horizontal_diagnostic_variables(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) ' <field_group id="fabm_scalar" grid_ref="grid_0">'
         DO jn=1,size(model%conserved_quantities)
            CALL write_variable_xml(xml_unit,model%conserved_quantities(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) '</field_definition>'

         CLOSE(xml_unit)
      END IF
      IF( lk_mpp )   CALL mppsync !Ensure field_def_fabm is ready.

1000 FORMAT (A)

   END SUBROUTINE nemo_fabm_init

   SUBROUTINE write_variable_xml(xml_unit,variable,flag_grid_ref)
      INTEGER,INTENT(IN) :: xml_unit
      INTEGER,INTENT(IN),OPTIONAL :: flag_grid_ref
      CLASS (type_external_variable),INTENT(IN) :: variable

      CHARACTER(LEN=20) :: missing_value,string_dimensions
      INTEGER :: number_dimensions

      ! Check variable dimension for grid_ref specificaiton.
      ! Default is to not specify the grid_ref in the field definition.
      IF (present(flag_grid_ref)) THEN
          number_dimensions=flag_grid_ref
      ELSE
          number_dimensions=-1 !default, don't specify grid_ref
      ENDIF

      WRITE (missing_value,'(E9.3)') variable%missing_value
      WRITE (string_dimensions,'(I1)') number_dimensions
      SELECT CASE (number_dimensions)
      CASE (3)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
      CASE (2)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_2D"/>'
      CASE (0)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="1point"/>'
      CASE (-1)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
      CASE default
         IF(lwp) WRITE(numout,*) ' trc_ini_fabm: Failing to initialise output of variable '//TRIM(variable%name)//': Output of '//TRIM(ADJUSTL(string_dimensions))//'-dimensional variables not supported!!!'
      END SELECT
          
   END SUBROUTINE write_variable_xml

   SUBROUTINE trc_ini_fabm
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_fabm  ***  
      !!
      !! ** Purpose :   initialization for FABM model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------

      TYPE (type_version),POINTER :: version
      INTEGER :: jn

      !                       ! Allocate FABM arrays
      IF( trc_sms_fabm_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_fabm: unable to allocate FABM arrays' )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_fabm: initialisation of FABM model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) ' FABM version:   ',fabm_commit_id,' (',fabm_branch_name,' branch)'

      call fabm_initialize_library()
      version => first_module_version
      
      do while (associated(version))
         IF(lwp) WRITE(numout,*)  ' '//trim(version%module_name)//' version:   ',trim(version%version_string)
         version => version%next
      end do
      
      ! Log mapping of FABM states:
      IF (lwp) THEN
         IF (jp_fabm.gt.0) WRITE(numout,*) " FABM tracers:"
         DO jn=1,jp_fabm
            WRITE(numout,*) "   State",jn,":",trim(model%state_variables(jn)%name), &
               " (",trim(model%state_variables(jn)%long_name), &
               ") [",trim(model%state_variables(jn)%units),"]"
         ENDDO
         IF (jp_fabm_surface.gt.0) WRITE(numout,*) "FABM seasurface states:"
         DO jn=1,jp_fabm_surface
            WRITE(numout,*) "   State",jn,":",trim(model%surface_state_variables(jn)%name), &
               " (",trim(model%surface_state_variables(jn)%long_name), &
               ") [",trim(model%surface_state_variables(jn)%units),"]"
         ENDDO
         IF (jp_fabm_bottom.gt.0) WRITE(numout,*) "FABM seafloor states:"
         DO jn=1,jp_fabm_bottom
            WRITE(numout,*) "   State",jn,":",trim(model%bottom_state_variables(jn)%name), &
               " (",trim(model%bottom_state_variables(jn)%long_name), &
               ") [",trim(model%bottom_state_variables(jn)%units),"]"
         ENDDO
      ENDIF
      
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
