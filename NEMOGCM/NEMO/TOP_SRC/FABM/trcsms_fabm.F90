MODULE trcsms_fabm
   !!======================================================================
   !!                         ***  MODULE trcsms_fabm  ***
   !! TOP :   Main module of the FABM tracers
   !!======================================================================
   !! History :   1.0  !  2015-04  (PML) Original code
   !!----------------------------------------------------------------------
#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                               FABM tracers
   !!----------------------------------------------------------------------
   !! trc_sms_fabm       : FABM model main routine
   !! trc_sms_fabm_alloc : allocate arrays specific to FABM sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trcbc
   USE trd_oce
   USE trdtrc

   USE oce, only: tsn  ! Needed?
   USE dom_oce
   USE iom

   USE st2D_fabm

   USE fldread         !  time interpolation

   USE fabm

   IMPLICIT NONE

#  include "domzgr_substitute.h90"

   PRIVATE

   PUBLIC   trc_sms_fabm       ! called by trcsms.F90 module
   PUBLIC   trc_sms_fabm_alloc ! called by trcini_fabm.F90 module
   PUBLIC   trc_sms_fabm_check_mass
   PUBLIC   model
   PUBLIC   st2d_fabm_nxt ! 2D state intergration

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: flux    ! Cross-interface flux of pelagic variables (# m-2 s-1)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   :: ext     ! Light extinction coefficient (m-1)

   ! Work arrays for vertical advection (residual movement/sinking/floating)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:) :: w_ct
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   :: w_if
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   :: zwgt_if
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   :: flux_if
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   :: flux_ct
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   :: current_total

   ! Arrays for environmental variables
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:) :: prn,rho

   REAL(wp), PUBLIC :: daynumber_in_year

   TYPE type_input_variable
      TYPE (type_horizontal_variable_id)   :: horizontal_id
      TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf
      INTEGER                              :: ntimes
      TYPE(type_input_variable), POINTER   :: next => null()
   END TYPE
   TYPE (type_input_variable), POINTER, SAVE :: first_input_variable => NULL()

   TYPE type_river_data
      INTEGER   :: jp_pos=0 ! position of linked state variable in trc fields
      TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf
      INTEGER                              :: ntimes
      TYPE(type_river_data), POINTER   :: next => null()
      REAL(wp) :: rn_trrnfac=1._wp ! unit conversion factor
   END TYPE

   TYPE (type_river_data), POINTER, SAVE :: first_river_data => NULL()

   TYPE (type_bulk_variable_id),SAVE :: swr_id

   TYPE (type_model) :: model

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_fabm( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_fabm  ***
      !!
      !! ** Purpose :   main routine of FABM model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER :: jn
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztrfabm

!!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_fabm')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_fabm:  FABM model, iteration',kt
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      CALL update_inputs( kt )

      CALL compute_rates

      CALL compute_vertical_movement

      CALL st2d_fabm_nxt( kt )

      IF( l_trdtrc )  CALL wrk_alloc( jpi, jpj, jpk, ztrfabm )

      CALL trc_bc_read  ( kt )       ! tracers: surface and lateral Boundary Conditions
      CALL trc_rnf_fabm ( kt ) ! River forcings

      IF( l_trdtrc ) THEN      ! Save the trends in the mixed layer
          DO jn = jp_fabm0, jp_fabm1
            ztrfabm(:,:,:) = tra(:,:,:,jn)
            CALL trd_trc( ztrfabm, jn, jptra_sms, kt )   ! save trends
          END DO
          CALL wrk_dealloc( jpi, jpj, jpk, ztrfabm )
      END IF
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_fabm')

   END SUBROUTINE trc_sms_fabm

   SUBROUTINE compute_rates()
      INTEGER :: ji,jj,jk,jn
      LOGICAL :: valid_state
      REAL(wp) :: zalfg

      ! Validate current model state (setting argument to .TRUE. enables repair=clipping)
      valid_state = check_state(.FALSE.)

      ! Compute the now hydrostatic pressure
      ! copied from istate.F90
      ! ------------------------------------

      zalfg = 0.5e-4 * grav ! FABM wants dbar, convert from Pa
      
      rho = rau0 * ( 1. + rhd )

      prn(:,:,1) = 10.1325 + zalfg * fse3t(:,:,1) * rho(:,:,1)

      daynumber_in_year=(fjulday-fjulstartyear+1)*1._wp

      DO jk = 2, jpk                                              ! Vertical integration from the surface
         prn(:,:,jk) = prn(:,:,jk-1) + zalfg * ( &
                     fse3t(:,:,jk-1) * rho(:,:,jk-1)  &
                     + fse3t(:,:,jk) * rho(:,:,jk) )
      END DO  

      ! Compute light extinction
      DO jk=1,jpk
          DO jj=1,jpj
            call fabm_get_light_extinction(model,1,jpi,jj,jk,ext)
         END DO
      END DO

      ! Compute light field (stored among FABM's internal diagnostics)
      DO jj=1,jpj
          DO ji=1,jpi
            call fabm_get_light(model,1,jpk,ji,jj)
         END DO
      END DO

      ! TODO: retrieve 3D shortwave and store in etot3

      ! Zero rate array of interface-attached state variables
      fabm_st2Da = 0._wp

      ! Compute interfacial source terms and fluxes
      DO jj=1,jpj
         ! Process bottom (fabm_do_bottom increments rather than sets, so zero flux array first)
         flux = 0._wp
         CALL fabm_do_bottom(model,1,jpi,jj,flux,fabm_st2Da(:,jj,jp_fabm_surface+1:))
         DO jn=1,jp_fabm
             DO ji=1,jpi
                 ! Divide bottom fluxes by height of bottom layer and add to source terms.
                 ! TODO: is there perhaps an existing variable for fse3t(ji,jj,mbkt(ji,jj))??
                 tra(ji,jj,mbkt(ji,jj),jp_fabm_m1+jn) = tra(ji,jj,mbkt(ji,jj),jp_fabm_m1+jn) + flux(ji,jn)/fse3t(ji,jj,mbkt(ji,jj))
             END DO
         END DO

         ! Process surface (fabm_do_surface increments rather than sets, so zero flux array first)
         flux = 0._wp
         CALL fabm_do_surface(model,1,jpi,jj,flux,fabm_st2Da(:,jj,1:jp_fabm_surface))
         DO jn=1,jp_fabm
             ! Divide surface fluxes by height of surface layer and add to source terms.
             tra(:,jj,1,jp_fabm_m1+jn) = tra(:,jj,1,jp_fabm_m1+jn) + flux(:,jn)/fse3t(:,jj,1)
         END DO
      END DO

      ! Compute interior source terms (NB fabm_do increments rather than sets)
      DO jk=1,jpk
          DO jj=1,jpj
              CALL fabm_do(model,1,jpi,jj,jk,tra(:,jj,jk,jp_fabm0:jp_fabm1))
          END DO
      END DO
   END SUBROUTINE compute_rates

   SUBROUTINE compute_vertical_movement()
      INTEGER :: ji,jj,jk,jn,k_floor

      ! Compute interior vertical velocities and include them in source array.
      DO jj=1,jpj
         ! Get vertical velocities at layer centres (entire 1:jpi,1:jpk slice).
         DO jk=1,jpk
            CALL fabm_get_vertical_movement(model,1,jpi,jj,jk,w_ct(:,jk,:))
         END DO

         DO ji=1,jpi
            ! Only process this horizontal point (ji,jj) if number of layers exceeds 1
            if (mbkt(ji,jj)>1) THEN
               k_floor=mbkt(ji,jj)
               ! Linearly interpolate to velocities at the interfaces between layers
               zwgt_if(1:k_floor-1,:)=spread(&
                   fse3t(ji,jj,2:k_floor)/ (fse3t(ji,jj,1:k_floor-1)+fse3t(ji,jj,2:k_floor))&
                   ,2,jp_fabm)
               w_if(1:k_floor-1,:) = zwgt_if(1:k_floor-1,:)*w_ct(ji,1:k_floor-1,:)&
                  +(1._wp-zwgt_if(1:k_floor-1,:))*w_ct(ji,2:k_floor,:)

               ! Convert velocities in m s-1 to mass fluxes (e.g., mol m-2 s-1) - both at interfaces between layers.
               DO jn=1,jp_fabm
                  DO jk=1,k_floor-1
                     IF (w_if(jk,jn)<0) THEN
                        ! Downward - use concentration from same level (velocity is defined on bottom interface)
                        flux_if(jk,jn) = w_if(jk,jn)*trn(ji,jj,jk,jp_fabm_m1+jn)
                     ELSE
                        ! Upward - use concentration from level below (velocity is defined on bottom interface)
                        flux_if(jk,jn) = w_if(jk,jn)*trn(ji,jj,jk+1,jp_fabm_m1+jn)
                     END IF
                  END DO
               END DO

               ! Combine interfacial mass fluxes (top + bottom) into total fluxes per layer
               flux_ct(1,              :) = flux_if(1,              :)
               flux_ct(2:k_floor-1,:) = flux_if(2:k_floor-1,:) - flux_if(1:k_floor-2,:)
               flux_ct(k_floor,    :) =                            - flux_if(k_floor-1,  :)                              

               ! Convert mass fluxes (m-2) into source terms (m-3) by dividing by layer height
               DO jn=1,jp_fabm
                  tra(ji,jj,1:k_floor,jp_fabm_m1+jn) = tra(ji,jj,1:k_floor,jp_fabm_m1+jn) + flux_ct(1:k_floor,jn)/fse3t(ji,jj,1:k_floor)
               END DO
            END IF
         END DO
      END DO

   END SUBROUTINE compute_vertical_movement

   FUNCTION check_state(repair) RESULT(valid)
      LOGICAL, INTENT(IN) :: repair
      LOGICAL             :: valid

      INTEGER             :: jj,jk
      LOGICAL             :: valid_int,valid_sf,valid_bt

      valid = .TRUE.
      DO jk=1,jpk
         DO jj=1,jpj
            CALL fabm_check_state(model,1,jpi,jj,jk,repair,valid_int)
            IF (.NOT.(valid_int.OR.repair)) valid = .FALSE.
         END DO
      END DO
      DO jj=1,jpj
         CALL fabm_check_surface_state(model,1,jpi,jj,repair,valid_sf)
         CALL fabm_check_bottom_state(model,1,jpi,jj,repair,valid_bt)
         IF (.NOT.(valid_sf.AND.valid_bt).AND..NOT.repair) valid = .FALSE.
      END DO
   END FUNCTION

   SUBROUTINE trc_sms_fabm_check_mass()
      REAL(wp) :: total(SIZE(model%conserved_quantities))
      INTEGER :: jk,jj,jn

      total = 0._wp

      DO jk=1,jpk
         DO jj=1,jpj
            CALL fabm_get_conserved_quantities(model,1,jpi,jj,jk,current_total)
            DO jn=1,SIZE(model%conserved_quantities)
               total(jn) = total(jn) + SUM(cvol(:,jj,jk)*current_total(:,jn)*tmask_i(:,jj))
            END DO
         END DO
      END DO

      DO jj=1,jpj
         CALL fabm_get_horizontal_conserved_quantities(model,1,jpi,jj,current_total)
         DO jn=1,SIZE(model%conserved_quantities)
            total(jn) = total(jn) + SUM(e1e2t(:,jj)*current_total(:,jn)*tmask_i(:,jj))
         END DO
      END DO

      IF( lk_mpp ) CALL mpp_sum(total,SIZE(model%conserved_quantities))

      DO jn=1,SIZE(model%conserved_quantities)
         IF(lwp) WRITE(numout,*) 'FABM '//TRIM(model%conserved_quantities(jn)%name),total(jn),TRIM(model%conserved_quantities(jn)%units)//'*m3'
      END DO

   END SUBROUTINE trc_sms_fabm_check_mass

   SUBROUTINE st2d_fabm_nxt( kt )
      !!----------------------------------------------------------------------
      !!                     ***  st2d_fabm_nxt  ***
      !!
      !! ** Purpose :   routine to integrate 2d states in time
      !!
      !! ** Method  :   based on integration of 3D passive tracer fields
      !!                implemented in TOP_SRC/TRP/trcnxt.F90, plus
      !!                tra_nxt_fix in OPA_SRC/TRA/tranxt.F90. Similar to
      !!                time integration of sea surface height in
      !!                OPA_SRC/DYN/sshwzv.F90.
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      REAL(wp) :: z2dt
!!----------------------------------------------------------------------
      !
      IF( neuler == 0 .AND. kt == nittrc000 ) THEN
          z2dt = rdt                  ! set time step size (Euler)
      ELSE
          z2dt = 2._wp * rdt          ! set time step size (Leapfrog)
      ENDIF

      ! Forward Euler time step to compute "now"
      fabm_st2Da(:,:,:) = fabm_st2db(:,:,:) + z2dt * fabm_st2da(:,:,:)

      IF( neuler == 0 .AND. kt == nittrc000 ) THEN        ! Euler time-stepping at first time-step
         !                                                ! (only swap)
         fabm_st2Dn(:,:,:) = fabm_st2Da(:,:,:)
         !                                              
      ELSE
         ! Update now state + Asselin filter time stepping
         fabm_st2Db(:,:,:) = (1._wp - 2._wp*atfp) * fabm_st2Dn(:,:,:) + &
             atfp * ( fabm_st2Db(:,:,:) + fabm_st2Da(:,:,:) )
         fabm_st2Dn(:,:,:) = fabm_st2Da(:,:,:)
      ENDIF

   END SUBROUTINE st2d_fabm_nxt

   INTEGER FUNCTION trc_sms_fabm_alloc()
      INTEGER :: jj,jk,jn
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_fabm_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to FABM
      ALLOCATE( lk_rad_fabm(jp_fabm))
      ALLOCATE( prn(jpi, jpj, jpk))
      ALLOCATE( rho(jpi, jpj, jpk))
      ! ALLOCATE( tab(...) , STAT=trc_sms_fabm_alloc )

      ! Allocate arrays to hold state for surface-attached and bottom-attached state variables
      ALLOCATE(fabm_st2Dn(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))
      ALLOCATE(fabm_st2Da(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))
      ALLOCATE(fabm_st2Db(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))

      ! Work array to hold surface and bottom fluxes
      ALLOCATE(flux(jpi,jp_fabm))

      ! Work array to hold extinction coefficients
      ALLOCATE(ext(jpi))
      ext=0._wp

      ! Allocate work arrays for vertical movement
      ALLOCATE(w_ct(jpi,jpk,jp_fabm))
      ALLOCATE(w_if(jpk-1,jp_fabm))
      ALLOCATE(zwgt_if(jpk-1,jp_fabm))
      ALLOCATE(flux_if(jpk-1,jp_fabm))
      ALLOCATE(flux_ct(jpk,jp_fabm))
      ALLOCATE(current_total(jpi,SIZE(model%conserved_quantities)))

      trc_sms_fabm_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_fabm_alloc /= 0 ) CALL ctl_warn('trc_sms_fabm_alloc : failed to allocate arrays')
      !

      ! Make FABM aware of diagnostics that are not needed [not included in output]
      DO jn=1,size(model%diagnostic_variables)
          !model%diagnostic_variables(jn)%save = iom_use(model%diagnostic_variables(jn)%name)
      END DO
      DO jn=1,size(model%horizontal_diagnostic_variables)
          !model%horizontal_diagnostic_variables(jn)%save = iom_use(model%horizontal_diagnostic_variables(jn)%name)
      END DO

      ! Provide FABM with domain extents [after this, the save attribute of diagnostic variables can no longe change!]
      call fabm_set_domain(model,jpi, jpj, jpk)

      ! Provide FABM with the vertical indices of the surface and bottom, and the land-sea mask.
      call model%set_bottom_index(mbkt)  ! NB mbkt extents should match dimension lengths provided to fabm_set_domain
      call model%set_surface_index(1)
      call fabm_set_mask(model,tmask,tmask(:,:,1)) ! NB tmask extents should match dimension lengths provided to fabm_set_domain

      ! Send pointers to state data to FABM
      do jn=1,jp_fabm
         call fabm_link_bulk_state_data(model,jn,trn(:,:,:,jp_fabm_m1+jn))
      end do
      DO jn=1,jp_fabm_surface
         CALL fabm_link_surface_state_data(model,jn,fabm_st2Dn(:,:,jn))
      END DO
      DO jn=1,jp_fabm_bottom
         CALL fabm_link_bottom_state_data(model,jn,fabm_st2Dn(:,:,jp_fabm_surface+jn))
      END DO

      ! Send pointers to environmental data to FABM
      call fabm_link_bulk_data(model,standard_variables%temperature,tsn(:,:,:,jp_tem))
      call fabm_link_bulk_data(model,standard_variables%practical_salinity,tsn(:,:,:,jp_sal))
      call fabm_link_bulk_data(model,standard_variables%density,rho(:,:,:))
      call fabm_link_bulk_data(model,standard_variables%pressure,prn)
      ! correct target for cell thickness depends on NEMO configuration:
#ifdef key_vvl
      call fabm_link_bulk_data(model,standard_variables%cell_thickness,e3t_n)
#else
      call fabm_link_bulk_data(model,standard_variables%cell_thickness,e3t_0)
#endif
      call fabm_link_horizontal_data(model,standard_variables%latitude,gphit)
      call fabm_link_horizontal_data(model,standard_variables%longitude,glamt)
      !call fabm_link_horizontal_data(model,standard_variables%bottom_stress,?) for now ignored.
      call fabm_link_scalar_data(model,standard_variables%number_of_days_since_start_of_the_year,daynumber_in_year)
      call fabm_link_horizontal_data(model,standard_variables%wind_speed,wndm(:,:))
      call fabm_link_horizontal_data(model,standard_variables%surface_downwelling_shortwave_flux,qsr(:,:))
      call fabm_link_horizontal_data(model,standard_variables%bottom_depth_below_geoid,bathy(:,:))

      swr_id = model%get_bulk_variable_id(standard_variables%downwelling_shortwave_flux)

      ! Obtain user-specified input variables (read from NetCDF file)
      call initialize_inputs
      call update_inputs(nit000)

      ! Check whether FABM has all required data
      call fabm_check_ready(model)
      
      ! Initialize state
      DO jj=1,jpj
         CALL fabm_initialize_surface_state(model,1,jpi,jj)
         CALL fabm_initialize_bottom_state(model,1,jpi,jj)
      END DO
      DO jk=1,jpk
         DO jj=1,jpj
            CALL fabm_initialize_state(model,1,jpi,jj,jk)
         END DO
      END DO

      ! Set mask for negativity corrections to the relevant states
      DO jn=1,jp_fabm
        IF (model%state_variables(jn)%minimum.ge.0) THEN
          lk_rad_fabm(jn)=.TRUE.
          IF(lwp) WRITE(numout,*) 'FABM clipping for '//TRIM(model%state_variables(jn)%name)//' activated.'   
        END IF
      END DO

      ! Mask land points in states with zeros, not nice, but coherent 
      ! with NEMO "convention":
      DO jn=jp_fabm0,jp_fabm1
        WHERE (tmask==0._wp)
          trn(:,:,:,jn)=0._wp
        END WHERE
      END DO
      DO jn=1,jp_fabm_surface+jp_fabm_bottom
        WHERE (tmask(:,:,1)==0._wp)
          fabm_st2Dn(:,:,jn)=0._wp
        END WHERE
      END DO

      ! Copy initial condition for interface-attached state variables to "previous" state field
      ! NB NEMO does this itself for pelagic state variables (trb) in TOP_SRC/trcini.F90.
      fabm_st2Db = fabm_st2Dn

      ! Ensure that all FABM diagnostics have a valid value.
      wndm=0._wp !uninitiased field at this point
      qsr=0._wp !uninitiased field at this point
      CALL compute_rates

   END FUNCTION trc_sms_fabm_alloc

   SUBROUTINE initialize_inputs
      TYPE(FLD_N)        :: sn, sn_empty
      CHARACTER(LEN=256) :: name
      REAL(wp) :: rfac
      NAMELIST /variable/ name,sn
      NAMELIST /riverdata/ name,sn,rfac
      LOGICAL :: l_ext
      INTEGER :: num, ierr, nmlunit
      TYPE (type_input_variable),POINTER :: input_variable
      TYPE (type_river_data),POINTER :: river_data
      INTEGER :: jn

      ! Check if fabm_input.nml exists - if not, do nothing and return.
      INQUIRE( FILE='fabm_input.nml', EXIST=l_ext )
      IF (.NOT.l_ext) return

      ! Open fabm_input.nml
      CALL ctl_opn( nmlunit, 'fabm_input.nml', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, num, .FALSE. )

      ! Read any number of "variable" namelists
      DO
         ! Initialize namelist variables
         name = ''
         sn = sn_empty

         ! Read the namelist
         READ(nmlunit,nml=variable,err=98,end=99)

         ! Transfer namelist settings to new input_variable object
         ALLOCATE(input_variable, STAT=ierr)
         IF( ierr > 0 ) CALL ctl_stop( 'STOP', 'trcsms_fabm:initialize_inputs: unable to allocate input_variable object for variable '//TRIM(name) )
         input_variable%horizontal_id = fabm_get_horizontal_variable_id(model,name)
         IF (.NOT.fabm_is_variable_used(input_variable%horizontal_id)) THEN
            ! This variable was not found among FABM's horizontal variables (at least, those that are read by one or more FABM modules)
            CALL ctl_stop('STOP', 'trcsms_fabm:initialize_inputs: variable "'//TRIM(name)//'" was not found among horizontal FABM variables.')
         END IF
         ALLOCATE(input_variable%sf(1), STAT=ierr)
         IF( ierr > 0 ) CALL ctl_stop( 'STOP', 'trcsms_fabm:initialize_inputs: unable to allocate sf structure for variable '//TRIM(name) )
         CALL fld_fill(input_variable%sf, (/sn/), '', 'trcsms_fabm:initialize_inputs', 'FABM variable '//TRIM(name), 'variable' )
         ALLOCATE( input_variable%sf(1)%fnow(jpi,jpj,1)   )
         IF( sn%ln_tint ) ALLOCATE( input_variable%sf(1)%fdta(jpi,jpj,1,2) )

         ! Provide FABM with pointer to field that will receive prescribed data.
         ! NB source=data_source_user guarantees that the prescribed data takes priority over any data FABM may already have for that variable.
         CALL fabm_link_horizontal_data(model,input_variable%horizontal_id,input_variable%sf(1)%fnow(:,:,1),source=data_source_user)

         ! Prepend new input variable to list.
         input_variable%next => first_input_variable
         first_input_variable => input_variable
      END DO

98    CALL ctl_stop('STOP', 'trcsms_fabm:initialize_inputs: unable to read namelist "riverdata"')

99    REWIND(nmlunit)

      ! Read any number of "riverdata" namelists
      DO
         ! Initialize namelist variables
         name = ''
         sn = sn_empty
         rfac = 1._wp

         ! Read the namelist
         READ(nmlunit,nml=riverdata,err=198,end=199)

         ! Transfer namelist settings to new river_data object
         ALLOCATE(river_data, STAT=ierr)
         IF( ierr > 0 ) CALL ctl_stop( 'STOP', 'trcsms_fabm:initialize_inputs: unable to allocate river_data object for variable '//TRIM(name) )
         ! Check if river data name is in FABM states and 
         ! provide NEMO with position of the respective state variable 
         ! within tracer field
         DO jn=1,jp_fabm
           IF (TRIM(name) == TRIM(model%state_variables(jn)%name)) THEN
             river_data%jp_pos = jp_fabm_m1+jn
           END IF
         END DO
         IF (river_data%jp_pos == 0) THEN
           ! This variable was not found among FABM's state variables 
           ! passed to NEMO!
           CALL ctl_stop('STOP', 'trcsms_fabm:initialize_inputs: variable "'//TRIM(name)//'" was not found among FABM state variables.')
         END IF

         ALLOCATE(river_data%sf(1), STAT=ierr)
         IF( ierr > 0 ) CALL ctl_stop( 'STOP', 'trcsms_fabm:initialize_inputs: unable to allocate sf structure for variable '//TRIM(name) )
         CALL fld_fill(river_data%sf, (/sn/), '', 'trcsms_fabm:initialize_inputs', 'FABM variable '//TRIM(name), 'riverdata' )
         ALLOCATE( river_data%sf(1)%fnow(jpi,jpj,1)   )
         IF( sn%ln_tint ) ALLOCATE( river_data%sf(1)%fdta(jpi,jpj,1,2) )

         ! Load unit conversion factor:
         river_data%rn_trrnfac=rfac

         ! Prepend new input variable to list.
         river_data%next => first_river_data
         first_river_data => river_data
      END DO

198   CALL ctl_stop('STOP', 'trcsms_fabm:initialize_inputs: unable to read namelist "riverdata"')

199   RETURN

   END SUBROUTINE initialize_inputs

   SUBROUTINE update_inputs( kt )
      INTEGER, INTENT(IN) :: kt
      TYPE (type_input_variable),POINTER :: input_variable
      TYPE (type_river_data),POINTER :: river_data

      input_variable => first_input_variable
      DO WHILE (ASSOCIATED(input_variable))
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. input_variable%ntimes > 1 ) ) CALL fld_read( kt, 1, input_variable%sf )
         input_variable => input_variable%next
      END DO

      river_data => first_river_data
      DO WHILE (ASSOCIATED(river_data))
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. river_data%ntimes > 1 ) ) CALL fld_read( kt, 1, river_data%sf )
         river_data => river_data%next
      END DO

   END SUBROUTINE update_inputs

   SUBROUTINE trc_rnf_fabm( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_rnf_fabm  ***
      !!
      !! ** Purpose :   Add river loadings of biogeochemistry to states
      !!
      !! ** Action  :   tra (sms) updated with loadings at time-step kt
      !!
      !! This routines assumes river loadings to be given in
      !! state variable units * m3 / sec
      !!--------------------------------------------------------------------

      INTEGER, INTENT(in) ::   kt          ! ocean time step
      REAL(wp) :: zcoef
      INTEGER :: ji,jj,jk
      !
      TYPE (type_river_data),POINTER :: river_data

      river_data => first_river_data
      DO WHILE (ASSOCIATED(river_data))
        IF( kt == nit000 .OR. ( kt /= nit000 ) ) THEN
            DO jj = 1, jpj
              DO ji = 1, jpi
                ! convert units and divide by surface area
                ! loading / cell volume * vertical fraction of riverload
                ! dtrc / dt (river) = riverload / e1e2t / e3t * e3t * h_rnf
                !                    = riverload / e1e2t / h_rnf
                zcoef = river_data%rn_trrnfac / e1e2t(ji,jj) / h_rnf(ji,jj)
                DO jk = 1,nk_rnf(ji,jj)
                  ! Add river loadings
                  tra(ji,jj,jk,river_data%jp_pos) = tra(ji,jj,jk,river_data%jp_pos) + river_data%sf(1)%fnow(ji,jj,1)*zcoef
                END DO
              END DO
            END DO
        END IF
        river_data => river_data%next
      END DO

   END SUBROUTINE trc_rnf_fabm

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No FABM model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_fabm( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_fabm: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_fabm
#endif

   !!======================================================================
END MODULE trcsms_fabm
