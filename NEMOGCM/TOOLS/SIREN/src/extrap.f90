!----------------------------------------------------------------------
! NEMO system team, System and Interface for oceanic RElocable Nesting
!----------------------------------------------------------------------
!
! MODULE: extrap
!
! DESCRIPTION:
!> @brief 
!> This module manage extrapolation.
!>
!> @details
!>    Extrapolation method to be used is specify inside variable 
!>    strcuture, as array of string character.<br/> 
!>    - td_var\%c_extrap(1) string character is the interpolation name choose between:
!>       - 'dist_weight'
!>       - 'min_error'
!>
!>    @note Extrapolation method could be specify for each variable in namelist _namvar_,
!>    defining string character _cn\_varinfo_. By default _dist_weight_.<br/>
!>    Example:
!>       - cn_varinfo='varname1:dist_weight', 'varname2:min_error'
!>
!>    to detect point to be extrapolated:<br/>
!> @code
!>    il_detect(:,:,:)=extrap_detect(td_var, [td_level], [id_offset,] [id_rho,] [id_ext]) 
!> @endcode
!>       - il_detect(:,:,:) is 3D array of point to be extrapolated
!>       - td_var  is coarse grid variable to be extrapolated
!>       - td_level is fine grid array of level (see vgrid_get_level) [optional]
!>       - id_offset is array of offset between fine and coarse grid [optional]
!>       - id_rho    is array of refinment factor [optional]
!>       - id_ext    is array of number of points to be extrapolated [optional]
!>
!>    to extrapolate variable:<br/>
!> @code
!>    CALL extrap_fill_value( td_var, [td_level], [id_offset], [id_rho], [id_iext], [id_jext], [id_kext], [id_radius], [id_maxiter])
!> @endcode
!>       - td_var  is coarse grid variable to be extrapolated
!>       - td_level is fine grid array of level (see vgrid_get_level) [optional]
!>       - id_offset is array of offset between fine and coarse grid [optional]
!>       - id_rho    is array of refinment factor [optional]
!>       - id_iext   is number of points to be extrapolated in i-direction [optional]
!>       - id_jext   is number of points to be extrapolated in j-direction [optional]
!>       - id_kext   is number of points to be extrapolated in k-direction [optional]
!>       - id_radius is radius of the halo used to compute extrapolation [optional]
!>       - id_maxiter is maximum number of iteration [optional]
!>
!>    to add extraband to the variable (to be extrapolated):<br/>
!> @code
!>    CALL extrap_add_extrabands(td_var, [id_isize,] [id_jsize] )
!> @endcode
!>       - td_var is variable structure
!>       - id_isize : i-direction size of extra bands [optional]
!>       - id_jsize : j-direction size of extra bands [optional]
!>
!>    to delete extraband of a variable:<br/>
!> @code
!>    CALL extrap_del_extrabands(td_var, [id_isize,] [id_jsize] )
!> @endcode
!>       - td_var is variable structure
!>       - id_isize : i-direction size of extra bands [optional]
!>       - id_jsize : j-direction size of extra bands [optional]
!>
!>    to compute first derivative of 1D array:<br/>
!> @code
!>    dl_value(:)=extrap_deriv_1D( dd_value(:), dd_fill, [ld_discont] )
!> @endcode
!>       - dd_value is 1D array of variable
!>       - dd_fill is FillValue of variable
!>       - ld_discont is logical to take into account longitudinal East-West discontinuity [optional]
!>
!>    to compute first derivative of 2D array:<br/>
!> @code
!>    dl_value(:,:)=extrap_deriv_2D( dd_value(:,:), dd_fill, cd_dim, [ld_discont] )
!> @endcode
!>       - dd_value is 2D array of variable
!>       - dd_fill is FillValue of variable
!>       - cd_dim is character to compute derivative on first (I) or second (J) dimension
!>       - ld_discont is logical to take into account longitudinal East-West discontinuity [optional]
!>
!>    to compute first derivative of 3D array:<br/>
!> @code
!>    dl_value(:,:,:)=extrap_deriv_3D( dd_value(:,:,:), dd_fill, cd_dim, [ld_discont] )
!> @endcode
!>       - dd_value is 3D array of variable
!>       - dd_fill is FillValue of variable
!>       - cd_dim is character to compute derivative on first (I), second (J), or third (K) dimension
!>       - ld_discont is logical to take into account longitudinal East-West discontinuity [optional]
!>
!> @warning _FillValue must not be zero (use var_chg_FillValue())
!>
!> @author
!> J.Paul
! REVISION HISTORY:
!> @date Nov, 2013 - Initial Version
!> @date September, 2014
!> - add header
!>
!> @todo
!> - create module for each extrapolation method
!>
!> @note Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!----------------------------------------------------------------------
MODULE extrap
   USE netcdf                          ! nf90 library
   USE kind                            ! F90 kind parameter
   USE phycst                          ! physical constant
   USE global                          ! global variable
   USE fct                             ! basic useful function
   USE date                            ! date manager
   USE logger                          ! log file manager
   USE att                             ! attribute manager
   USE dim                             ! dimension manager
   USE var                             ! variable manager

   IMPLICIT NONE
   ! NOTE_avoid_public_variables_if_possible

   ! type and variable
   PRIVATE :: im_maxiter   !< default maximum number of iteration 
   PRIVATE :: im_minext    !< default minumum number of point to extrapolate
   PRIVATE :: im_mincubic  !< default minumum number of point to extrapolate for cubic interpolation

   ! function and subroutine
   PUBLIC :: extrap_detect         !< detected point to be extrapolated
   PUBLIC :: extrap_fill_value     !< extrapolate value over detected point 
   PUBLIC :: extrap_add_extrabands !< add extraband to the variable (to be extrapolated) 
   PUBLIC :: extrap_del_extrabands !< delete extraband of the variable 
   PUBLIC :: extrap_deriv_1D       !< compute first derivative of 1D array 
   PUBLIC :: extrap_deriv_2D       !< compute first derivative of 2D array 
   PUBLIC :: extrap_deriv_3D       !< compute first derivative of 3D array

   PRIVATE :: extrap__detect_wrapper      ! detected point to be extrapolated wrapper
   PRIVATE :: extrap__detect              ! detected point to be extrapolated
   PRIVATE :: extrap__fill_value_wrapper  ! extrapolate value over detected point wrapper 
   PRIVATE :: extrap__fill_value          ! extrapolate value over detected point 
   PRIVATE :: extrap__3D                  !
   PRIVATE :: extrap__3D_min_error_coef   !
   PRIVATE :: extrap__3D_min_error_fill   !
   PRIVATE :: extrap__3D_dist_weight_coef !
   PRIVATE :: extrap__3D_dist_weight_fill ! 

   INTEGER(i4), PARAMETER :: im_maxiter = 10 !< default maximum number of iteration
   INTEGER(i4), PARAMETER :: im_minext  = 2  !< default minumum number of point to extrapolate
   INTEGER(i4), PARAMETER :: im_mincubic= 4  !< default minumum number of point to extrapolate for cubic interpolation

   INTERFACE extrap_detect
      MODULE PROCEDURE extrap__detect_wrapper     !< detected point to be extrapolated
   END INTERFACE extrap_detect

   INTERFACE extrap_fill_value
      MODULE PROCEDURE extrap__fill_value_wrapper !< detected point to be interpolated
   END INTERFACE extrap_fill_value
   
CONTAINS
   !-------------------------------------------------------------------
   !> @brief
   !> This function detected point to be extrapolated, given variable structure.
   !> 
   !> @details 
   !> optionaly, you could sepcify fine grid level, refinment factor (default 1), 
   !> offset between fine and coarse grid (default compute from refinment factor
   !> as offset=(rho-1)/2), number of point to be extrapolated in each direction
   !> (default im_minext).<br/>
   !>
   !> First coarsening fine grid level, if need be, then select point near 
   !> grid point already inform.
   !>
   !> @note point to be extrapolated are selected using FillValue, 
   !> so to avoid mistake FillValue should not be zero (use var_chg_FillValue)
   !> 
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] td_var0   coarse grid variable to extrapolate
   !> @param[in] td_level1 fine grid array of level
   !> @param[in] id_offset array of offset between fine and coarse grid 
   !> @param[in] id_rho    array of refinment factor 
   !> @param[in] id_ext    array of number of points to be extrapolated
   !> @return array of point to be extrapolated
   !-------------------------------------------------------------------
   FUNCTION extrap__detect( td_var0, td_level1, &
   &                        id_offset, id_rho, id_ext )
      IMPLICIT NONE
      ! Argument
      TYPE(TVAR) ,                 INTENT(IN   ) :: td_var0
      TYPE(TVAR) , DIMENSION(:)  , INTENT(IN   ), OPTIONAL :: td_level1
      INTEGER(i4), DIMENSION(:,:), INTENT(IN   ), OPTIONAL :: id_offset
      INTEGER(i4), DIMENSION(:)  , INTENT(IN   ), OPTIONAL :: id_rho
      INTEGER(i4), DIMENSION(:)  , INTENT(IN   ), OPTIONAL :: id_ext

      ! function
      INTEGER(i4), DIMENSION(td_var0%t_dim(1)%i_len,&
      &                      td_var0%t_dim(2)%i_len,&
      &                      td_var0%t_dim(3)%i_len ) :: extrap__detect

      ! local variable
      CHARACTER(LEN=lc)                                :: cl_level

      INTEGER(i4)                                      :: il_ind
      INTEGER(i4)      , DIMENSION(:,:,:), ALLOCATABLE :: il_detect
      INTEGER(i4)      , DIMENSION(:,:,:), ALLOCATABLE :: il_tmp
      INTEGER(i4)      , DIMENSION(:,:)  , ALLOCATABLE :: il_offset
      INTEGER(i4)      , DIMENSION(:,:)  , ALLOCATABLE :: il_level1
      INTEGER(i4)      , DIMENSION(:,:)  , ALLOCATABLE :: il_level1_G0
      INTEGER(i4)      , DIMENSION(:,:)  , ALLOCATABLE :: il_extra
      INTEGER(i4)      , DIMENSION(:)    , ALLOCATABLE :: il_ext
      INTEGER(i4)      , DIMENSION(:)    , ALLOCATABLE :: il_rho
      INTEGER(i4)      , DIMENSION(:)    , ALLOCATABLE :: il_dim0

      TYPE(TVAR)                                       :: tl_var1

      ! loop indices
      INTEGER(i4) :: ji0
      INTEGER(i4) :: jj0
      INTEGER(i4) :: jk0
      INTEGER(i4) :: ji1
      INTEGER(i4) :: jj1
      INTEGER(i4) :: ji1m
      INTEGER(i4) :: jj1m
      INTEGER(i4) :: ji1p
      INTEGER(i4) :: jj1p
      !----------------------------------------------------------------

      ! init
      extrap__detect(:,:,:)=0

      ALLOCATE( il_dim0(3) )
      il_dim0(:)=td_var0%t_dim(1:3)%i_len

      ! optional argument
      ALLOCATE( il_rho(ip_maxdim) )
      il_rho(:)=1
      IF( PRESENT(id_rho) ) il_rho(1:SIZE(id_rho(:)))=id_rho(:)

      ALLOCATE( il_offset(ip_maxdim,2) )
      il_offset(:,:)=0
      IF( PRESENT(id_offset) )THEN
         il_offset(1:SIZE(id_offset(:,:),DIM=1),&
         &         1:SIZE(id_offset(:,:),DIM=2) )= id_offset(:,:)
      ELSE
         il_offset(jp_I,:)=FLOOR(REAL(il_rho(jp_I)-1,dp)*0.5)
         il_offset(jp_J,:)=FLOOR(REAL(il_rho(jp_J)-1,dp)*0.5)
      ENDIF

      ALLOCATE( il_ext(ip_maxdim) )
      il_ext(:)=im_minext
      IF( PRESENT(id_ext) ) il_ext(1:SIZE(id_ext(:)))=id_ext(:)

      ALLOCATE( il_detect(il_dim0(1),&
      &                   il_dim0(2),&
      &                   il_dim0(3)) )
      il_detect(:,:,:)=0

      ! select point already inform
      DO jk0=1,td_var0%t_dim(3)%i_len
         DO jj0=1,td_var0%t_dim(2)%i_len
            DO ji0=1,td_var0%t_dim(1)%i_len
               IF( td_var0%d_value(ji0,jj0,jk0,1) /= td_var0%d_fill ) il_detect(ji0,jj0,jk0)=1
            ENDDO
         ENDDO
      ENDDO
 
      IF( PRESENT(td_level1) )THEN
         SELECT CASE(TRIM(td_var0%c_point))
            CASE DEFAULT !'T'
               cl_level='tlevel'
            CASE('U')
               cl_level='ulevel'
            CASE('V')
               cl_level='vlevel'
            CASE('F')
               cl_level='flevel'
         END SELECT

         il_ind=var_get_index(td_level1(:),TRIM(cl_level))
         IF( il_ind == 0 )THEN
            CALL logger_error("EXTRAP DETECT: can not compute point to be "//&
            &     "extrapolated for variable "//TRIM(td_var0%c_name)//&
            &      ". can not find "//&
            &     "level for variable point "//TRIM(TRIM(td_var0%c_point)))
         ELSE
            tl_var1=var_copy(td_level1(il_ind))

            ALLOCATE( il_level1_G0( il_dim0(1), il_dim0(2)) )
            IF( ALL(tl_var1%t_dim(1:2)%i_len == il_dim0(1:2)) )THEN

               ! variable to be extrapolated use same resolution than level
               il_level1_G0(:,:)=INT(tl_var1%d_value(:,:,1,1),i4)
               
            ELSE
               ! variable to be extrapolated do not use same resolution than level
               ALLOCATE( il_level1(tl_var1%t_dim(1)%i_len, &
               &                   tl_var1%t_dim(2)%i_len) )
               ! match fine grid vertical level with coarse grid
               il_level1(:,:)=INT(tl_var1%d_value(:,:,1,1),i4)/il_rho(jp_K)

               ALLOCATE( il_extra(ip_maxdim,2) )
               ! coarsening fine grid level
               il_extra(jp_I,1)=CEILING(REAL(il_rho(jp_I)-1,dp)*0.5_dp)
               il_extra(jp_I,2)=FLOOR(REAL(il_rho(jp_I)-1,dp)*0.5_dp)

               il_extra(jp_J,1)=CEILING(REAL(il_rho(jp_J)-1,dp)*0.5_dp)
               il_extra(jp_J,2)=FLOOR(REAL(il_rho(jp_J)-1,dp)*0.5_dp)

               DO jj0=1,td_var0%t_dim(2)%i_len
                  
                  jj1=(jj0-1)*il_rho(jp_J)+1-il_offset(jp_J,1)

                  jj1m=MAX( jj1-il_extra(jp_J,1), 1 )
                  jj1p=MIN( jj1+il_extra(jp_J,2), &
                  &         tl_var1%t_dim(2)%i_len-il_offset(jp_J,2) )
                  
                  DO ji0=1,td_var0%t_dim(1)%i_len

                     ji1=(ji0-1)*il_rho(jp_I)+1-id_offset(jp_I,1)

                     ji1m=MAX( ji1-il_extra(jp_I,1), 1 )
                     ji1p=MIN( ji1+il_extra(jp_I,2), &
                     &         tl_var1%t_dim(1)%i_len-id_offset(jp_I,2) )
               
                     il_level1_G0(ji0,jj0)=MAXVAL(il_level1(ji1m:ji1p,jj1m:jj1p))

                  ENDDO
               ENDDO

               ! clean
               DEALLOCATE( il_extra )
               DEALLOCATE( il_level1 )

            ENDIF

            ! look for sea point
            DO jk0=1,td_var0%t_dim(3)%i_len
               WHERE( il_level1_G0(:,:) >= jk0)
                  il_detect(:,:,jk0)=1
               END WHERE
            ENDDO

            ! clean
            DEALLOCATE( il_level1_G0 )
            CALL var_clean(tl_var1)

         ENDIF
      ENDIF

      ! clean
      DEALLOCATE( il_offset )
 
      ALLOCATE( il_tmp(il_dim0(1),&
      &                il_dim0(2),&
      &                il_dim0(3)) )
      il_tmp(:,:,:)=il_detect(:,:,:)
      ! select extra point depending on interpolation method
      ! compute point near grid point already inform
      DO jk0=1,il_dim0(3)
         DO jj0=1,il_dim0(2)
            DO ji0=1,il_dim0(1)

               IF( il_tmp(ji0,jj0,jk0) == 1 )THEN
                  il_detect( &
                  &  MAX(1,ji0-il_ext(jp_I)):MIN(ji0+il_ext(jp_I),il_dim0(1)),&
                  &  MAX(1,jj0-il_ext(jp_J)):MIN(jj0+il_ext(jp_J),il_dim0(2)),&
                  &  MAX(1,jk0-il_ext(jp_K)):MIN(jk0+il_ext(jp_K),il_dim0(3)) &
                  &  ) = 1 
               ENDIF

            ENDDO
         ENDDO
      ENDDO
      
      ! clean
      DEALLOCATE( il_tmp )

      ! do not compute grid point already inform
      DO jk0=1,td_var0%t_dim(3)%i_len
         DO jj0=1,td_var0%t_dim(2)%i_len
            DO ji0=1,td_var0%t_dim(1)%i_len
               IF( td_var0%d_value(ji0,jj0,jk0,1) /= td_var0%d_fill ) il_detect(ji0,jj0,jk0)=0
            ENDDO
         ENDDO
      ENDDO

      ! save result
      extrap__detect(:,:,:)=il_detect(:,:,:)

      ! clean
      DEALLOCATE( il_dim0 )
      DEALLOCATE( il_ext )
      DEALLOCATE( il_detect )
      DEALLOCATE( il_rho )

   END FUNCTION extrap__detect
   !-------------------------------------------------------------------
   !> @brief
   !> This function sort variable to be extrapolated, depending on number of
   !> dimentsion, then detected point to be extrapolated.
   !> 
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !>
   !> @param[in] td_var    coarse grid variable to extrapolate
   !> @param[in] td_level  fine grid array of level
   !> @param[in] id_offset array of offset between fine and coarse grid 
   !> @param[in] id_rho    array of refinment factor 
   !> @param[in] id_ext    array of number of points to be extrapolated
   !> @return 3D array of point to be extrapolated
   !-------------------------------------------------------------------
   FUNCTION extrap__detect_wrapper( td_var, td_level, &
   &                                id_offset, id_rho, id_ext )

      IMPLICIT NONE
      ! Argument
      TYPE(TVAR) ,                 INTENT(IN   ) :: td_var
      TYPE(TVAR) , DIMENSION(:)  , INTENT(IN   ), OPTIONAL :: td_level
      INTEGER(i4), DIMENSION(:,:), INTENT(IN   ), OPTIONAL :: id_offset
      INTEGER(i4), DIMENSION(:)  , INTENT(IN   ), OPTIONAL :: id_rho
      INTEGER(i4), DIMENSION(:)  , INTENT(IN   ), OPTIONAL :: id_ext

      ! function
      INTEGER(i4), DIMENSION(td_var%t_dim(1)%i_len,&
      &                      td_var%t_dim(2)%i_len,&
      &                      td_var%t_dim(3)%i_len ) :: extrap__detect_wrapper

      ! local variable
      ! loop indices
      !----------------------------------------------------------------
      ! init
      extrap__detect_wrapper(:,:,:)=0

      IF( .NOT. ANY(td_var%t_dim(1:3)%l_use) )THEN
         ! no dimension I-J-K used
         CALL logger_debug(" EXTRAP DETECT: nothing done for variable"//&
         &              TRIM(td_var%c_name) )
      ELSE IF( ALL(td_var%t_dim(1:3)%l_use) )THEN
         
         ! detect point to be extrapolated on I-J-K
         CALL logger_debug(" EXTRAP DETECT: detect point "//&
         &              " for variable "//TRIM(td_var%c_name) )
         
         extrap__detect_wrapper(:,:,:)=extrap__detect( td_var, td_level, &
         &                                             id_offset, &
         &                                             id_rho,    &
         &                                             id_ext     )

      ELSE IF( ALL(td_var%t_dim(1:2)%l_use) )THEN

         ! detect point to be extrapolated on I-J
         CALL logger_debug(" EXTRAP DETECT: detect horizontal point "//&
         &              " for variable "//TRIM(td_var%c_name) )
         
         extrap__detect_wrapper(:,:,1:1)=extrap__detect( td_var , td_level,&
         &                                               id_offset, &
         &                                               id_rho,    &
         &                                               id_ext     )

      ELSE IF( td_var%t_dim(3)%l_use )THEN
         
         ! detect point to be extrapolated on K
         CALL logger_debug(" EXTRAP DETECT: detect vertical point "//&
         &              " for variable "//TRIM(td_var%c_name) )
         
         extrap__detect_wrapper(1:1,1:1,:)=extrap__detect( td_var , td_level, &
         &                                                 id_offset, &
         &                                                 id_rho,    &
         &                                                 id_ext     )

      ENDIF              

      CALL logger_debug(" EXTRAP DETECT: "//&
      &  TRIM(fct_str(SUM(extrap__detect_wrapper(:,:,:))))//&
      &  " points to be extrapolated" )
      
   END FUNCTION extrap__detect_wrapper
   !-------------------------------------------------------------------
   !> @brief
   !> This subroutine select method to be used for extrapolation.
   !> If need be, increase number of points to be extrapolated.
   !> Finally launch extrap__fill_value.
   !> 
   !> @details
   !> optionaly, you could specify :<br/>
   !>  - refinment factor (default 1)
   !>  - offset between fine and coarse grid (default compute from refinment factor
   !> as offset=(rho-1)/2) 
   !>  - number of point to be extrapolated in each direction (default im_minext)
   !>  - radius of the halo used to compute extrapolation
   !>  - maximum number of iteration
   !>
   !> @author J.Paul
   !> - Nov, 2013- Initial Version
   !
   !> @param[inout] td_var    variable structure
   !> @param[in] td_level     fine grid array of level
   !> @param[in] id_offset    array of offset between fine and coarse grid 
   !> @param[in] id_rho       array of refinment factor 
   !> @param[in] id_iext      number of points to be extrapolated in i-direction
   !> @param[in] id_jext      number of points to be extrapolated in j-direction
   !> @param[in] id_kext      number of points to be extrapolated in k-direction
   !> @param[in] id_radius    radius of the halo used to compute extrapolation 
   !> @param[in] id_maxiter   maximum number of iteration
   !-------------------------------------------------------------------
   SUBROUTINE extrap__fill_value_wrapper( td_var, td_level, &
   &                                      id_offset,        &
   &                                      id_rho,           &
   &                                      id_iext, id_jext, id_kext, &
   &                                      id_radius, id_maxiter )
      IMPLICIT NONE
      ! Argument
      TYPE(TVAR) ,                  INTENT(INOUT) :: td_var
      TYPE(TVAR) ,  DIMENSION(:)  , INTENT(IN   ), OPTIONAL :: td_level
      INTEGER(i4),  DIMENSION(:,:), INTENT(IN   ), OPTIONAL :: id_offset
      INTEGER(i4),  DIMENSION(:)  , INTENT(IN   ), OPTIONAL :: id_rho
      INTEGER(i4),                  INTENT(IN   ), OPTIONAL :: id_iext
      INTEGER(i4),                  INTENT(IN   ), OPTIONAL :: id_jext
      INTEGER(i4),                  INTENT(IN   ), OPTIONAL :: id_kext
      INTEGER(i4),                  INTENT(IN   ), OPTIONAL :: id_radius
      INTEGER(i4),                  INTENT(IN   ), OPTIONAL :: id_maxiter

      ! local variable
      INTEGER(i4) :: il_iext
      INTEGER(i4) :: il_jext
      INTEGER(i4) :: il_kext
      INTEGER(i4) :: il_radius
      INTEGER(i4) :: il_maxiter

      CHARACTER(LEN=lc) :: cl_method

      ! loop indices
      !----------------------------------------------------------------
      IF( .NOT. ASSOCIATED(td_var%d_value) )THEN
         CALL logger_error("EXTRAP FILL VALUE: no value "//&
         &  "associted to variable "//TRIM(td_var%c_name) )
      ELSE

         SELECT CASE(TRIM(td_var%c_extrap(1)))
            CASE('min_error')
               cl_method='min_error'
            CASE DEFAULT
               cl_method='dist_weight'

               !update variable structure
               td_var%c_extrap(1)='dist_weight'
         END SELECT

         il_iext=im_minext
         IF( PRESENT(id_iext) ) il_iext=id_iext
         il_jext=im_minext
         IF( PRESENT(id_jext) ) il_jext=id_jext
         il_kext=0
         IF( PRESENT(id_kext) ) il_kext=id_kext

         IF( TRIM(td_var%c_interp(1)) == 'cubic')THEN
            IF( il_iext > 0 .AND. il_iext < im_mincubic ) il_iext=im_mincubic
            IF( il_jext > 0 .AND. il_jext < im_mincubic ) il_jext=im_mincubic
         ENDIF

         IF( il_iext < 0 )THEN
            CALL logger_error("EXTRAP FILL VALUE: invalid "//&
            &  " number of points to be extrapolated in i-direction "//&
            &  "("//TRIM(fct_str(il_iext))//")")
         ENDIF

         IF( il_jext < 0 )THEN
            CALL logger_error("EXTRAP FILL VALUE: invalid "//&
            &  " number of points to be extrapolated in j-direction "//&
            &  "("//TRIM(fct_str(il_jext))//")")
         ENDIF

         IF( il_kext < 0 )THEN
            CALL logger_error("EXTRAP FILL VALUE: invalid "//&
            &  " number of points to be extrapolated in k-direction "//&
            &  "("//TRIM(fct_str(il_kext))//")")
         ENDIF

         IF( (il_iext /= 0 .AND. td_var%t_dim(1)%l_use) .OR. &
         &   (il_jext /= 0 .AND. td_var%t_dim(2)%l_use) .OR. &
         &   (il_kext /= 0 .AND. td_var%t_dim(3)%l_use) )THEN

            ! number of point use to compute box
            il_radius=1
            IF( PRESENT(id_radius) ) il_radius=id_radius
            IF( il_radius < 0 )THEN
               CALL logger_error("EXTRAP FILL VALUE: invalid "//&
               &  " radius of the box used to compute extrapolation "//&
               &  "("//TRIM(fct_str(il_radius))//")")
            ENDIF

            ! maximum number of iteration
            il_maxiter=im_maxiter
            IF( PRESENT(id_maxiter) ) il_maxiter=id_maxiter
            IF( il_maxiter < 0 )THEN
               CALL logger_error("EXTRAP FILL VALUE: invalid "//&
               &  " maximum nuber of iteration "//&
               &  "("//TRIM(fct_str(il_maxiter))//")")
            ENDIF

            CALL logger_info("EXTRAP FILL: extrapolate "//TRIM(td_var%c_name)//&
            &  " using "//TRIM(cl_method)//" method." )

            CALL extrap__fill_value( td_var, cl_method, &
            &                        il_iext, il_jext, il_kext,   &
            &                        il_radius, il_maxiter,       &
            &                        td_level,                    &
            &                        id_offset, id_rho )
 
         ENDIF
 
      ENDIF

   END SUBROUTINE extrap__fill_value_wrapper
   !-------------------------------------------------------------------
   !> @brief
   !> This subroutine compute point to be extrapolated, then extrapolate point.
   !> 
   !> @details 
   !> optionaly, you could specify :<br/>
   !>  - refinment factor (default 1)
   !>  - offset between fine and coarse grid (default compute from refinment factor
   !> as offset=(rho-1)/2) 
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[inout] td_var    variable structure
   !> @param[in] cd_method    extrapolation method
   !> @param[in] id_iext      number of points to be extrapolated in i-direction
   !> @param[in] id_jext      number of points to be extrapolated in j-direction
   !> @param[in] id_kext      number of points to be extrapolated in k-direction
   !> @param[in] id_radius    radius of the halo used to compute extrapolation
   !> @param[in] id_maxiter   maximum number of iteration
   !> @param[in] td_level     fine grid array of level
   !> @param[in] id_offset    array of offset between fine and coarse grid 
   !> @param[in] id_rho       array of refinment factor
   !-------------------------------------------------------------------
   SUBROUTINE extrap__fill_value( td_var, cd_method, &
   &                              id_iext, id_jext, id_kext, &
   &                              id_radius, id_maxiter, &
   &                              td_level,          &
   &                              id_offset,         &
   &                              id_rho )
      IMPLICIT NONE
      ! Argument
      TYPE(TVAR)      ,                 INTENT(INOUT) :: td_var
      CHARACTER(LEN=*),                 INTENT(IN   ) :: cd_method
      INTEGER(i4)     ,                 INTENT(IN   ) :: id_iext
      INTEGER(i4)     ,                 INTENT(IN   ) :: id_jext
      INTEGER(i4)     ,                 INTENT(IN   ) :: id_kext
      INTEGER(i4)     ,                 INTENT(IN   ) :: id_radius
      INTEGER(i4)     ,                 INTENT(IN   ) :: id_maxiter
      TYPE(TVAR)      , DIMENSION(:)  , INTENT(IN   ), OPTIONAL :: td_level
      INTEGER(i4)     , DIMENSION(:,:), INTENT(IN   ), OPTIONAL :: id_offset
      INTEGER(i4)     , DIMENSION(:)  , INTENT(IN   ), OPTIONAL :: id_rho

      ! local variable
      CHARACTER(LEN=lc)                            :: cl_extrap

      INTEGER(i4), DIMENSION(:,:,:)  , ALLOCATABLE :: il_detect

      TYPE(TATT)                                   :: tl_att

      ! loop indices
      !----------------------------------------------------------------

      !1- detect point to be extrapolated
      ALLOCATE( il_detect( td_var%t_dim(1)%i_len, &
      &                    td_var%t_dim(2)%i_len, &
      &                    td_var%t_dim(3)%i_len) )

      il_detect(:,:,:) = extrap_detect( td_var, td_level, &
      &                                 id_offset,        &
      &                                 id_rho,           &
      &                                 id_ext=(/id_iext, id_jext, id_kext/) )
      !2- add attribute to variable
      cl_extrap=fct_concat(td_var%c_extrap(:))
      tl_att=att_init('extrapolation',cl_extrap)
      CALL var_move_att(td_var, tl_att)

      CALL att_clean(tl_att)

      CALL logger_info(" EXTRAP FILL: "//&
         &              TRIM(fct_str(SUM(il_detect(:,:,:))))//&
         &              " point(s) to extrapolate " )

      !3- extrapolate
      CALL extrap__3D(td_var%d_value(:,:,:,:), td_var%d_fill,    &
      &               il_detect(:,:,:),                           &
      &               cd_method, id_radius, id_maxiter  )

      DEALLOCATE(il_detect)

   END SUBROUTINE extrap__fill_value
   !-------------------------------------------------------------------
   !> @brief
   !> This subroutine compute point to be extrapolated in 3D array.
   !> 
   !> @details 
   !> in case of 'min_error' method:<br/>
   !>    - compute derivative in i-, j- and k- direction
   !>    - compute minimum error coefficient (distance to center of halo)
   !>    - compute extrapolatd values by calculated minimum error using taylor expansion 
   !> in case of 'dist_weight' method:<br/>
   !>    - compute distance weight coefficient (inverse of distance to center of halo) 
   !>    - compute extrapolatd values using Inverse Distance Weighting
   !>
   !> @author J.Paul
   !> - Nov, 2013- Initial Version
   !
   !> @param[inout] dd_value  3D array of variable to be extrapolated
   !> @param[in] dd_fill      FillValue of variable
   !> @param[inout] id_detect array of point to extrapolate
   !> @param[in] cd_method    extrapolation method
   !> @param[in] id_radius    radius of the halo used to compute extrapolation
   !-------------------------------------------------------------------
   SUBROUTINE extrap__3D( dd_value, dd_fill, id_detect,&
   &                      cd_method, id_radius, id_maxiter )
      IMPLICIT NONE
      ! Argument
      REAL(dp)   , DIMENSION(:,:,:,:), INTENT(INOUT) :: dd_value
      REAL(dp)   ,                   INTENT(IN   ) :: dd_fill
      INTEGER(i4), DIMENSION(:,:,:), INTENT(INOUT) :: id_detect
      CHARACTER(LEN=*),              INTENT(IN   ) :: cd_method
      INTEGER(i4),                   INTENT(IN   ) :: id_radius
      INTEGER(i4),                   INTENT(IN   ) :: id_maxiter

      ! local variable
      INTEGER(i4) :: il_imin
      INTEGER(i4) :: il_imax
      INTEGER(i4) :: il_jmin
      INTEGER(i4) :: il_jmax
      INTEGER(i4) :: il_kmin
      INTEGER(i4) :: il_kmax
      INTEGER(i4) :: il_iter
      INTEGER(i4) :: il_radius

      INTEGER(i4), DIMENSION(4) :: il_shape
      INTEGER(i4), DIMENSION(3) :: il_dim

      INTEGER(i4), DIMENSION(:,:,:), ALLOCATABLE :: il_detect

      REAL(dp)   , DIMENSION(:,:,:), ALLOCATABLE :: dl_dfdx
      REAL(dp)   , DIMENSION(:,:,:), ALLOCATABLE :: dl_dfdy
      REAL(dp)   , DIMENSION(:,:,:), ALLOCATABLE :: dl_dfdz
      REAL(dp)   , DIMENSION(:,:,:), ALLOCATABLE :: dl_coef

      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: jj
      INTEGER(i4) :: jk
      INTEGER(i4) :: jl
      !----------------------------------------------------------------

      il_shape(:)=SHAPE(dd_value)

      ALLOCATE( il_detect( il_shape(1), il_shape(2), il_shape(3)) )

      SELECT CASE(TRIM(cd_method))
      CASE('min_error')
         DO jl=1,il_shape(4)

            ! intitialise table of poitn to be extrapolated
            il_detect(:,:,:)=id_detect(:,:,:)

            il_iter=1
            DO WHILE( ANY(il_detect(:,:,:)==1) )
               ! change extend value to minimize number of iteration
               il_radius=id_radius+(il_iter/id_maxiter)

               ALLOCATE( dl_dfdx(il_shape(1), il_shape(2), il_shape(3)) ) 
               ALLOCATE( dl_dfdy(il_shape(1), il_shape(2), il_shape(3)) ) 
               ALLOCATE( dl_dfdz(il_shape(1), il_shape(2), il_shape(3)) ) 

               ! compute derivative in i-direction
               dl_dfdx(:,:,:)=dd_fill
               IF( il_shape(1) > 1 )THEN
                  dl_dfdx(:,:,:)=extrap_deriv_3D( dd_value(:,:,:,jl), dd_fill, 'I' )
               ENDIF

               ! compute derivative in j-direction
               dl_dfdy(:,:,:)=dd_fill
               IF( il_shape(2) > 1 )THEN
                  dl_dfdy(:,:,:)=extrap_deriv_3D( dd_value(:,:,:,jl), dd_fill, 'J' )
               ENDIF
 
               ! compute derivative in k-direction
               dl_dfdz(:,:,:)=dd_fill
               IF( il_shape(3) > 1 )THEN
                  dl_dfdz(:,:,:)=extrap_deriv_3D( dd_value(:,:,:,jl), dd_fill, 'K' )
               ENDIF
 
               il_dim(1)=2*il_radius+1
               IF( il_shape(1) < 2*il_radius+1 ) il_dim(1)=1
               il_dim(2)=2*il_radius+1
               IF( il_shape(2) < 2*il_radius+1 ) il_dim(2)=1
               il_dim(3)=2*il_radius+1
               IF( il_shape(3) < 2*il_radius+1 ) il_dim(3)=1

               ALLOCATE( dl_coef(il_dim(1), il_dim(2), il_dim(3)) ) 

               dl_coef(:,:,:)=extrap__3D_min_error_coef(dd_value( 1:il_dim(1), &
               &                                                  1:il_dim(2), &
               &                                                  1:il_dim(3), &
               &                                                  jl ))

               DO jk=1,il_shape(3)
                  IF( ALL(il_detect(:,:,jk) == 0) ) CYCLE
                  DO jj=1,il_shape(2)
                     IF( ALL(il_detect(:,jj,jk) == 0) ) CYCLE
                     DO ji=1,il_shape(1)

                        IF( il_detect(ji,jj,jk) == 1 )THEN
                          
                           il_imin=MAX(ji-il_radius,1)
                           il_imax=MIN(ji+il_radius,il_shape(1))
                           IF( il_dim(1) == 1 )THEN
                              il_imin=ji
                              il_imax=ji
                           ENDIF

                           il_jmin=MAX(jj-il_radius,1)
                           il_jmax=MIN(jj+il_radius,il_shape(2))
                           IF( il_dim(2) == 1 )THEN
                              il_jmin=jj
                              il_jmax=jj
                           ENDIF

                           il_kmin=MAX(jk-il_radius,1)
                           il_kmax=MIN(jk+il_radius,il_shape(3))
                           IF( il_dim(3) == 1 )THEN
                              il_kmin=jk
                              il_kmax=jk
                           ENDIF

                           dd_value(ji,jj,jk,jl)=extrap__3D_min_error_fill( &
                           &  dd_value( il_imin:il_imax, &
                           &            il_jmin:il_jmax, &
                           &            il_kmin:il_kmax,jl ), dd_fill, il_radius, &
                           &  dl_dfdx(  il_imin:il_imax, &
                           &            il_jmin:il_jmax, &
                           &            il_kmin:il_kmax ), &
                           &  dl_dfdy(  il_imin:il_imax, &
                           &            il_jmin:il_jmax, &
                           &            il_kmin:il_kmax ), &
                           &  dl_dfdz(  il_imin:il_imax, &
                           &            il_jmin:il_jmax, &
                           &            il_kmin:il_kmax ), &
                           &  dl_coef(:,:,:) )

                           IF( dd_value(ji,jj,jk,jl) /= dd_fill )THEN
                              il_detect(ji,jj,jk)= 0
                           ENDIF

                        ENDIF

                     ENDDO
                  ENDDO
               ENDDO

               DEALLOCATE( dl_dfdx )
               DEALLOCATE( dl_dfdy )
               DEALLOCATE( dl_dfdz )
               DEALLOCATE( dl_coef )

               il_iter=il_iter+1
            ENDDO
         ENDDO

      CASE DEFAULT ! 'dist_weight'
         DO jl=1,il_shape(4)

            ! intitialise table of poitn to be extrapolated
            il_detect(:,:,:)=id_detect(:,:,:)

            il_iter=1
            DO WHILE( ANY(il_detect(:,:,:)==1) )
               ! change extend value to minimize number of iteration
               il_radius=id_radius+(il_iter/id_maxiter)

               il_dim(1)=2*il_radius+1
               IF( il_shape(1) < 2*il_radius+1 ) il_dim(1)=1
               il_dim(2)=2*il_radius+1
               IF( il_shape(2) < 2*il_radius+1 ) il_dim(2)=1
               il_dim(3)=2*il_radius+1
               IF( il_shape(3) < 2*il_radius+1 ) il_dim(3)=1
               
               ALLOCATE( dl_coef(il_dim(1), il_dim(2), il_dim(3)) )

               dl_coef(:,:,:)=extrap__3D_dist_weight_coef(dd_value(1:il_dim(1), &
               &                                                   1:il_dim(2), &
               &                                                   1:il_dim(3), &
               &                                                   jl ) )

               DO jk=1,il_shape(3)
                  IF( ALL(il_detect(:,:,jk) == 0) ) CYCLE
                  DO jj=1,il_shape(2)
                     IF( ALL(il_detect(:,jj,jk) == 0) ) CYCLE
                     DO ji=1,il_shape(1)

                        IF( il_detect(ji,jj,jk) == 1 )THEN
                           
                           il_imin=MAX(ji-il_radius,1)
                           il_imax=MIN(ji+il_radius,il_shape(1))
                           IF( il_dim(1) == 1 )THEN
                              il_imin=ji
                              il_imax=ji
                           ENDIF

                           il_jmin=MAX(jj-il_radius,1)
                           il_jmax=MIN(jj+il_radius,il_shape(2))
                           IF( il_dim(2) == 1 )THEN
                              il_jmin=jj
                              il_jmax=jj
                           ENDIF

                           il_kmin=MAX(jk-il_radius,1)
                           il_kmax=MIN(jk+il_radius,il_shape(3))
                           IF( il_dim(3) == 1 )THEN
                              il_kmin=jk
                              il_kmax=jk
                           ENDIF

                           dd_value(ji,jj,jk,jl)=extrap__3D_dist_weight_fill( &
                           &  dd_value( il_imin:il_imax, &
                           &            il_jmin:il_jmax, &
                           &            il_kmin:il_kmax, &
                           &            jl), dd_fill, il_radius, &
                           &  dl_coef(:,:,:) )

                           IF( dd_value(ji,jj,jk,jl) /= dd_fill )THEN
                              il_detect(ji,jj,jk)= 0
                           ENDIF

                        ENDIF

                     ENDDO
                  ENDDO
               ENDDO

               DEALLOCATE( dl_coef )
               il_iter=il_iter+1
            ENDDO
         ENDDO            
      END SELECT

      DEALLOCATE( il_detect )

   END SUBROUTINE extrap__3D
   !-------------------------------------------------------------------
   !> @brief
   !> This function compute derivative of 1D array.
   !> 
   !> @details 
   !> optionaly you could specify to take into account east west discontinuity
   !> (-180° 180° or 0° 360° for longitude variable)
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] dd_value     1D array of variable to be extrapolated
   !> @param[in] dd_fill      FillValue of variable
   !> @param[in] ld_discont   logical to take into account east west discontinuity 
   !-------------------------------------------------------------------
   PURE FUNCTION extrap_deriv_1D( dd_value, dd_fill, ld_discont )

      IMPLICIT NONE
      ! Argument
      REAL(dp)   , DIMENSION(:), INTENT(IN) :: dd_value
      REAL(dp)                 , INTENT(IN) :: dd_fill
      LOGICAL                  , INTENT(IN), OPTIONAL :: ld_discont

      ! function
      REAL(dp), DIMENSION(SIZE(dd_value,DIM=1) ) :: extrap_deriv_1D

      ! local variable
      INTEGER(i4)                            :: il_imin
      INTEGER(i4)                            :: il_imax
      INTEGER(i4), DIMENSION(1)              :: il_shape

      REAL(dp)                               :: dl_min
      REAL(dp)                               :: dl_max
      REAL(dp)   , DIMENSION(:), ALLOCATABLE :: dl_value

      LOGICAL                                :: ll_discont

      ! loop indices
      INTEGER(i4) :: ji

      INTEGER(i4) :: i1
      INTEGER(i4) :: i2
      !----------------------------------------------------------------
      ! init
      extrap_deriv_1D(:)=dd_fill

      ll_discont=.FALSE.
      IF( PRESENT(ld_discont) ) ll_discont=ld_discont

      il_shape(:)=SHAPE(dd_value(:))

      ALLOCATE( dl_value(3))

      ! compute derivative in i-direction
      DO ji=1,il_shape(1)
         
            il_imin=MAX(ji-1,1)
            il_imax=MIN(ji+1,il_shape(1))

            IF( il_imin==ji-1 .AND. il_imax==ji+1 )THEN
               i1=1  ; i2=3
            ELSEIF( il_imin==ji .AND. il_imax==ji+1 )THEN
               i1=1  ; i2=2
            ELSEIF( il_imin==ji-1 .AND. il_imax==ji )THEN
               i1=2  ; i2=3
            ENDIF

            dl_value(i1:i2)=dd_value(il_imin:il_imax)
            IF( il_imin == 1 )THEN
               dl_value(:)=EOSHIFT( dl_value(:), &
               &                    DIM=1,         &
               &                    SHIFT=-1,      &
               &                    BOUNDARY=dl_value(1) )
            ENDIF
            IF( il_imax == il_shape(1) )THEN
               dl_value(:)=EOSHIFT( dl_value(:), &
               &                    DIM=1,         &
               &                    SHIFT=1,       &
               &                    BOUNDARY=dl_value(3))
            ENDIF

            IF( ll_discont )THEN
               dl_min=MINVAL( dl_value(:), dl_value(:)/=dd_fill )
               dl_max=MAXVAL( dl_value(:), dl_value(:)/=dd_fill )
               IF( dl_min < -170_dp .AND. dl_max > 170_dp )THEN
                  WHERE( dl_value(:) < 0._dp ) 
                     dl_value(:) = dl_value(:)+360._dp
                  END WHERE
               ELSEIF( dl_min < 10_dp .AND. dl_max > 350_dp )THEN
                  WHERE( dl_value(:) > 180._dp ) 
                     dl_value(:) = dl_value(:)-180._dp
                  END WHERE
               ENDIF
            ENDIF

         IF( dl_value( 2) /= dd_fill .AND. & ! ji
         &   dl_value( 3) /= dd_fill .AND. & ! ji+1
         &   dl_value( 1) /= dd_fill )THEN   ! ji-1

            extrap_deriv_1D(ji)=&
            &  ( dl_value(3) - dl_value(1) ) / &
            &  REAL( il_imax-il_imin ,dp)

         ENDIF

      ENDDO

      DEALLOCATE( dl_value )

   END FUNCTION extrap_deriv_1D
   !-------------------------------------------------------------------
   !> @brief
   !> This function compute derivative of 2D array.
   !> you have to specify in which direction derivative have to be computed:
   !> first (I) or second (J) dimension. 
   !>
   !> @details 
   !> optionaly you could specify to take into account east west discontinuity
   !> (-180° 180° or 0° 360° for longitude variable)
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] dd_value     2D array of variable to be extrapolated
   !> @param[in] dd_fill      FillValue of variable
   !> @param[in] cd_dim       compute derivative on first (I) or second (J) dimension 
   !> @param[in] ld_discont   logical to take into account east west discontinuity 
   !-------------------------------------------------------------------
   FUNCTION extrap_deriv_2D( dd_value, dd_fill, cd_dim, ld_discont )

      IMPLICIT NONE
      ! Argument
      REAL(dp)   , DIMENSION(:,:), INTENT(IN) :: dd_value
      REAL(dp)                   , INTENT(IN) :: dd_fill
      CHARACTER(LEN=*)           , INTENT(IN) :: cd_dim
      LOGICAL                    , INTENT(IN), OPTIONAL :: ld_discont

      ! function
      REAL(dp), DIMENSION(SIZE(dd_value,DIM=1), &
      &                   SIZE(dd_value,DIM=2) ) :: extrap_deriv_2D

      ! local variable
      INTEGER(i4)                              :: il_imin
      INTEGER(i4)                              :: il_imax
      INTEGER(i4)                              :: il_jmin
      INTEGER(i4)                              :: il_jmax
      INTEGER(i4), DIMENSION(2)                :: il_shape

      REAL(dp)                                 :: dl_min
      REAL(dp)                                 :: dl_max
      REAL(dp)   , DIMENSION(:,:), ALLOCATABLE :: dl_value

      LOGICAL                                  :: ll_discont

      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: jj

      INTEGER(i4) :: i1
      INTEGER(i4) :: i2

      INTEGER(i4) :: j1
      INTEGER(i4) :: j2
      !----------------------------------------------------------------
      ! init
      extrap_deriv_2D(:,:)=dd_fill

      ll_discont=.FALSE.
      IF( PRESENT(ld_discont) ) ll_discont=ld_discont

      il_shape(:)=SHAPE(dd_value(:,:))

      SELECT CASE(TRIM(fct_upper(cd_dim)))

      CASE('I')

         ALLOCATE( dl_value(3,il_shape(2)) )
         ! compute derivative in i-direction
         DO ji=1,il_shape(1)

            ! init
            dl_value(:,:)=dd_fill
            
            il_imin=MAX(ji-1,1)
            il_imax=MIN(ji+1,il_shape(1))

            IF( il_imin==ji-1 .AND. il_imax==ji+1 )THEN
               i1=1  ; i2=3
            ELSEIF( il_imin==ji .AND. il_imax==ji+1 )THEN
               i1=1  ; i2=2
            ELSEIF( il_imin==ji-1 .AND. il_imax==ji )THEN
               i1=2  ; i2=3
            ENDIF

            dl_value(i1:i2,:)=dd_value(il_imin:il_imax,:)
            IF( il_imin == 1 )THEN
               dl_value(:,:)=EOSHIFT( dl_value(:,:), &
               &                      DIM=1,         &
               &                      SHIFT=-1,      &
               &                      BOUNDARY=dl_value(1,:) )
            ENDIF
            IF( il_imax == il_shape(1) )THEN
               dl_value(:,:)=EOSHIFT( dl_value(:,:), &
               &                      DIM=1,         &
               &                      SHIFT=1,       &
               &                      BOUNDARY=dl_value(3,:))
            ENDIF

            IF( ll_discont )THEN
               dl_min=MINVAL( dl_value(:,:), dl_value(:,:)/=dd_fill )
               dl_max=MAXVAL( dl_value(:,:), dl_value(:,:)/=dd_fill )
               IF( dl_min < -170_dp .AND. dl_max > 170_dp )THEN
                  WHERE( dl_value(:,:) < 0_dp ) 
                     dl_value(:,:) = dl_value(:,:)+360._dp
                  END WHERE
               ELSEIF( dl_min < 10_dp .AND. dl_max > 350_dp )THEN
                  WHERE( dl_value(:,:) > 180 ) 
                     dl_value(:,:) = dl_value(:,:)-180._dp
                  END WHERE
               ENDIF
            ENDIF
            
            WHERE( dl_value(2,:) /= dd_fill .AND. &  ! ji
            &      dl_value(3,:) /= dd_fill .AND. &  ! ji+1
            &      dl_value(1,:) /= dd_fill )        ! ji-1

               extrap_deriv_2D(ji,:)=&
               &  ( dl_value(3,:) - dl_value(1,:) ) / &
               &    REAL( il_imax-il_imin,dp)

            END WHERE

         ENDDO

      CASE('J')

         ALLOCATE( dl_value(il_shape(1),3) )
         ! compute derivative in j-direction
         DO jj=1,il_shape(2)
         
            il_jmin=MAX(jj-1,1)
            il_jmax=MIN(jj+1,il_shape(2))

            IF( il_jmin==jj-1 .AND. il_jmax==jj+1 )THEN
               j1=1  ; j2=3
            ELSEIF( il_jmin==jj .AND. il_jmax==jj+1 )THEN
               j1=1  ; j2=2
            ELSEIF( il_jmin==jj-1 .AND. il_jmax==jj )THEN
               j1=2  ; j2=3
            ENDIF

            dl_value(:,j1:j2)=dd_value(:,il_jmin:il_jmax)
            IF( il_jmin == 1 )THEN
               dl_value(:,:)=EOSHIFT( dl_value(:,:), &
               &                      DIM=2,         &
               &                      SHIFT=-1,      &
               &                      BOUNDARY=dl_value(:,1))
            ENDIF
            IF( il_jmax == il_shape(2) )THEN
               dl_value(:,:)=EOSHIFT( dl_value(:,:), &
               &                      DIM=2,         &
               &                      SHIFT=1,       &
               &                      BOUNDARY=dl_value(:,3))
            ENDIF

            IF( ll_discont )THEN
               dl_min=MINVAL( dl_value(:,:), dl_value(:,:)/=dd_fill )
               dl_max=MAXVAL( dl_value(:,:), dl_value(:,:)/=dd_fill )
               IF( dl_min < -170_dp .AND. dl_max > 170_dp )THEN
                  WHERE( dl_value(:,:) < 0_dp ) 
                     dl_value(:,:) = dl_value(:,:)+360._dp
                  END WHERE
               ELSEIF( dl_min < 10_dp .AND. dl_max > 350_dp )THEN
                  WHERE( dl_value(:,:) > 180 ) 
                     dl_value(:,:) = dl_value(:,:)-180._dp
                  END WHERE
               ENDIF
            ENDIF

            WHERE( dl_value(:, 2) /= dd_fill .AND. & ! jj
            &      dl_value(:, 3) /= dd_fill .AND. & ! jj+1
            &      dl_value(:, 1) /= dd_fill )       ! jj-1

               extrap_deriv_2D(:,jj)=&
               &  ( dl_value(:,3) - dl_value(:,1) ) / &
               &   REAL(il_jmax-il_jmin,dp)         

            END WHERE

         ENDDO
         
      END SELECT

      DEALLOCATE( dl_value )

   END FUNCTION extrap_deriv_2D
   !-------------------------------------------------------------------
   !> @brief
   !> This function compute derivative of 3D array.
   !> you have to specify in which direction derivative have to be computed:
   !> first (I), second (J) or third (K) dimension.
   !> 
   !> @details 
   !> optionaly you could specify to take into account east west discontinuity
   !> (-180° 180° or 0° 360° for longitude variable)
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[inout] dd_value  3D array of variable to be extrapolated
   !> @param[in] dd_fill      FillValue of variable
   !> @param[in] cd_dim       compute derivative on first (I) second (J) or third (K) dimension   
   !> @param[in] ld_discont   logical to take into account east west discontinuity
   !-------------------------------------------------------------------
   PURE FUNCTION extrap_deriv_3D( dd_value, dd_fill, cd_dim, ld_discont )

      IMPLICIT NONE
      ! Argument
      REAL(dp)        , DIMENSION(:,:,:), INTENT(IN) :: dd_value
      REAL(dp)                          , INTENT(IN) :: dd_fill
      CHARACTER(LEN=*)                  , INTENT(IN) :: cd_dim
      LOGICAL                           , INTENT(IN), OPTIONAL :: ld_discont

      ! function
      REAL(dp), DIMENSION(SIZE(dd_value,DIM=1), &
      &                   SIZE(dd_value,DIM=2), &
      &                   SIZE(dd_value,DIM=3)) :: extrap_deriv_3D

      ! local variable
      INTEGER(i4)                                :: il_imin
      INTEGER(i4)                                :: il_imax
      INTEGER(i4)                                :: il_jmin
      INTEGER(i4)                                :: il_jmax
      INTEGER(i4)                                :: il_kmin
      INTEGER(i4)                                :: il_kmax
      INTEGER(i4), DIMENSION(3)                  :: il_shape

      REAL(dp)                                   :: dl_min
      REAL(dp)                                   :: dl_max
      REAL(dp)   , DIMENSION(:,:,:), ALLOCATABLE :: dl_value

      LOGICAL                                    :: ll_discont

      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: jj
      INTEGER(i4) :: jk

      INTEGER(i4) :: i1
      INTEGER(i4) :: i2

      INTEGER(i4) :: j1
      INTEGER(i4) :: j2
      
      INTEGER(i4) :: k1
      INTEGER(i4) :: k2      
      !----------------------------------------------------------------
      ! init
      extrap_deriv_3D(:,:,:)=dd_fill

      ll_discont=.FALSE.
      IF( PRESENT(ld_discont) ) ll_discont=ld_discont

      il_shape(:)=SHAPE(dd_value(:,:,:))


      SELECT CASE(TRIM(fct_upper(cd_dim)))

      CASE('I')

         ALLOCATE( dl_value(3,il_shape(2),il_shape(3)) )
         ! compute derivative in i-direction
         DO ji=1,il_shape(1)
            
            il_imin=MAX(ji-1,1)
            il_imax=MIN(ji+1,il_shape(1))

            IF( il_imin==ji-1 .AND. il_imax==ji+1 )THEN
               i1=1  ; i2=3
            ELSEIF( il_imin==ji .AND. il_imax==ji+1 )THEN
               i1=1  ; i2=2
            ELSEIF( il_imin==ji-1 .AND. il_imax==ji )THEN
               i1=2  ; i2=3
            ENDIF

            dl_value(i1:i2,:,:)=dd_value(il_imin:il_imax,:,:)
            IF( il_imin == 1 )THEN
               dl_value(:,:,:)=EOSHIFT( dl_value(:,:,:), &
               &                      DIM=1,         &
               &                      SHIFT=-1,      &
               &                      BOUNDARY=dl_value(1,:,:) )
            ENDIF
            IF( il_imax == il_shape(1) )THEN
               dl_value(:,:,:)=EOSHIFT( dl_value(:,:,:), &
               &                      DIM=1,         &
               &                      SHIFT=1,       &
               &                      BOUNDARY=dl_value(3,:,:))
            ENDIF

            IF( ll_discont )THEN
               dl_min=MINVAL( dl_value(:,:,:), dl_value(:,:,:)/=dd_fill )
               dl_max=MAXVAL( dl_value(:,:,:), dl_value(:,:,:)/=dd_fill )
               IF( dl_min < -170_dp .AND. dl_max > 170_dp )THEN
                  WHERE( dl_value(:,:,:) < 0_dp ) 
                     dl_value(:,:,:) = dl_value(:,:,:)+360._dp
                  END WHERE
               ELSEIF( dl_min < 10_dp .AND. dl_max > 350_dp )THEN
                  WHERE( dl_value(:,:,:) > 180 ) 
                     dl_value(:,:,:) = dl_value(:,:,:)-180._dp
                  END WHERE
               ENDIF
            ENDIF

            WHERE( dl_value(2,:,:) /= dd_fill .AND. & ! ji
            &      dl_value(3,:,:) /= dd_fill .AND. & !ji+1 
            &      dl_value(1,:,:) /= dd_fill )       !ji-1

               extrap_deriv_3D(ji,:,:)= &
               &  ( dl_value(3,:,:) - dl_value(1,:,:) ) / &
               &  REAL( il_imax-il_imin ,dp)

            END WHERE

         ENDDO

      CASE('J')

         ALLOCATE( dl_value(il_shape(1),3,il_shape(3)) )
         ! compute derivative in j-direction
         DO jj=1,il_shape(2)
         
            il_jmin=MAX(jj-1,1)
            il_jmax=MIN(jj+1,il_shape(2))

            IF( il_jmin==jj-1 .AND. il_jmax==jj+1 )THEN
               j1=1  ; j2=3
            ELSEIF( il_jmin==jj .AND. il_jmax==jj+1 )THEN
               j1=1  ; j2=2
            ELSEIF( il_jmin==jj-1 .AND. il_jmax==jj )THEN
               j1=2  ; j2=3
            ENDIF

            dl_value(:,j1:j2,:)=dd_value(:,il_jmin:il_jmax,:)
            IF( il_jmin == 1 )THEN
               dl_value(:,:,:)=EOSHIFT( dl_value(:,:,:), &
               &                      DIM=2,         &
               &                      SHIFT=-1,      &
               &                      BOUNDARY=dl_value(:,1,:) )
            ENDIF
            IF( il_jmax == il_shape(2) )THEN
               dl_value(:,:,:)=EOSHIFT( dl_value(:,:,:), &
               &                      DIM=2,         &
               &                      SHIFT=1,       &
               &                      BOUNDARY=dl_value(:,3,:))
            ENDIF

            IF( ll_discont )THEN
               dl_min=MINVAL( dl_value(:,:,:), dl_value(:,:,:)/=dd_fill )
               dl_max=MAXVAL( dl_value(:,:,:), dl_value(:,:,:)/=dd_fill )
               IF( dl_min < -170_dp .AND. dl_max > 170_dp )THEN
                  WHERE( dl_value(:,:,:) < 0_dp ) 
                     dl_value(:,:,:) = dl_value(:,:,:)+360._dp
                  END WHERE
               ELSEIF( dl_min < 10_dp .AND. dl_max > 350_dp )THEN
                  WHERE( dl_value(:,:,:) > 180 ) 
                     dl_value(:,:,:) = dl_value(:,:,:)-180._dp
                  END WHERE
               ENDIF
            ENDIF

            WHERE( dl_value(:, 2,:) /= dd_fill .AND. & ! jj
               &   dl_value(:, 3,:) /= dd_fill .AND. & ! jj+1
            &      dl_value(:, 1,:) /= dd_fill )       ! jj-1

               extrap_deriv_3D(:,jj,:)=&
               &  ( dl_value(:,3,:) - dl_value(:,1,:) ) / &
               &  REAL( il_jmax - il_jmin ,dp)         

            END WHERE

         ENDDO
         
      CASE('K')
         ! compute derivative in k-direction
         DO jk=1,il_shape(3)

            il_kmin=MAX(jk-1,1)
            il_kmax=MIN(jk+1,il_shape(3))

            IF( il_kmin==jk-1 .AND. il_kmax==jk+1 )THEN
               k1=1  ; k2=3
            ELSEIF( il_kmin==jk .AND. il_kmax==jk+1 )THEN
               k1=1  ; k2=2
            ELSEIF( il_kmin==jk-1 .AND. il_kmax==jk )THEN
               k1=2  ; k2=3
            ENDIF

            dl_value(:,:,k1:k2)=dd_value(:,:,il_kmin:il_kmax)
            IF( il_kmin == 1 )THEN
               dl_value(:,:,:)=EOSHIFT( dl_value(:,:,:), &
               &                      DIM=3,         &
               &                      SHIFT=-1,      &
               &                      BOUNDARY=dl_value(:,:,1) )
            ENDIF
            IF( il_kmax == il_shape(3) )THEN
               dl_value(:,:,:)=EOSHIFT( dl_value(:,:,:), &
               &                        DIM=3,         &
               &                        SHIFT=1,       &
               &                        BOUNDARY=dl_value(:,:,3))
            ENDIF

            IF( ll_discont )THEN
               dl_min=MINVAL( dl_value(:,:,:), dl_value(:,:,:)/=dd_fill )
               dl_max=MAXVAL( dl_value(:,:,:), dl_value(:,:,:)/=dd_fill )
               IF( dl_min < -170_dp .AND. dl_max > 170_dp )THEN
                  WHERE( dl_value(:,:,:) < 0_dp ) 
                     dl_value(:,:,:) = dl_value(:,:,:)+360._dp
                  END WHERE
               ELSEIF( dl_min < 10_dp .AND. dl_max > 350_dp )THEN
                  WHERE( dl_value(:,:,:) > 180 ) 
                     dl_value(:,:,:) = dl_value(:,:,:)-180._dp
                  END WHERE
               ENDIF
            ENDIF         

            WHERE( dl_value(:,:, 2) /= dd_fill .AND. & ! jk
               &   dl_value(:,:, 3) /= dd_fill .AND. & ! jk+1
               &   dl_value(:,:, 1) /= dd_fill )       ! jk-1

               extrap_deriv_3D(:,:,jk)=&
               &  ( dl_value(:,:,3) - dl_value(:,:,1) ) / &
               &  REAL( il_kmax-il_kmin,dp)         

            END WHERE

         ENDDO

      END SELECT

      DEALLOCATE( dl_value )

   END FUNCTION extrap_deriv_3D
   !-------------------------------------------------------------------
   !> @brief
   !> This function compute coefficient for min_error extrapolation.
   !> 
   !> @details 
   !> coefficients are  "grid distance" to the center of the box choosed to compute
   !> extrapolation.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] dd_value  3D array of variable to be extrapolated
   !> @return 3D array of coefficient for minimum error extrapolation
   !-------------------------------------------------------------------
   PURE FUNCTION extrap__3D_min_error_coef( dd_value )

      IMPLICIT NONE
      ! Argument
      REAL(dp)   , DIMENSION(:,:,:), INTENT(IN) :: dd_value

      ! function
      REAL(dp), DIMENSION(SIZE(dd_value(:,:,:),DIM=1), &
      &                   SIZE(dd_value(:,:,:),DIM=2), &
      &                   SIZE(dd_value(:,:,:),DIM=3) ) :: extrap__3D_min_error_coef

      ! local variable
      INTEGER(i4), DIMENSION(3) :: il_shape

      INTEGER(i4) :: il_imid
      INTEGER(i4) :: il_jmid
      INTEGER(i4) :: il_kmid

      REAL(dp)   , DIMENSION(:,:,:), ALLOCATABLE :: dl_dist

      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: jj
      INTEGER(i4) :: jk
      !----------------------------------------------------------------

      ! init
      extrap__3D_min_error_coef(:,:,:)=0

      il_shape(:)=SHAPE(dd_value(:,:,:))

      il_imid=INT(REAL(il_shape(1),sp)*0.5+1)
      il_jmid=INT(REAL(il_shape(2),sp)*0.5+1)
      il_kmid=INT(REAL(il_shape(3),sp)*0.5+1)

      ALLOCATE( dl_dist(il_shape(1),il_shape(2),il_shape(3)) )

      DO jk=1,il_shape(3)
         DO jj=1,il_shape(2)
            DO ji=1,il_shape(1)

               ! compute distance
               dl_dist(ji,jj,jk) = (ji-il_imid)**2 + &
               &                   (jj-il_jmid)**2 + &
               &                   (jk-il_kmid)**2

               IF( dl_dist(ji,jj,jk) /= 0 )THEN
                  dl_dist(ji,jj,jk)=SQRT( dl_dist(ji,jj,jk) )
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      WHERE( dl_dist(:,:,:) /= 0 )
         extrap__3D_min_error_coef(:,:,:)=dl_dist(:,:,:)
      END WHERE

      DEALLOCATE( dl_dist )

   END FUNCTION extrap__3D_min_error_coef
   !-------------------------------------------------------------------
   !> @brief
   !> This function compute extrapolatd value by calculated minimum error using
   !> taylor expansion
   !> 
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !>
   !> @param[in] dd_value  3D array of variable to be extrapolated
   !> @param[in] dd_fill   FillValue of variable
   !> @param[in] id_radius radius of the halo used to compute extrapolation
   !> @param[in] dd_dfdx   derivative of function in i-direction
   !> @param[in] dd_dfdy   derivative of function in j-direction
   !> @param[in] dd_dfdz   derivative of function in k-direction
   !> @param[in] dd_coef   array of coefficient for min_error extrapolation 
   !> @return extrapolatd value
   !-------------------------------------------------------------------
   PURE FUNCTION extrap__3D_min_error_fill( dd_value, dd_fill, id_radius, &
   &                                        dd_dfdx, dd_dfdy, dd_dfdz, &
   &                                        dd_coef )
      IMPLICIT NONE
      ! Argument
      REAL(dp)   , DIMENSION(:,:,:), INTENT(IN) :: dd_value
      REAL(dp)   ,                   INTENT(IN) :: dd_fill
      INTEGER(i4),                   INTENT(IN) :: id_radius
      REAL(dp)   , DIMENSION(:,:,:), INTENT(IN) :: dd_dfdx
      REAL(dp)   , DIMENSION(:,:,:), INTENT(IN) :: dd_dfdy
      REAL(dp)   , DIMENSION(:,:,:), INTENT(IN) :: dd_dfdz
      REAL(dp)   , DIMENSION(:,:,:), INTENT(IN) :: dd_coef

      ! function
      REAL(dp) :: extrap__3d_min_error_fill

      ! local variable
      INTEGER(i4), DIMENSION(3) :: il_shape
      INTEGER(i4), DIMENSION(3) :: il_ind

      INTEGER(i4), DIMENSION(:,:,:), ALLOCATABLE :: il_mask
      REAL(dp)   , DIMENSION(:,:,:), ALLOCATABLE :: dl_error

      INTEGER(i4) :: il_min
      ! loop indices

      !----------------------------------------------------------------

      ! init
      extrap__3D_min_error_fill=dd_fill

      il_min=MAX(1,(SIZE(dd_value(:,:,:)))/(1+id_radius*2))

      IF( COUNT(dd_value(:,:,:) /= dd_fill) >= il_min )THEN

         il_shape(:)=SHAPE(dd_value(:,:,:))
         ALLOCATE( il_mask( il_shape(1),il_shape(2),il_shape(3)) )
         ALLOCATE( dl_error(il_shape(1),il_shape(2),il_shape(3)) )

         ! compute error
         dl_error(:,:,:)=0.
         il_mask(:,:,:)=0
         WHERE( dd_dfdx(:,:,:) /= dd_fill )
            dl_error(:,:,:)=dd_coef(:,:,:)*dd_dfdx(:,:,:)
            il_mask(:,:,:)=1
         END WHERE
         WHERE( dd_dfdy(:,:,:) /= dd_fill  )
            dl_error(:,:,:)=(dl_error(:,:,:)+dd_coef(:,:,:)*dd_dfdy(:,:,:))
            il_mask(:,:,:)=1
         END WHERE
         WHERE( dd_dfdz(:,:,:) /= dd_fill  )
            dl_error(:,:,:)=(dl_error(:,:,:)+dd_coef(:,:,:)*dd_dfdz(:,:,:))
            il_mask(:,:,:)=1
         END WHERE

         ! get minimum error indices
         il_ind(:)=MINLOC(dl_error(:,:,:),il_mask(:,:,:)==1)

         ! return value
         IF( ALL(il_ind(:)/=0) )THEN
            extrap__3D_min_error_fill=dd_value(il_ind(1),il_ind(2),il_ind(3))
         ENDIF

         DEALLOCATE( il_mask )
         DEALLOCATE( dl_error )

      ENDIF

   END FUNCTION extrap__3D_min_error_fill
   !-------------------------------------------------------------------
   !> @brief
   !> This function compute coefficient for inverse distance weighted method 
   !> 
   !> @details 
   !> coefficients are inverse "grid distance" to the center of the box choosed to compute
   !> extrapolation.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] dd_value  3D array of variable to be extrapolated
   !> @return 3D array of coefficient for inverse distance weighted extrapolation
   !-------------------------------------------------------------------
   PURE FUNCTION extrap__3D_dist_weight_coef( dd_value )

      IMPLICIT NONE
      ! Argument
      REAL(dp)   , DIMENSION(:,:,:), INTENT(IN) :: dd_value

      ! function
      REAL(dp), DIMENSION(SIZE(dd_value(:,:,:),DIM=1), &
      &                   SIZE(dd_value(:,:,:),DIM=2), &
      &                   SIZE(dd_value(:,:,:),DIM=3) ) :: extrap__3D_dist_weight_coef

      ! local variable
      INTEGER(i4), DIMENSION(3) :: il_shape

      INTEGER(i4) :: il_imid
      INTEGER(i4) :: il_jmid
      INTEGER(i4) :: il_kmid

      REAL(dp)   , DIMENSION(:,:,:), ALLOCATABLE :: dl_dist

      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: jj
      INTEGER(i4) :: jk
      !----------------------------------------------------------------

      ! init
      extrap__3D_dist_weight_coef(:,:,:)=0

      il_shape(:)=SHAPE(dd_value(:,:,:))

      il_imid=INT(REAL(il_shape(1),sp)*0.5+1,i4)
      il_jmid=INT(REAL(il_shape(2),sp)*0.5+1,i4)
      il_kmid=INT(REAL(il_shape(3),sp)*0.5+1,i4)

      ALLOCATE( dl_dist(il_shape(1),il_shape(2),il_shape(3)) )

      DO jk=1,il_shape(3)
         DO jj=1,il_shape(2)
            DO ji=1,il_shape(1)

               ! compute distance
               dl_dist(ji,jj,jk) = (ji-il_imid)**2 + &
               &                   (jj-il_jmid)**2 + &
               &                   (jk-il_kmid)**2

               IF( dl_dist(ji,jj,jk) /= 0 )THEN
                  dl_dist(ji,jj,jk)=SQRT( dl_dist(ji,jj,jk) )
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      WHERE( dl_dist(:,:,:) /= 0 ) 
         extrap__3D_dist_weight_coef(:,:,:)=1./dl_dist(:,:,:)
      END WHERE

      DEALLOCATE( dl_dist )

   END FUNCTION extrap__3D_dist_weight_coef
   !-------------------------------------------------------------------
   !> @brief
   !> This function compute extrapolatd value using inverse distance weighted
   !> method
   !> 
   !> @details 
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] dd_value  3D array of variable to be extrapolated
   !> @param[in] dd_fill   FillValue of variable
   !> @param[in] id_radius radius of the halo used to compute extrapolation
   !> @param[in] dd_coef   3D array of coefficient for inverse distance weighted extrapolation 
   !> @return extrapolatd value
   !-------------------------------------------------------------------
   FUNCTION extrap__3D_dist_weight_fill( dd_value, dd_fill, id_radius, &
   &                                     dd_coef )
      IMPLICIT NONE
      ! Argument
      REAL(dp)   , DIMENSION(:,:,:), INTENT(IN) :: dd_value
      REAL(dp)   ,                   INTENT(IN) :: dd_fill
      INTEGER(i4),                   INTENT(IN) :: id_radius
      REAL(dp)   , DIMENSION(:,:,:), INTENT(IN) :: dd_coef

      ! function
      REAL(dp) :: extrap__3D_dist_weight_fill

      ! local variable
      INTEGER(i4), DIMENSION(3) :: il_shape

      REAL(dp)   , DIMENSION(:,:,:), ALLOCATABLE :: dl_value
      REAL(dp)   , DIMENSION(:,:,:), ALLOCATABLE :: dl_coef

      INTEGER(i4) :: il_min
      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: jj
      INTEGER(i4) :: jk

      !----------------------------------------------------------------

      ! init
      extrap__3D_dist_weight_fill=dd_fill

      il_min=MAX(1,(SIZE(dd_value(:,:,:)))/(1+id_radius*2))

      IF( COUNT(dd_value(:,:,:)/= dd_fill) >= il_min )THEN

         il_shape(:)=SHAPE(dd_value(:,:,:))
         ALLOCATE( dl_value( il_shape(1),il_shape(2),il_shape(3)) )
         ALLOCATE( dl_coef( il_shape(1),il_shape(2),il_shape(3)) )

         dl_value(:,:,:)=0
         dl_coef(:,:,:)=0

         DO jk=1,il_shape(3)
            DO jj=1,il_shape(2)
               DO ji=1,il_shape(1)

                  IF( dd_value(ji,jj,jk) /= dd_fill )THEN
                     ! compute factor
                     dl_value(ji,jj,jk)=dd_coef(ji,jj,jk)*dd_value(ji,jj,jk)
                     dl_coef(ji,jj,jk)=dd_coef(ji,jj,jk)
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

         ! return value
         IF( SUM( dl_coef(:,:,:) ) /= 0 )THEN
            extrap__3D_dist_weight_fill = &
            &  SUM( dl_value(:,:,:) )/SUM( dl_coef(:,:,:) )
         ENDIF

         DEALLOCATE( dl_value )
         DEALLOCATE( dl_coef )

      ENDIF

   END FUNCTION extrap__3D_dist_weight_fill
   !-------------------------------------------------------------------
   !> @brief
   !> This subroutine add to the variable (to be extrapolated) an 
   !> extraband of N points at north,south,east and west boundaries.
   !> 
   !> @details
   !> optionaly you could specify size of extra bands in i- and j-direction
   !>
   !> @author J.Paul
   !> - November, 2013-Initial version
   !
   !> @param[inout] td_var variable 
   !> @param[in] id_isize  i-direction size of extra bands (default=im_minext)
   !> @param[in] id_jsize  j-direction size of extra bands (default=im_minext)
   !> @todo
   !> - invalid special case for grid with north fold
   !-------------------------------------------------------------------
   SUBROUTINE extrap_add_extrabands(td_var, id_isize, id_jsize )
      IMPLICIT NONE
      ! Argument
      TYPE(TVAR) , INTENT(INOUT)  :: td_var
      INTEGER(i4), INTENT(IN   ), OPTIONAL :: id_isize
      INTEGER(i4), INTENT(IN   ), OPTIONAL :: id_jsize

      ! local variable
      REAL(dp), DIMENSION(:,:,:,:) , ALLOCATABLE :: dl_value

      INTEGER(i4) :: il_isize
      INTEGER(i4) :: il_jsize
      INTEGER(i4) :: il_tmp

      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: ij
      !----------------------------------------------------------------
      il_isize=im_minext
      IF(PRESENT(id_isize)) il_isize=id_isize
      IF( il_isize < im_minext .AND. &
      &   TRIM(td_var%c_interp(1)) == 'cubic' )THEN
         CALL logger_warn("EXTRAP ADD EXTRABANDS: size of extrabands "//&
         &  "should be at least "//TRIM(fct_str(im_minext))//" for "//&
         &  " cubic interpolation ")
      ENDIF

      il_jsize=im_minext
      IF(PRESENT(id_jsize)) il_jsize=id_jsize
      IF( il_jsize < im_minext .AND. &
      &   TRIM(td_var%c_interp(1)) == 'cubic' )THEN
         CALL logger_warn("EXTRAP ADD EXTRABANDS: size of extrabands "//&
         &  "should be at least "//TRIM(fct_str(im_minext))//" for "//&
         &  " cubic interpolation ")
      ENDIF

      IF( .NOT. td_var%t_dim(1)%l_use ) il_isize=0
      IF( .NOT. td_var%t_dim(2)%l_use ) il_jsize=0

      CALL logger_trace( "EXTRAP ADD EXTRABANDS: dimension change "//&
      &              "in variable "//TRIM(td_var%c_name) )

      ! add extrabands in variable
      ALLOCATE(dl_value( td_var%t_dim(1)%i_len, &
      &                  td_var%t_dim(2)%i_len, &
      &                  td_var%t_dim(3)%i_len, &
      &                  td_var%t_dim(4)%i_len ))

      dl_value(:,:,:,:)=td_var%d_value(:,:,:,:)


      td_var%t_dim(1)%i_len = td_var%t_dim(1)%i_len + 2*il_isize
      td_var%t_dim(2)%i_len = td_var%t_dim(2)%i_len + 2*il_jsize

      DEALLOCATE(td_var%d_value)
      ALLOCATE( td_var%d_value(td_var%t_dim(1)%i_len, &
      &                        td_var%t_dim(2)%i_len, &
      &                        td_var%t_dim(3)%i_len, &
      &                        td_var%t_dim(4)%i_len ) )

      ! intialise
      td_var%d_value(:,:,:,:)=td_var%d_fill

      ! fill center
      td_var%d_value( 1+il_isize:td_var%t_dim(1)%i_len-il_isize, &
      &               1+il_jsize:td_var%t_dim(2)%i_len-il_jsize, &
      &                :,:) = dl_value(:,:,:,:)

      ! special case for overlap
      IF( td_var%i_ew >= 0 .AND. il_isize /= 0 )THEN
         DO ji=1,il_isize
            ! from east to west
            il_tmp=td_var%t_dim(1)%i_len-td_var%i_ew+ji-2*il_isize
            td_var%d_value(ji,:,:,:) = td_var%d_value(il_tmp,:,:,:)

            ! from west to east
            ij=td_var%t_dim(1)%i_len-ji+1
            il_tmp=td_var%i_ew-ji+2*il_isize+1
            td_var%d_value(ij,:,:,:) = td_var%d_value(il_tmp,:,:,:)
         ENDDO
      ENDIF

      DEALLOCATE( dl_value )

   END SUBROUTINE extrap_add_extrabands
   !-------------------------------------------------------------------
   !> @brief
   !> This subroutine remove of the variable an extraband 
   !> of N points at north,south,east and west boundaries.
   !> 
   !> @details
   !> optionaly you could specify size of extra bands in i- and j-direction
   !>
   !> @author J.Paul
   !> - November, 2013-Initial version
   !>
   !> @param[inout] td_var variable 
   !> @param[in] id_isize  i-direction size of extra bands (default=im_minext)
   !> @param[in] id_jsize  j-direction size of extra bands (default=im_minext)
   !-------------------------------------------------------------------
   SUBROUTINE extrap_del_extrabands(td_var, id_isize, id_jsize )
      IMPLICIT NONE
      ! Argument
      TYPE(TVAR) , INTENT(INOUT) :: td_var
      INTEGER(i4), INTENT(IN   ), OPTIONAL :: id_isize
      INTEGER(i4), INTENT(IN   ), OPTIONAL :: id_jsize

      ! local variable
      REAL(dp), DIMENSION(:,:,:,:) , ALLOCATABLE :: dl_value

      INTEGER(i4) :: il_isize
      INTEGER(i4) :: il_jsize
 
      INTEGER(i4) :: il_imin
      INTEGER(i4) :: il_imax
      INTEGER(i4) :: il_jmin
      INTEGER(i4) :: il_jmax

      ! loop indices
      !----------------------------------------------------------------
      il_isize=im_minext
      IF(PRESENT(id_isize)) il_isize=id_isize

      il_jsize=im_minext
      IF(PRESENT(id_jsize)) il_jsize=id_jsize

      IF( .NOT. td_var%t_dim(1)%l_use ) il_isize=0
      IF( .NOT. td_var%t_dim(2)%l_use ) il_jsize=0

      CALL logger_trace( "EXTRAP DEL EXTRABANDS: dimension change "//&
      &              "in variable "//TRIM(td_var%c_name) )

      ! add extrabands in variable
      ALLOCATE(dl_value( td_var%t_dim(1)%i_len, &
      &                  td_var%t_dim(2)%i_len, &
      &                  td_var%t_dim(3)%i_len, &
      &                  td_var%t_dim(4)%i_len ))

      dl_value(:,:,:,:)=td_var%d_value(:,:,:,:)

      ! fill center
      il_imin=1+il_isize
      il_imax=td_var%t_dim(1)%i_len-il_isize

      il_jmin=1+il_jsize
      il_jmax=td_var%t_dim(2)%i_len-il_jsize
      
      td_var%t_dim(1)%i_len = td_var%t_dim(1)%i_len - 2*il_isize
      td_var%t_dim(2)%i_len = td_var%t_dim(2)%i_len - 2*il_jsize

      DEALLOCATE(td_var%d_value)
      ALLOCATE( td_var%d_value(td_var%t_dim(1)%i_len, &
      &                        td_var%t_dim(2)%i_len, &
      &                        td_var%t_dim(3)%i_len, &
      &                        td_var%t_dim(4)%i_len ) )

      ! intialise
      td_var%d_value(:,:,:,:)=td_var%d_fill

      td_var%d_value(:,:,:,:)=dl_value(il_imin:il_imax,&
      &                                il_jmin:il_jmax,&
      &                                :,:)

      DEALLOCATE( dl_value )

   END SUBROUTINE extrap_del_extrabands
END MODULE extrap
