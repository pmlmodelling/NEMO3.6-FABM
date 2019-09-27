MODULE vertical_movement_fabm
   !!======================================================================
   !!                         ***  MODULE vertical_movement_fabm  ***
   !! TOP :   Module for the vertical movement of the FABM tracers
   !!======================================================================

#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                               FABM tracers
   !!----------------------------------------------------------------------
   !! compute_vertical_movement : compute vertical movement of FABM fields
   !!----------------------------------------------------------------------
   USE par_trc
   USE oce_trc
   USE trc
   USE par_fabm
   USE fabm
   USE dom_oce
#if defined key_trdtrc && defined key_iomput
   USE iom
   USE trdtrc_oce
#endif

   IMPLICIT NONE

#  include "domzgr_substitute.h90"

   PRIVATE

   PUBLIC compute_vertical_movement

   ! Work arrays for vertical advection (residual movement/sinking/floating)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:) :: w_ct
#if defined key_trdtrc && defined key_iomput
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:,:) :: tr_vmv
#endif

   CONTAINS

   SUBROUTINE compute_vertical_movement( kt, method )
      !!----------------------------------------------------------------------
      !!                     ***  compute_vertical_movement  ***
      !!
      !! ** Purpose : compute vertical movement of FABM tracers through the water
      !!              (sinking/floating/active movement)
      !!
      !! ** Method  : Retrieves additional vertical velocity field and applies
      !!              advection scheme.
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt     ! ocean time-step index
      INTEGER, INTENT(in) ::   method ! advection method (1: 1st order upstream, 3: 3rd order TVD with QUICKEST limiter)

      INTEGER :: ji,jj,jk,jn,k_floor
      REAL(wp) :: zwgt_if(0:jpk), flux(0:jpk), dc(1:jpk), w_if(0:jpk)
#if defined key_trdtrc
      CHARACTER (len=20) :: cltra
#endif

#if defined key_trdtrc && defined key_iomput
      IF( lk_trdtrc ) tr_vmv = 0.0_wp
#endif

      ! Compute interior vertical velocities and include them in source array.
      DO jj=1,jpj ! j-loop
         ! Get vertical velocities at layer centres (entire 1:jpi,1:jpk slice).
         DO jk=1,jpk
            CALL fabm_get_vertical_movement(model,1,jpi,jj,jk,w_ct(:,jk,:))
         END DO
         DO ji=1,jpi ! i-loop
            ! Only process this horizontal point (ji,jj) if number of layers exceeds 1
            k_floor = mbkt(ji,jj)
            IF (k_floor > 1) THEN
               ! Linearly interpolate to velocities at the interfaces between layers
               ! Note:
               !    - interface k sits between cell centre k and k+1 (k=0 for surface)
               !    - k [1,jpk] increases downwards
               !    - upward velocity is positive, downward velocity is negative
               zwgt_if(0) = 0._wp ! surface
               zwgt_if(1:k_floor-1) = fse3t(ji,jj,2:k_floor) / (fse3t(ji,jj,1:k_floor-1) + fse3t(ji,jj,2:k_floor))
               zwgt_if(k_floor:) = 0._wp ! sea floor and below
               ! Advect:
               DO jn=1,jp_fabm ! State loop
                  IF (ALL(w_ct(ji,1:k_floor,jn) == 0)) CYCLE

                  ! Compute velocities at interfaces
                  w_if(0) = 0._wp ! surface
                  w_if(1:k_floor-1) = zwgt_if(1:k_floor-1) * w_ct(ji,1:k_floor-1,jn) + (1._wp - zwgt_if(1:k_floor-1)) * w_ct(ji,2:k_floor,jn)
                  w_if(k_floor:) = 0._wp ! sea floor and below

                  IF (method == 1) THEN
                     CALL advect_1(k_floor, trn(ji,jj,1:k_floor,jp_fabm_m1+jn), w_if(0:k_floor), fse3t(ji,jj,1:k_floor), dc(1:k_floor))
                  ELSE
                     CALL advect_3(k_floor, trn(ji,jj,1:k_floor,jp_fabm_m1+jn), w_if(0:k_floor), fse3t(ji,jj,1:k_floor), dc(1:k_floor))
                  END IF

                  ! Compute change (per volume) due to vertical movement per layer
                  tra(ji,jj,1:k_floor,jp_fabm_m1+jn) = tra(ji,jj,1:k_floor,jp_fabm_m1+jn) + dc(1:k_floor)

#if defined key_trdtrc && defined key_iomput
                  ! Store change due to vertical movement as diagnostic
                  IF( lk_trdtrc .AND. ln_trdtrc( jp_fabm_m1+jn)) tr_vmv(ji,jj,1:k_floor,jn) = dc(1:k_floor)
#endif
               END DO ! State loop
            END IF ! Level check
         END DO ! i-loop
      END DO ! j-loop

#if defined key_trdtrc && defined key_iomput
      DO jn=1,jp_fabm ! State loop
        IF( lk_trdtrc .AND. ln_trdtrc(jp_fabm_m1+jn) ) THEN
          cltra = 'VMV_'//TRIM(ctrcnm(jp_fabm_m1+jn))
          CALL iom_put( cltra,  tr_vmv(:,:,:,jn) )
        END IF
      ENDDO
#endif

   END SUBROUTINE compute_vertical_movement

   SUBROUTINE advect_1(nk, c, w, h, trend)
      INTEGER,  INTENT(IN)  :: nk
      REAL(wp), INTENT(IN)  :: c(1:nk)
      REAL(wp), INTENT(IN)  :: w(0:nk)
      REAL(wp), INTENT(IN)  :: h(1:nk)
      REAL(wp), INTENT(OUT) :: trend(1:nk)

      REAL(wp) :: flux(0:nk)
      INTEGER  :: jk
      ! Compute fluxes (per surface area) over at interfaces (remember: positive for upwards)
      flux(0) = 0._wp
      DO jk=1,nk-1 ! k-loop
         IF (w_if(jk) > 0) THEN
            ! Upward movement (source layer is jk+1)
            flux(jk) = w(jk) * c(jk+1)
         ELSE
            ! Downward movement (source layer is jk)
            flux(jk) = w(jk) * c(jk)
         END IF
      END DO
      flux(nk) = 0._wp
      trend = (flux(1:nk) - flux(0:nk-1)) / h
   END SUBROUTINE

   SUBROUTINE advect_3(nk, c, w, h, trend)
      INTEGER,  INTENT(IN)  :: nk
      REAL(wp), INTENT(IN)  :: c(1:nk)
      REAL(wp), INTENT(IN)  :: w(0:nk)
      REAL(wp), INTENT(IN)  :: h(1:nk)
      REAL(wp), INTENT(OUT) :: trend(1:nk)

      INTEGER, PARAMETER :: n_itermax=100
      REAL(wp) :: z2dt,cmax_no
      REAL(wp) :: cfl(1:nk-1)
      INTEGER  :: n_iter, n_count, jk
      REAL(wp) :: c_new(1:nk)
      REAL(wp) :: tr_u(1:nk)
      REAL(wp) :: tr_c(1:nk)
      REAL(wp) :: tr_d(1:nk)
      REAL(wp) :: ratio(1:nk-1)
      REAL(wp) :: x_fac(1:nk-1)
      REAL(wp) :: phi_lim(1:nk-1)
      REAL(wp) :: limiter(1:nk-1)
      REAL(wp) :: flux_if(1:nk)

      IF( neuler == 0 .AND. kt == nittrc000 ) THEN
          z2dt = rdt                  ! set time step size (Euler)
      ELSE
          z2dt = 2._wp * rdt          ! set time step size (Leapfrog)
      ENDIF

      ! get maximum Courant number:
      cfl = abs(w(1:nk-1)) * z2dt / ( 0.5_wp*(h(2:nk) + h(1:nk-1)))
      cmax_no = MAXVAL(cfl)

      ! number of iterations:
      n_iter = min(n_itermax, int(cmax_no)+1)
      IF (ln_ctl.AND.(n_iter .gt. 1)) THEN
         WRITE(numout,*) 'compute_vertical_movement::advect_3():'
         WRITE(numout,*) '   Maximum Courant number is ',cmax_no,'.'
         WRITE(numout,*) '   ',n_iter,' iterations used for vertical advection.'
      ENDIF

      ! effective Courant number:
      cfl = cfl/n_iter

      DO n_count=1,n_iter ! Iterative loop
         !Compute slope ratio
         IF (nk.gt.2) THEN !More than 2 vertical wet points
            IF (nk.gt.3) THEN
               WHERE (w(2:nk-2).ge.0._wp) !upward movement
                  tr_u(2:nk-2)=c(4:nk)
                  tr_c(2:nk-2)=c(3:nk-1)
                  tr_d(2:nk-2)=c(2:nk-2)
               ELSEWHERE !downward movement
                  tr_u(2:nk-2)=c(1:nk-3)
                  tr_c(2:nk-2)=c(2:nk-2)
                  tr_d(2:nk-2)=c(3:nk-1)
               ENDWHERE
            ENDIF
            IF (w(1).ge.0._wp) THEN
               tr_u(1)=c(3)
               tr_c(1)=c(2)
               tr_d(1)=c(1)
            ELSE
               tr_u(1)=c(1)
               tr_c(1)=c(1)
               tr_d(1)=c(2)
            ENDIF
            IF (w(nk-1).ge.0._wp) THEN
               tr_u(nk-1)=c(nk)
               tr_c(nk-1)=c(nk)
               tr_d(nk-1)=c(nk-1)
            ELSE
               tr_u(nk-1)=c(nk-2)
               tr_c(nk-1)=c(nk-1)
               tr_d(nk-1)=c(nk)
            ENDIF
         ELSE !only 2 vertical wet points, i.e. only 1 interface
            IF (w(1).ge.0._wp) THEN
               tr_u(1)=c(2)
               tr_c(1)=c(2)
               tr_d(1)=c(1)
            ELSE
               tr_u(1)=c(1)
               tr_c(1)=c(1)
               tr_d(1)=c(2)
            ENDIF
         ENDIF
         WHERE (abs(tr_d - tr_c) .gt. 1.e-10_wp)
            ratio = (tr_c - tr_u) / (tr_d - tr_c)
         ELSEWHERE
            ratio = (tr_c - tr_u) * SIGN(1._wp, w(1:nk-1)) * 1.e10_wp
         ENDWHERE

         !QUICKEST flux limiter:
         x_fac = (1._wp - 2._wp * cfl) / 6._wp
         phi_lim = (0.5_wp + x_fac) + (0.5_wp - x_fac) * ratio
         limiter = max(0._wp, min( phi_lim, 2._wp / (1._wp - cfl), 2._wp * ratio / (cfl + 1.e-10_wp)))

         ! Compute limited flux:
         flux_if = w(1:nk-1) * (tr_c + 0.5_wp * limiter * (1._wp - cfl) * (tr_d - tr_c))

         ! Compute pseudo update for trend aggregation:
         c_new = c
         c_new(1:nk-1) = c_new(1:nk-1) + z2dt / real(n_iter, wp) / h(1:nk-1) * flux_if(2:nk)
         c_new(2:nk)   = c_new(2:nk)   - z2dt / real(n_iter, wp) / h(2:nk)   * flux_if(2:nk)

      ENDDO ! Iterative loop

      ! Estimate rate of change from pseudo state updates (source
      ! splitting):
      trend = (c_new - c) / z2dt
   END SUBROUTINE

#endif
END MODULE
