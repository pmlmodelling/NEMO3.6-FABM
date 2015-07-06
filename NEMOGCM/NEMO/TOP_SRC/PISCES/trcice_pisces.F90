MODULE trcice_pisces
   !!======================================================================
   !!                         ***  MODULE trcice_pisces  ***
   !! TOP :   initialisation of the PISCES biochemical model
   !!======================================================================
   !! History :  3.5  ! 2013    (M. Vancoppenolle, O. Aumont, G. Madec), original code
   !! Comment ! probably not properly done when the second particle export
   !! scheme (kriest) is used
   !!----------------------------------------------------------------------
#if defined key_pisces || defined key_pisces_reduced
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_ice_pisces   : PISCES fake sea ice model setting
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE par_pisces      ! PISCES parameters
   USE oce_trc         ! Shared variables between ocean and passive tracers
   USE trc             ! Passive tracers common variables 
   USE phycst          ! Ocean physics parameters
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE in_out_manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ice_ini_pisces ! called by trcini.F90 module

CONTAINS

   SUBROUTINE trc_ice_ini_pisces
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ice_ini_pisces ***
      !!
      !! ** Purpose :   PISCES fake sea ice model setting
      !!    Method  :   Assign prescribe values to tracer concentrations in sea ice
      !!
      !! For levitating sea ice, constant ocean tracer concentrations also have to be defined. 
      !! This is done specifically for Global, Arctic, Antarctic and Baltic regions
      !! 
      !! Sea ice concentrations are by default prescribed as follows
      !!  trc_i = zratio * trc_o
      !!
      !! This formulation is modulated by the namelist parameter trc_ice_ratio
      !!
      !! trc_ice_ratio  * betw 0 and 1: prescribed ice/ocean tracer concentration ratio
      !!                * -1 => the ice-ocean tracer concentration ratio follows the
      !!                         ice-ocean salinity ratio
      !!                * -2 => no ice-ocean tracer concentration is used
      !!                        instead, the tracer concentration in sea ice 
      !!                        is prescribed to trc_ice_prescr
      !! 
      !! cn_trc_o  specifies which disinctions are made for prescribed tracer concentration
      !!                * 'GL' use global ocean values making distinction for Baltic Sea only
      !!                * 'AA' use Arctic/Antarctic contrasted values, + Baltic
      !!
      !!----------------------------------------------------------------------

                                        !--- Dummy variables
      REAL(wp), DIMENSION(jptra,2) &
               ::  zratio            ! effective ice-ocean tracer cc ratio
      REAL(wp), DIMENSION(2) :: zrs  ! ice-ocean salinity ratio, 1 - global, 2- Baltic
      REAL(wp) :: zsice_bal          ! prescribed ice salinity in the Baltic
      REAL(wp) :: zsoce_bal          ! prescribed ocean salinity in the Baltic
      REAL(wp) :: zfeoce_glo         ! prescribed iron concentration in the global ocean
      REAL(wp) :: zfeoce_bal         ! prescribed iron concentration in the global ocean
      INTEGER  :: jn                 ! dummy loop index

      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ice_ini_pisces: Prescribed sea ice biogeochemistry '
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~'

      !--------------------------------------------
      ! Initialize ocean prescribed concentrations
      !--------------------------------------------
      ! values taken from a 500 yr equilibrium run
      ! used only in the levitating sea ice case with virtual salt / tracer
      ! fluxes

      !--- Global case 
      IF ( cn_trc_o(jpdic) == 'GL ' ) trc_o(:,:,jpdic) =  1.99e-3_wp 
      IF ( cn_trc_o(jpdoc) == 'GL ' ) trc_o(:,:,jpdoc) =  2.04e-5_wp 
      IF ( cn_trc_o(jptal) == 'GL ' ) trc_o(:,:,jptal) =  2.31e-3_wp 
      IF ( cn_trc_o(jpoxy) == 'GL ' ) trc_o(:,:,jpoxy) =  2.47e-4_wp
      IF ( cn_trc_o(jpcal) == 'GL ' ) trc_o(:,:,jpcal) =  1.04e-8_wp
      IF ( cn_trc_o(jppo4) == 'GL ' ) trc_o(:,:,jppo4) =  5.77e-7_wp / po4r 
      IF ( cn_trc_o(jppoc) == 'GL ' ) trc_o(:,:,jppoc) =  1.27e-6_wp  
#  if ! defined key_kriest
      IF ( cn_trc_o(jpgoc) == 'GL ' ) trc_o(:,:,jpgoc) =  5.23e-8_wp  
      IF ( cn_trc_o(jpbfe) == 'GL ' ) trc_o(:,:,jpbfe) =  9.84e-13_wp 
#  else
      IF ( cn_trc_o(jpnum) == 'GL ' ) trc_o(:,:,jpnum) = 0. ! could not get this value since did not use it
#  endif
      IF ( cn_trc_o(jpsil) == 'GL ' ) trc_o(:,:,jpsil) =  7.36e-6_wp  
      IF ( cn_trc_o(jpdsi) == 'GL ' ) trc_o(:,:,jpdsi) =  1.07e-7_wp 
      IF ( cn_trc_o(jpgsi) == 'GL ' ) trc_o(:,:,jpgsi) =  1.53e-8_wp
      IF ( cn_trc_o(jpphy) == 'GL ' ) trc_o(:,:,jpphy) =  9.57e-8_wp
      IF ( cn_trc_o(jpdia) == 'GL ' ) trc_o(:,:,jpdia) =  4.24e-7_wp
      IF ( cn_trc_o(jpzoo) == 'GL ' ) trc_o(:,:,jpzoo) =  6.07e-7_wp
      IF ( cn_trc_o(jpmes) == 'GL ' ) trc_o(:,:,jpmes) =  3.44e-7_wp
      IF ( cn_trc_o(jpfer) == 'GL ' ) trc_o(:,:,jpfer) =  4.06e-10_wp
      IF ( cn_trc_o(jpsfe) == 'GL ' ) trc_o(:,:,jpsfe) =  2.51e-11_wp
      IF ( cn_trc_o(jpdfe) == 'GL ' ) trc_o(:,:,jpdfe) =  6.57e-12_wp
      IF ( cn_trc_o(jpnfe) == 'GL ' ) trc_o(:,:,jpnfe) =  1.76e-11_wp
      IF ( cn_trc_o(jpnch) == 'GL ' ) trc_o(:,:,jpnch) =  1.67e-7_wp
      IF ( cn_trc_o(jpdch) == 'GL ' ) trc_o(:,:,jpdch) =  1.02e-7_wp
      IF ( cn_trc_o(jpno3) == 'GL ' ) trc_o(:,:,jpno3) =  5.79e-6_wp / rno3 
      IF ( cn_trc_o(jpnh4) == 'GL ' ) trc_o(:,:,jpnh4) =  3.22e-7_wp / rno3

      !--- Arctic specificities (dissolved inorganic & DOM)
      IF ( cn_trc_o(jpdic) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpdic) =  1.98e-3_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpdoc) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpdoc) =  6.00e-6_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jptal) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jptal) =  2.13e-3_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpoxy) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpoxy) =  3.65e-4_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpcal) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpcal) =  1.50e-9_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jppo4) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jppo4) =  4.09e-7_wp / po4r ; END WHERE ; ENDIF
      IF ( cn_trc_o(jppoc) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jppoc) =  4.05e-7_wp  ; END WHERE ; ENDIF
#  if ! defined key_kriest
      IF ( cn_trc_o(jpgoc) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpgoc) =  2.84e-8_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpbfe) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpbfe) =  7.03e-13_wp ; END WHERE ; ENDIF
#  else
      IF ( cn_trc_o(jpnum) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpnum) =  0.00e-00_wp ; END WHERE ; ENDIF
#  endif
      IF ( cn_trc_o(jpsil) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpsil) =  6.87e-6_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpdsi) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpdsi) =  1.73e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpgsi) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpgsi) =  7.93e-9_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpphy) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpphy) =  5.25e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpdia) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpdia) =  7.75e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpzoo) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpzoo) =  3.34e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpmes) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpmes) =  2.49e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpfer) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpfer) =  1.43e-9_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpsfe) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpsfe) =  2.21e-11_wp ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpdfe) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpdfe) =  2.04e-11_wp ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpnfe) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpnfe) =  1.75e-11_wp ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpnch) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpnch) =  1.46e-07_wp ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpdch) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpdch) =  2.36e-07_wp ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpno3) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpno3) =  3.51e-06_wp / rno3 ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpnh4) == 'AA ' ) THEN ; WHERE( gphit(:,:) >= 00._wp ) ; trc_o(:,:,jpnh4) =  6.15e-08_wp / rno3 ; END WHERE ; ENDIF

      !--- Antarctic specificities (dissolved inorganic & DOM)
      IF ( cn_trc_o(jpdic) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpdic) =  2.20e-3_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpdoc) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpdoc) =  7.02e-6_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jptal) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jptal) =  2.37e-3_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpoxy) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpoxy) =  3.42e-4_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpcal) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpcal) =  3.17e-9_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jppo4) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jppo4) =  1.88e-6_wp / po4r  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jppoc) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jppoc) =  1.13e-6_wp  ; END WHERE ; ENDIF
#  if ! defined key_kriest
      IF ( cn_trc_o(jpgoc) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpgoc) =  2.89e-8_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpbfe) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpbfe) =  5.63e-13_wp ; END WHERE ; ENDIF
#  else
      IF ( cn_trc_o(jpnum) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpnum) =  0.00e-00_wp ; END WHERE ; ENDIF
#  endif
      IF ( cn_trc_o(jpsil) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpsil) =  4.96e-5_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpdsi) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpdsi) =  5.63e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpgsi) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpgsi) =  5.35e-8_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpphy) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpphy) =  8.10e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpdia) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpdia) =  5.77e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpzoo) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpzoo) =  6.68e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpmes) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpmes) =  3.55e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpfer) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpfer) =  1.62e-10_wp ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpsfe) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpsfe) =  2.29e-11_wp ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpdfe) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpdfe) =  8.75e-12_wp ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpnfe) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpnfe) =  1.48e-11_wp ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpnch) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpnch) =  2.02e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpdch) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpdch) =  1.60e-7_wp  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpno3) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpno3) =  2.64e-5_wp / rno3  ; END WHERE ; ENDIF
      IF ( cn_trc_o(jpnh4) == 'AA ' ) THEN ; WHERE( gphit(:,:) <  00._wp ) ; trc_o(:,:,jpnh4) =  3.39e-7_wp / rno3  ; END WHERE ; ENDIF

      !--- Baltic Sea particular case for ORCA configurations
      IF( cp_cfg == "orca" ) THEN            ! Baltic mask
         WHERE( 14._wp <= glamt(:,:) .AND. glamt(:,:) <= 32._wp .AND.    &
                54._wp <= gphit(:,:) .AND. gphit(:,:) <= 66._wp )
         trc_o(:,:,jpdic) = 1.14e-3_wp
         trc_o(:,:,jpdoc) = 1.06e-5_wp
         trc_o(:,:,jptal) = 1.16e-3_wp
         trc_o(:,:,jpoxy) = 3.71e-4_wp
         trc_o(:,:,jpcal) = 1.51e-9_wp
         trc_o(:,:,jppo4) = 2.85e-9_wp / po4r
         trc_o(:,:,jppoc) = 4.84e-7_wp
#  if ! defined key_kriest
         trc_o(:,:,jpgoc) = 1.05e-8_wp
         trc_o(:,:,jpbfe) = 4.97e-13_wp
#  else
         trc_o(:,:,jpnum) = 0. ! could not get this value
#  endif
         trc_o(:,:,jpsil) = 4.91e-5_wp
         trc_o(:,:,jpdsi) = 3.25e-7_wp
         trc_o(:,:,jpgsi) = 1.93e-8_wp
         trc_o(:,:,jpphy) = 6.64e-7_wp
         trc_o(:,:,jpdia) = 3.41e-7_wp
         trc_o(:,:,jpzoo) = 3.83e-7_wp
         trc_o(:,:,jpmes) = 0.225e-6_wp
         trc_o(:,:,jpfer) = 2.45e-9_wp
         trc_o(:,:,jpsfe) = 3.89e-11_wp
         trc_o(:,:,jpdfe) = 1.33e-11_wp
         trc_o(:,:,jpnfe) = 2.62e-11_wp
         trc_o(:,:,jpnch) = 1.17e-7_wp
         trc_o(:,:,jpdch) = 9.69e-8_wp
         trc_o(:,:,jpno3) = 5.36e-5_wp / rno3
         trc_o(:,:,jpnh4) = 7.18e-7_wp / rno3
         END WHERE
      ENDIF ! cfg

      !-----------------------------
      ! Assign ice-ocean cc ratios 
      !-----------------------------
      ! 0 means zero concentration in sea ice
      ! 1 means same concentration in the sea ice as in the ocean

      ! Ice ocean salinity ratio
      zsoce_bal   = 4. ; zsice_bal   = 2. !! Baltic ocean and sea ice salinities
      zrs(1) = sice / soce                !! ice-ocean salinity ratio, global case
      zrs(2) = zsice_bal / zsoce_bal      !! ice-ocean salinity ratio, Baltic case

      DO jn = jp_pcs0, jp_pcs1
         IF ( trc_ice_ratio(jn) >= 0._wp )  zratio(jn,:) = trc_ice_ratio(jn)
         IF ( trc_ice_ratio(jn) == -1._wp ) zratio(jn,:) = zrs(:)
         IF ( trc_ice_ratio(jn) == -2._wp ) zratio(jn,:) = -9999.99_wp
      END DO

      !-------------------------------
      ! Sea ice tracer concentrations
      !-------------------------------
      DO jn = jp_pcs0, jp_pcs1
         !-- Everywhere but in the Baltic
         IF ( trc_ice_ratio(jn) >= -1._wp ) THEN !! no prescribed concentration
                                              !! (typically everything but iron) 
            trc_i(:,:,jn) = zratio(jn,1) * trc_o(:,:,jn) 
         ELSE                                 !! prescribed concentration
            trc_i(:,:,jn) = trc_ice_prescr(jn)
         ENDIF
       
         !-- Baltic
         IF( cp_cfg == "orca" ) THEN !! Baltic treated seperately for ORCA configs
            IF ( trc_ice_ratio(jn) >= - 1._wp ) THEN !! no prescribed concentration
                                                 !! (typically everything but iron) 
               WHERE( 14._wp <= glamt(:,:) .AND. glamt(:,:) <= 32._wp .AND.    &
                      54._wp <= gphit(:,:) .AND. gphit(:,:) <= 66._wp )
                     trc_i(:,:,jn) = zratio(jn,2) * trc_o(:,:,jn) 
               END WHERE
            ELSE                                 !! prescribed tracer concentration in ice
               WHERE( 14._wp <= glamt(:,:) .AND. glamt(:,:) <= 32._wp .AND.    &
                   54._wp <= gphit(:,:) .AND. gphit(:,:) <= 66._wp )
                     trc_i(:,:,jn) = trc_ice_prescr(jn)
               END WHERE
            ENDIF ! trc_ice_ratio
         ENDIF
      !
      END DO ! jn

   END SUBROUTINE trc_ice_ini_pisces

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                            No PISCES biochemical model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ice_ini_pisces         ! Empty routine
   END SUBROUTINE trc_ice_ini_pisces
#endif

   !!======================================================================
END MODULE trcice_pisces
