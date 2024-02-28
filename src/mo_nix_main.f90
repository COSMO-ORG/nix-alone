!> Swiss snow cover scheme NIX (i.e. latin for snow)
!!
!! --------------------------------------------------------------------
!! -------------------------------------------------------------------
!!
!! @par Description:
!!  This module 'mo_snow_main.f90' performs all snow cover related
!!  calculations. It contains calls to all required physics of the
!!  Swiss  muli-layer snow cover scheme NIX (formerly known as SNOWPOLINO) developed
!!  in a joint project between WSL Institute for Snow  and Avaalnche Research SLF
!!  and MeteoSwiss. Model physics are outlined in seperate module and subroutines
!!  in /src/snow
!!
!!  All paramteric scalar and array data for this snow model routine are
!!  defined in tha data module 'mo_snow_data.f90'. All global fields are
!!  passed to the snow model by argument list
!!
!!  All global scalar variables of the module that are used by the snow
!!  model routine are imported by USE statements.
!!
!!  Potential initiations are done within the module 'mo_snow_init.f90'
!!
!! @author: S. Bellaire (MeteoSwiss)
!!
!! @par Reference ADD PUBLICATIONS ONCE AVAILABLE
!!
!! @par Revision History
!!  implemented into ICON by S. Bellaire (2022-10)
!!
!! @par Copyright and License
!!  This code is subject to the DWD and MPI-M-Software-License-Agreement in
!!  its most recent form.
!!  Please see the file LICENSE in the root of the source tree for this code.
!!  Where software is supplied by third parties, it is indicated in the
!!  headers of the routines.
!!
!! ----------------------------------------------------------------------
!! ----------------------------------------------------------------------

! ------------------------------------------------------------------------------
! Begin of module nix_main
! ------------------------------------------------------------------------------

MODULE mo_nix_main

   USE mo_kind,                ONLY: wp

   USE sfc_terra_data,         ONLY: eps_div

   USE mo_physical_constants,  ONLY: rho_w => rhoh2o  ! density if liquid water

   USE mo_nix_config,          ONLY: min_height_layer, max_height_layer

   USE mo_nix_meteo_util,      ONLY: calc_precip                   , &
   &                               calc_wind                     , &
   &                               calculate_tch                 , &
   &                               calculate_turbulent_fluxes    , &
   &                               calculate_radiative_fluxes    , &
   &                               calculate_atmospheric_forcing

   USE mo_nix_snow_util,       ONLY: calc_hn_density , &
   &                               calc_hn_amounts , &
   &                               update_nix_state

   USE mo_nix_stratigraphy,    ONLY: nix_stratigraphy, &
      split_bottom_layer

   USE mo_nix_heat_equation,   ONLY: heat_equation_wrapper

   USE mo_nix_substrate,       ONLY: calculate_soil_properties

   USE mo_nix_phase_change,    ONLY: phase_change

   USE mo_nix_misc_util,       ONLY: save_nix_state

   USE mo_nix_water_transport, ONLY: surface_mass_flux, &
      water_transport

   USE mo_nix_settling,        ONLY: aggregate_layers , &
      settling

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------------------

   PUBLIC :: nix_core

CONTAINS

! =============================================================================
! + Begin subroutine: NIX core routine
! ============================================================================

   SUBROUTINE nix_core(             &
   ! Utility variables (IN)
   &     nvec                             , & ! IN array dimensions
   &     ivstart                          , & ! IN optional start/end indicies
   &     ivend                            , & ! IN optional start/end indicies
   &     iblock                           , & ! number of block
   &     ke_soil                          , & ! IN without lowermost (climat.) soil layer
   &     ke_snow                          , & ! IN maximum number of snow layers
   &     nclass_gscp                      , & ! IN number of hydrometeor classes
   &     dt                               , & ! IN time step
   ! NIX input variables (IN)
   &     prr_con                          , & ! IN precipitation rate of rain, convective       (kg/m2*s)
   &     prs_con                          , & ! IN precipitation rate of snow, convective       (kg/m2*s)
   &     prr_gsp                          , & ! IN precipitation rate of rain, grid-scale       (kg/m2*s)
   &     prs_gsp                          , & ! IN precipitation rate of snow, grid-scale       (kg/m2*s)
   &     prg_gsp                          , & ! IN precipitation rate of graupel, grid-scale    (kg/m2*s)
   &     u                                , & ! IN zonal wind speed
   &     v                                , & ! IN meridional wind speed
   &     t                                , & ! IN temperature                                  (  K  )
   &     qv                               , & ! IN specific water vapor content                 (kg/kg)
   &     ps                               , & ! IN surface pressure                             ( Pa  )

   ! NIX prognostic variables (INOUT)
   &     t_sn                           , & ! INOUT snow layer temperature (main level)
   &     t_sn_n                           , & ! INOUT snow temperature at NODES
   &     theta_i                          , & ! INOUT volumetric ice content
   &     theta_w                          , & ! INOUT volumetric water content
   &     theta_a                          , & ! INOUT volumetric air content
   &     dzm_sn                           , & ! INOUT snow layer thickness
   &     hn_sn                            , & ! INOUT new snow amount
   &     top_sn                           , & ! INOUT index of first (top) snow layer
   ! NIX specific fluxes and related variables
   &     t_sn_sfc                         , & ! INOUT snow surface temperature                    [ K ]
!
   &     swflx_sn_net                     , & ! INOUT short wave radiation - net                  [W/m**2]
   &     swflx_sn_dn                      , & ! INOUT short wave radiation - downward             [W/m**2]
   &     swflx_sn_up                      , & ! INOUT short wave radiation - upward               [W/m**2]
   &     swflx_sn_abs                     , & ! INOUT short wave radiation - absorbed             [W/m**2]
!
   &     lwflx_sn_up                      , & ! INOUT long wave radiation  - upward               [W/m**2]
   &     lwflx_sn_dn                      , & ! INOUT long wave radiation  - downward             [W/m**2]
!
   &     alpha_sn                         , & ! INOUT surface albedo (snow)                       [-]
!
   &     shflx_sn                         , & ! INOUT turbulent sensible heat flux                [W/m**2]
   &     lhflx_sn                         , & ! INOUT turbulent latent heat flux                  [W/m**2]
   &     tch_sn                           , & ! INOUT transfer coefficient                        [    ]
!
   &     for_sn                           , & ! INOUT total atmospheric forcing at snow surface   [W/m**2]
   &     qv0_sn                           , & ! INOUT specific humidity at snow surface
   ! Heat equation
   &     hcap_sn                          , & ! INOUT snow layer heat capacity                    [   ]
   &     hcon_sn                          , & ! INOUT snow layer heat conductivity                [   ]
   &     hdif_sn                          , & ! INOUT snow layer heat diffusivity                 [   ]
   ! Additional fields
   &     rho_sn                           , & ! INOUT snow layer density                          [kg/m**3]
   &     mass_sn                          , & ! INOUT snow layer mass                             [kg]
   &     hm_sn                            , & ! INOUT height (from bottom) of main snow level     [ m ]
   &     zm_sn                            , & ! INOUT depth (from top) of main snow level         [ m ]
   &     swe_sn                           , & ! INOUT snow water equivialent                      [   ]
   &     runoff_sn                        , &
   &     h_snow                           , & ! INOUT total melt runoff                           [ m ]
   ! Soil properties
   &     t_so                             )   ! IN processing soil level structure


      ! Subroutine Arguments
      INTEGER,                           INTENT(IN)   :: &
         nvec,        & ! < array dimensions
         ivstart,     & ! < start index for computations in the parallel program
         ivend,       & ! < end index for computations in the parallel program
         iblock,      & ! number of block
         ke_soil,     & ! < number of soil layers (minus climatological layer)
         ke_snow,     & ! < maximum number of snow layers (nlev_sn)
         nclass_gscp    ! < number of hydrometeor classes of grid scale microphysics


      REAL    (KIND = wp), INTENT(IN)  ::  &
         dt                   ! time step



      REAL (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
         prr_con    , & ! precipitation rate of rain, convective        (kg/m2*s)
         prs_con    , & ! precipitation rate of snow, convective        (kg/m2*s)
         prr_gsp    , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
         prs_gsp    , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
         prg_gsp    , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
         u          , & ! zonal wind speed                              (m/s)
         v          , & ! meridional wind speed                         (m/s)
         t          , & ! temperature                                   (K)
         qv         , & ! specific water vapor content                  (kg/kg)
         ps             ! surface pressure                              ( Pa  )

      REAL (KIND = wp), DIMENSION(nvec, ke_snow),    INTENT(INOUT) :: &
         t_sn                           , & ! snow layer temperature (main level)
         theta_i                          , & ! volumetric ice content
         theta_w                          , & ! volumetric water content
         theta_a                          , & ! volumetric air content
         dzm_sn                               ! snow layer thickness

      REAL (KIND=wp), DIMENSION(nvec, ke_snow+1), INTENT(INOUT) :: t_sn_n ! snow temperature at nodes

      INTEGER, DIMENSION(nvec), INTENT(INOUT) :: &
         top_sn


      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
         hn_sn                            , & ! new snow amount
         t_sn_sfc                         , & ! snow surface temperature                    [ K ]
         swflx_sn_net                     , & ! short wave radiation - net                  [W/m**2]
         swflx_sn_dn                      , & ! short wave radiation - downward             [W/m**2]
         swflx_sn_up                          ! short wave radiation - upward               [W/m**2]

      REAL (KIND = wp), DIMENSION(nvec, ke_snow),    INTENT(INOUT) :: &
         swflx_sn_abs                         ! INOUT short wave radiation - absorbed             [W/m**2]

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
         lwflx_sn_up                       , & ! long wave radiation  - upward               [W/m**2]
         lwflx_sn_dn                       , & ! long wave radiation  - downward             [W/m**2]
         alpha_sn                         , & ! surface albedo (snow)                       [-]
         shflx_sn                         , & ! turbulent sensible heat flux                [W/m**2]
         lhflx_sn                         , & ! turbulent latent heat flux                  [W/m**2]
         tch_sn                           , & ! transfer coefficient                        [    ]
         for_sn                           , & ! total atmospheric forcing at snow surface   [W/m**2]
         qv0_sn                               ! specific humidity ar snow surface

      REAL (KIND = wp), DIMENSION(nvec, ke_snow),    INTENT(INOUT) :: &
         hcap_sn                          , & ! snow layer heat capacity                    [   ]
         hcon_sn                          , & ! snow layer heat conductivity                [   ]
         hdif_sn                          , & ! snow layer heat diffusivity                 [   ]
         rho_sn                           , & ! snow layer density                          [kg/m**3]
         mass_sn                          , & ! snow layer mass                             [kg]
         hm_sn                            , & ! height (from bottom) of main snow level     [ m ]
         zm_sn                                ! depth (from top) of main snow level         [ m ]

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
!
         swe_sn                           , & ! snow water equivialent                      [   ]
         runoff_sn                            ! total melt runoff                           [   ]

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
         h_snow                               ! snow height                                 [ m ]

      REAL (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
         t_so                                 ! soil layer temperature

      ! Local variables
      INTEGER ::           &
         i                                , & ! loop index in x-direction
         ksn                                  ! loop index in z-drection (snow layers)

      REAL (KIND = wp), DIMENSION(nvec) :: &
         zsnow_rate                       , & ! snow fall rate                              [mm/s]
         zrain_rate                       , & ! rain fall rate                              [mm/s]
         zuv                              , & ! wind speed                                  [m/s]
         rho_hn                               ! new snow density                            [kg/m3]

      INTEGER, DIMENSION(nvec)          :: &
         top              ! index of first (top) snow layer  [-]

      REAL (KIND = wp), DIMENSION(nvec) :: hcon_so

      REAL (KIND=wp), DIMENSION(nvec,ke_snow) ::  &
         zm_sn_old                        , & ! old value of layer depth
         theta_i_old                      , & ! volumetric ice content
         theta_w_old                          ! volumetric water content


      INTEGER :: my_cart_id

      INTEGER :: my_thrd_id, mcid, mtid, mbid, mvid

      LOGICAL :: ldebug = .TRUE.


      CHARACTER(len=12)     :: str_module    = 'nix_main'  ! Output of module for 1 line debug

      ! imposing a soil conductivity for standalone applications.
      hcon_so(:) = 0.8_wp

      ! ----------------------
      ! Section 0.2 -  Initializations
      ! ----------------------

      !Initiate index of top snow layer - FIXME: If we use top_sn directly we can spare these lines
      !                                          top_sn(i) is not a REAL anymore
      DO i = ivstart, ivend
         top(i)    = top_sn(i)
      ENDDO


! ------------------------------------------------------------------------------
! Section 1: Stratigraphy
! ------------------------------------------------------------------------------


      ! ----------------------
      ! Section 1.1 -  Preparations
      ! ----------------------

      ! Calculate precipitations rates (snow and rain) - FIXME: Occasionally idividual precipitation rates can get zero. Ignore for now!
      CALL calc_precip(nvec, ivstart, ivend,       &
         nclass_gscp,                &
         zsnow_rate, zrain_rate,     &
         prr_con, prs_con, prr_gsp,  &
         prs_gsp, prg_gsp )

      ! Calculate wind speed
      CALL calc_wind(nvec, ivstart, ivend, zuv, u, v)

      ! Calculate new snow density
      CALL calc_hn_density(nvec, ivstart, ivend, rho_hn, zuv, t)

      ! Calculate new snow amounts
      CALL calc_hn_amounts(nvec, ivstart, ivend, hn_sn, zsnow_rate, dt)

      ! ----------------------
      ! Section 1.2 - Stratigraphy - Main routine
      ! ----------------------

      ! Main routine for stratigraphy
      CALL nix_stratigraphy(nvec, ivstart, ivend           , &
      &             top, ke_snow, dzm_sn, rho_sn   , &
      &             theta_i, theta_w, theta_a      , &
      &             t_sn, t_sn_n, hn_sn, rho_hn, t, alpha_sn )


      ! Update layering
      CALL update_nix_state(nvec, ivstart, ivend, top, ke_snow , &
      &                   dzm_sn, rho_sn                     , &
      &                   theta_i, theta_w, theta_a          , &
      &                   t_sn, t_sn_n, hm_sn, zm_sn                 , &
      &                   hcap_sn, hcon_sn,  mass_sn         , &
      &                   h_snow, t_sn_sfc )

! ------------------------------------------------------------------------------
! Section 2: Heat Equation
! ------------------------------------------------------------------------------

      ! ----------------------
      ! Section 2.1 -  Preparations
      ! ----------------------

      ! Calculate transfer coefficient
      CALL calculate_tch(nvec, ivstart, ivend, t_sn_sfc, t, tch_sn, zuv, qv, ps)


      ! Calcualte turbulent fluxes
      CALL calculate_turbulent_fluxes(nvec, ivstart, ivend, lhflx_sn, shflx_sn  , &
      &                             t, t_sn_sfc, qv, qv0_sn                   , &
      &                             tch_sn, zuv, ps)


      ! Calculate radiative fluxes
      CALL calculate_radiative_fluxes(nvec, ivstart, ivend, ke_snow, top       , &
      &                             swflx_sn_net, swflx_sn_dn, swflx_sn_up   , &
      &                             swflx_sn_abs, lwflx_sn_up, lwflx_sn_dn   , &
      &                             alpha_sn, t_sn_sfc                       , &
      &                             rho_sn, dzm_sn, dt )


      ! ----------------------
      ! Section 2.2 - Heat equation main routine
      ! ----------------------

      CALL heat_equation_wrapper(nvec, ivstart, ivend, &
      & ke_snow, top,  &
      & dzm_sn, hcon_sn, hcap_sn, hdif_sn   , &
      & t_sn, t_sn_n, &
      & swflx_sn_abs,lwflx_sn_dn, &
      & lwflx_sn_up,lhflx_sn, &
      & shflx_sn,  &
      & hcon_so, &
      & t_so,dt,t_sn_sfc,tch_sn,rho_sn, t, theta_w)


! ------------------------------------------------------------------------------
! Section 3: Phase Change
! ------------------------------------------------------------------------------

      ! Save nix state for later use to local variables
      CALL save_nix_state(nvec, ivstart, ivend, ke_soil, ke_snow   , &
      &                 zm_sn_old, theta_i_old, theta_w_old      , &
      &                 zm_sn    , theta_i    , theta_w )

      ! Call main phase change routine
      CALL phase_change(nvec, ivstart, ivend, ke_soil, ke_snow  , &
      &               top, dzm_sn, theta_i, theta_w, theta_a  , &
      &               t_sn, t_sn_n, hcap_sn, rho_sn, dt)


      ! Update layering
      CALL update_nix_state(nvec, ivstart, ivend, top, ke_snow , &
      &                   dzm_sn, rho_sn                     , &
      &                   theta_i, theta_w, theta_a          , &
      &                   t_sn, t_sn_n, hm_sn, zm_sn                 , &
      &                   hcap_sn, hcon_sn,  mass_sn         , &
      &                   h_snow, t_sn_sfc )


! ------------------------------------------------------------------------------
! Section 4: Water Transport
! ------------------------------------------------------------------------------

      ! Calculate sublimation/deposition/evaporation
      CALL surface_mass_flux(nvec, ivstart, ivend, ke_soil, ke_snow    , &
      &                 top, lhflx_sn, dt, dzm_sn, rho_sn   , &
      &                 theta_i, theta_w, t_sn_sfc)

      ! Update layering
      CALL update_nix_state(nvec, ivstart, ivend, top, ke_snow , &
      &                   dzm_sn, rho_sn                     , &
      &                   theta_i, theta_w, theta_a          , &
      &                   t_sn, t_sn_n,hm_sn, zm_sn                 , &
      &                   hcap_sn, hcon_sn,  mass_sn         , &
      &                   h_snow, t_sn_sfc )


      ! Call main water transport routine
      CALL water_transport(nvec, ivstart, ivend, ke_soil, ke_snow  , &
      &               top, lhflx_sn, dt, dzm_sn, rho_sn       , &
      &               theta_i, theta_w, theta_a, hcap_sn      , &
      &               runoff_sn, t_sn)



! ------------------------------------------------------------------------------
! Section 5: Settling
! ------------------------------------------------------------------------------

      ! Aggregate layers
      CALL aggregate_layers(nvec, ivstart, ivend, ke_soil, ke_snow     , &
      &                   top, dzm_sn, hm_sn, zm_sn                  , &
      &                   theta_i, theta_i_old, theta_w, theta_a     , &
      &                   hcap_sn, hcon_sn, hdif_sn                  , &
      &                   t_sn,t_sn_n,rho_sn, mass_sn, dt)

      ! Update nix state
      CALL update_nix_state(nvec, ivstart, ivend, top, ke_snow , &
      &                   dzm_sn, rho_sn                     , &
      &                   theta_i, theta_w, theta_a          , &
      &                   t_sn, t_sn_n,hm_sn, zm_sn                 , &
      &                   hcap_sn, hcon_sn,  mass_sn         , &
      &                   h_snow, t_sn_sfc )

      ! Calculate settling
      CALL settling(nvec, ivstart, ivend, ke_soil, ke_snow     , &
      &           top, dzm_sn, theta_i, theta_i_old, theta_w  , &
      &           t_sn, rho_sn, mass_sn, dt, zm_sn, h_snow)

! ------------------------------------------------------------------------------
! Section 6: Post-processing
! ------------------------------------------------------------------------------

      ! ----------------------
      ! Special cases
      ! ----------------------

      ! Split bottom layer if feasible and top(i) == 1
      CALL  split_bottom_layer(nvec, ivstart, ivend         , &
      &              top, ke_snow, dzm_sn         , &
      &              theta_i, theta_w, theta_a    , &
      &              t_sn, t_sn_n )

      ! Update nix state
      CALL update_nix_state(nvec, ivstart, ivend, top, ke_snow , &
      &                   dzm_sn, rho_sn                     , &
      &                   theta_i, theta_w, theta_a          , &
      &                   t_sn, t_sn_n,hm_sn, zm_sn                 , &
      &                   hcap_sn, hcon_sn,  mass_sn         , &
      &                   h_snow, t_sn_sfc )

      ! ----------------------
      ! Update prognostic fields with local fields
      ! ----------------------

      ! Snow height FIXME: Maybe do this in update_nix_state?
      DO i = ivstart, ivend

         IF (top(i) .GE. 1 ) THEN
            h_snow(i) = hm_sn(i,top(i))
         ELSE
            h_snow(i) = 0.0_wp
         ENDIF

      ENDDO

      ! Copy top to top_sn and convert FIXME: See comment above
      DO i = ivstart, ivend
         top_sn(i) = top(i)        ! index of first (top) snow layer
      END DO

      ! ----------------------
      ! Calculate snow water equivalent
      ! ----------------------

      DO i = ivstart, ivend ! FIXME probably also something for update_nix_state()

         swe_sn(i) = 0.0_wp  ! Initiate to zero

         DO ksn=1,ke_snow
            swe_sn(i) = swe_sn(i) + ( dzm_sn(i,ksn) * rho_sn(i,ksn) ) / rho_w
         ENDDO

      END DO


      ! ----------------------
      ! Output fields to standard output - FIXME: Add message level here and reduce output!
      ! ----------------------

      ! Do some plausibility printing/writing
      DO i = ivstart, ivend

        WRITE(*,*) top(i), h_snow(i), t_sn_sfc(i), swflx_sn_dn(i), swflx_sn_up(i)  , &
          &        lwflx_sn_up(i), lwflx_sn_dn(i)                                  , &
          &        shflx_sn(i), lhflx_sn(i)                                        , &
          &        dzm_sn(i,:), theta_i(i,:), theta_w(i,:)                        , &
          &        t_sn(i,:)

      ENDDO


! =============================================================================
! - END subroutine: NIX
! ============================================================================

   END SUBROUTINE nix_core




!------------------------------------------------------------------------------
! End of module sfc_snow
!------------------------------------------------------------------------------

END MODULE mo_nix_main




