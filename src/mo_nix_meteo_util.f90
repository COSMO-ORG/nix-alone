!> Swiss snow cover scheme NIX (i.e. latin for snow)
!!
!! --------------------------------------------------------------------
!! -------------------------------------------------------------------
!!
!! @par Description:
!!  This module ...
!!
!! @author:
!!
!! @par Reference ADD PUBLICATIONS ONCE AVAILABLE
!!
!! @par Revision History
!!
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
! Begin of module mo_nix_meteo_util
! ------------------------------------------------------------------------------

MODULE mo_nix_meteo_util

   USE mo_kind,                    ONLY: wp

   USE mo_nix_config,              ONLY: z0_sn

   USE mo_physical_constants,      ONLY: r_d     => rd      , & ! gas constant for dry air
      rvd_m_o => vtmpc1  , & ! r_v/r_d - 1
      lh_v    => alv     , & ! latent heat of vapourization
      cp_d    => cpd     , & ! specific heat of dry air at constant press
      stbo    => stbo    , & ! Stefan Boltzman Konstante
      t0_melt => tmelt       ! melting temperature of ice/snow

   USE mo_nix_constants,             ONLY: ctalb, eps_div

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: calc_precip
   PUBLIC :: calc_wind
   PUBLIC :: calculate_tch
   PUBLIC :: calculate_turbulent_fluxes
   PUBLIC :: calculate_radiative_fluxes
   PUBLIC :: calculate_atmospheric_forcing

CONTAINS

! =============================================================================
! + Begin subroutine: calc_precip
! ============================================================================

   SUBROUTINE calc_precip(nvec, ivstart, ivend       , &
      nclass_gscp                , &
      zsnow_rate, zrain_rate     , &
      prr_con, prs_con, prr_gsp  , &
      prs_gsp, prg_gsp )

      ! Subroutine arguments
      INTEGER                          :: &
      !
         nvec              , & ! < array dimensions
         ivstart           , & ! < start index for computations in the parallel program
         ivend             , & ! < end index for computations in the parallel program
         nclass_gscp           ! < number of hydrometeor classes of grid scale microphysics

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
      !
         zsnow_rate        , & ! < rate of snow fall
         zrain_rate           ! < rate of rain fall

      REAL (KIND = wp), DIMENSION(nvec) :: &
      !
         prr_con          , & ! precipitation rate of rain, convective        (kg/m2*s)
         prs_con          , & ! precipitation rate of snow, convective        (kg/m2*s)
         prr_gsp          , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
         prs_gsp          , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
         prg_gsp              ! precipitation rate of graupel, grid-scale     (kg/m2*s)


      ! Local variables
      INTEGER :: &
         i                    ! loop index in x-direction


      ! ------------------------------------------------------------------------------
      ! Calculate precipitations rate
      ! ------------------------------------------------------------------------------

      DO i = ivstart, ivend

         ! ------------------------------
         ! Snow rate
         ! ------------------------------

         IF ( nclass_gscp >= 6 ) THEN
            zsnow_rate(i) = prs_gsp(i)+prs_con(i)+prg_gsp(i)            ! [kg/m**2 s]
         ELSE
            zsnow_rate(i) = prs_gsp(i)+prs_con(i)                       ! [kg/m**2 s]
         ENDIF

         ! ------------------------------
         ! Rain rate
         ! ------------------------------

         zrain_rate(i) = prr_gsp(i)+prr_con(i)  ! [kg/m**2 s]

      END DO



! =============================================================================
! - END subroutine: calc_precip
! ============================================================================

   END SUBROUTINE calc_precip

! =============================================================================
! + Begin subroutine: calc_wind
! ============================================================================

   SUBROUTINE calc_wind(nvec, ivstart, ivend, zuv, u, v)

      ! Subroutine arguments
      INTEGER                          :: &
      !
         nvec              , & ! < array dimensions
         ivstart           , & ! < start index for computations in the parallel program
         ivend                 ! < end index for computations in the parallel program

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
      !
         zuv                   ! < wind speed

      REAL (KIND = wp), DIMENSION(nvec) :: &
      !
         u                 , & ! < zonal component of wind
         v                     ! < meridional component of wind


      ! Local variables
      INTEGER :: &
         i                    ! loop index in x-direction


      ! ------------------------------------------------------------------------------
      ! Calculate wind speed
      ! ------------------------------------------------------------------------------

      DO i = ivstart, ivend

         zuv(i)        = SQRT ( u(i)**2 + v(i)**2 )

      END DO



! =============================================================================
! - END subroutine: calc_wind
! ============================================================================

   END SUBROUTINE calc_wind


! =============================================================================
! + Begin subroutine: calculate_tch
! ============================================================================

   SUBROUTINE calculate_tch(nvec, ivstart, ivend, t_sn_sfc, t, tch_sn, zuv, qv, ps)

      ! Subroutine arguments
      INTEGER                          :: &
      !
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend                  ! < end index for computations in the parallel program

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
      !
         tch_sn                 ! < wind speed

      REAL (KIND = wp), DIMENSION(nvec) :: &
      !
         t_sn_sfc           , & ! < snow surface temperature
         zuv                , & ! < wind speed
         t                  , & ! < air temperature first atmospheric level
         qv                 , & ! < specific water vapour content
         ps

      ! Local variables
      INTEGER :: &
         i                      ! loop index in x-direction

      REAL (KIND = wp) ::  &

         z1                     ! < reference height of meteo values

      REAL  (KIND=wp), PARAMETER           ::  &

         kappa = 0.4_wp

      ! ------------------------------------------------------------------------------
      ! Calculate transfer coefficient
      ! ------------------------------------------------------------------------------

      ! Assign a few value - FIXME: This should be the height of first atmospheric level
      z1 = 10.0_wp


      DO i=ivstart, ivend  ! FIXME: Method shoudl be choosable via namelist

         ! ------------------
         ! Transfer coefficient after Schloegl et al. 2017
         ! ------------------

         ! Needs to be implemented if needed


         ! ------------------
         ! Assume neutral conditions
         ! -----------------

         ! write(*,*) 'ivstart,ivend,z1,z0_sn', ivstart,ivend,z1,z0_sn
         tch_sn = (kappa*kappa) /  ( LOG(z1/z0_sn) * LOG(z1/z0_sn) )


      ENDDO




! =============================================================================
! - END subroutine: calculate_tch
! ============================================================================

   END SUBROUTINE calculate_tch



! =============================================================================
! + Begin subroutine: calculate_turbulent_fluxes
! ============================================================================

   SUBROUTINE calculate_turbulent_fluxes(nvec, ivstart, ivend   , &
      lhflx_sn, shflx_sn     , &
      t1, t0, q1, q0         , &
      tch_sn, zuv, ps)

      ! Subroutine arguments
      INTEGER                          :: &
      !
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend                  ! < end index for computations in the parallel program

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
      !
         q0                  , & ! < specific water vapour at snow surface
         lhflx_sn            , & ! < latent heat flux
         shflx_sn                ! < sensible heat flux


      REAL (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
      !
         t1                 , & ! < air temperature               - upper value
         t0                 , & ! < snow surface temperature      - lower value
         q1                 , & ! < specific water vapour content - upper value
         tch_sn             , & ! < transfer coefficent
         zuv                , & ! < wind speed
         ps

      ! Local variables
      INTEGER :: &
         i                      ! loop index in x-direction

      ! Local arrays, vectors and scalars
      REAL    (KIND = wp) ::  &

         E_s                   , & ! saturation water pressure
         e_v                   , & ! water vapour pressure

         t_v                   , & ! virtuell temperature

         rho_atm               , & ! density of atmosphere
         low_uv                    ! lower limit of windspeed

      ! ------------------------------------------------------------------------------
      ! + Calculate turbulent fluxes
      ! ------------------------------------------------------------------------------

      DO i=ivstart, ivend

         ! -------------------------
         ! + Some pre-calculations
         ! -------------------------

         ! Specific humidity at surface assuming saturation
         E_s = 6.112_wp * EXP( (22.46_wp * (t0(i)-273.15_wp)) / (272.62_wp + t0(i)) )  ! Saturation vapour pressure via Magnus Equation
         e_v = 100.0_wp * E_s                                          ! Water vapour pressure from relative humidty
         q0(i)  = 0.622 * (e_v/ps(i))

         ! Virtuell temperature
         t_v = t1(i) * (1.0_wp + rvd_m_o * q1(i))

         ! Density of atmosphere
         rho_atm = ps(i) / ( r_d * t_v)

         ! Apply lower limit for wind - i.e. there is always some exchange
         low_uv = MAX(0.1_wp, zuv(i))

         ! -------------------------
         ! + Latent heat
         ! -------------------------

         lhflx_sn(i) = tch_sn(i) * low_uv * rho_atm * lh_v * (q1(i) - q0(i))

         ! Limit latent flux
         lhflx_sn(i) = MIN(150.0_wp, MAX(-150.0_wp,lhflx_sn(i)) )

         ! -------------------------
         ! + Sensible heat
         ! -------------------------

         shflx_sn(i) = tch_sn(i) * low_uv * rho_atm * cp_d * (t1(i) - t0(i))

         ! Limit sensible heat flux
         shflx_sn(i) = MIN(200.0_wp, MAX(-200.0_wp,shflx_sn(i)) )

      END DO


! =============================================================================
! - END subroutine: calculate_turbulent_fluxes
! ============================================================================

   END SUBROUTINE calculate_turbulent_fluxes


! =============================================================================
! + Begin subroutine: calculate_radiative_fluxes
! ============================================================================

   SUBROUTINE calculate_radiative_fluxes(nvec, ivstart, ivend, ke_snow, top       , &
   &                                 swflx_sn_net, swflx_sn_dn, swflx_sn_up   , &
   &                                 swflx_sn_abs, lwflx_sn_up, lwflx_sn_dn   , &
   &                                 alpha_sn, t_sn_sfc                       , &
   &                                 rho_sn, dzm_sn, dt )




      ! Subroutine arguments
      INTEGER, INTENT(IN)                          :: &
      !
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow                ! < number of snow layers

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
         swflx_sn_net       , & ! < short-wave radiation flux - netto              [W/m**2]
         swflx_sn_dn        , & ! < short-wave radiation flux - downward           [W/m**2]
         swflx_sn_up        , & ! < short-wave radiation flux - upward             [W/m**2]

         lwflx_sn_dn        , & ! < long-wave radiation flux - downward            [W/m**2]
         lwflx_sn_up        , & ! < long-wave radiation flux - upward              [W/m**2]

         alpha_sn           , & ! < snow surface albedo                            [K]
         t_sn_sfc               ! < snow surface temperature                       [K]

      REAL (KIND = wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
         swflx_sn_abs           ! < absorbed short wave radiation flux             [W/m**2]


      REAL (KIND = wp), DIMENSION(nvec,ke_snow), INTENT(IN) :: &
         rho_sn             , & ! < snow layer density                             [kg/m3]
         dzm_sn                 ! < snow layer thickness                           [m]

      INTEGER, DIMENSION(nvec), INTENT(IN) :: &
         top                    ! top layer index


      real (kind = wp), intent(in)  ::  &
         dt                     ! time step

      ! Local variables
      INTEGER :: &
         i                  , & ! < loop index in x-direction
         ksn                    ! < loop index in z-direction


      REAL (KIND=wp), PARAMETER            :: &

         a1 = 0.90_wp       , &    ! coefficients for the albedo parameterization
         a0 = 0.60_wp


      REAL (KIND=wp) :: &

         k_ext                     ! extinction coeficient

      ! ------------------------------------------------------------------------------
      ! + Calculate radiative fluxes
      ! ------------------------------------------------------------------------------

      DO i=ivstart,ivend

         ! -------------------------
         ! + Calaculate albedo
         ! -------------------------
         ! IF(t_sn_sfc(i) - 273.15_wp .LE. -2.0_wp) THEN
         !    alpha_sn(i) = a1
         ! ELSE
         !    alpha_sn(i)  = a1 - 0.15_wp*(a1-a0)*((MIN(t_sn_sfc(i), t0_melt) - 273.15_wp) + 2.0_wp)
         ! ENDIF

         if( t_sn_sfc(i) .ge. 273.0_wp ) then
            alpha_sn(i) = a0 + (alpha_sn(i) - a0 ) * exp( (-1.0_wp/(200.0_wp*3600.0_wp)) * dt )
         else
            alpha_sn(i) = a0 + (alpha_sn(i) - a0 ) * exp( (-1.0_wp/(480.0_wp*3600.0_wp)) * dt )
         endif
         alpha_sn(i) = min(a1, max(a0, alpha_sn(i)))

         ! -------------------------
         ! + Calculate upward short-wave radiation
         ! -------------------------

         swflx_sn_up(i) = swflx_sn_dn(i) * alpha_sn(i)

         ! -------------------------
         ! + Calculate net. short-wave radiation
         ! -------------------------

         swflx_sn_net(i) = swflx_sn_dn(i) - swflx_sn_up(i)

         ! -------------------------
         ! + Calculate absprbed short-wave radiation
         ! -------------------------

         swflx_sn_abs(i,:) = 0.0_wp
         tbloop: DO ksn = top(i), 1, -1

            k_ext               = ( rho_sn(i,ksn) / 3.0_wp ) + 50.0_wp
            swflx_sn_abs(i,ksn) = swflx_sn_net(i) * (1.0_wp - exp(-k_ext * dzm_sn(i,ksn)))
            swflx_sn_net(i)     = swflx_sn_net(i) - swflx_sn_abs(i,ksn)

            IF( swflx_sn_net(i) .LE. 0.0_wp ) EXIT tbloop

         ENDDO tbloop


         ! Put the remaining energy into the bottom layer
         IF(swflx_sn_net(i) .GT.  0.0_wp) THEN
            swflx_sn_abs(i,1) = swflx_sn_abs(i,1) + swflx_sn_net(i)
         ENDIF

         ! Recalculate swflx_sn_net to have it available for output
         swflx_sn_net(i) = swflx_sn_dn(i) - swflx_sn_up(i)

         ! -------------------------
         ! + Calculate upward long-wave radiation
         ! -------------------------

         lwflx_sn_up(i) = stbo * (1.0_wp - ctalb) * t_sn_sfc(i)**4


      ENDDO

! =============================================================================
! - END subroutine: calculate_radiative_fluxes
! ============================================================================

   END SUBROUTINE calculate_radiative_fluxes


! =============================================================================
! + Begin subroutine: calculate_atmospheric_forcing
! ============================================================================

   SUBROUTINE calculate_atmospheric_forcing(nvec, ivstart, ivend, ke_snow, top  , &
   &                                    for_sn, swflx_sn_net, swflx_sn_abs  , &
   &                                    lwflx_sn_up, lwflx_sn_dn            , &
   &                                    lhflx_sn, shflx_sn)


      ! Subroutine arguments
      INTEGER, INTENT(IN)                          :: &
      !
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow                ! < number of snow layers

      INTEGER, DIMENSION(nvec), INTENT(IN) :: &
         top                    ! top layer index


      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
         for_sn             , & ! < total atmospheric forcing                      [W/m**2]
         swflx_sn_net       , & ! < short-wave radiation flux - netto              [W/m**2]

         lwflx_sn_dn        , & ! < long-wave radiation flux - downward            [W/m**2]
         lwflx_sn_up        , & ! < long-wave radiation flux - upward              [W/m**2]
         lhflx_sn           , & ! < latent heat flux                               [W/m**2]
         shflx_sn               ! < sensible heat flux                             [W/m**2]


      REAL (KIND = wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
         swflx_sn_abs           ! < absorbed short wave radiation flux             [W/m**2]


      ! Local variables
      INTEGER :: &
         i                      ! < loop index in x-direction


      ! ------------------------------------------------------------------------------
      ! + Calculate atmospheric forcing
      ! ------------------------------------------------------------------------------

      DO i=ivstart,ivend

         ! -------------------------
         ! + With absorption of short wave radiation
         ! -------------------------

         IF(top(i) .GE. 1) THEN

            for_sn(i)  = swflx_sn_abs(i,top(i)) + (lwflx_sn_dn(i) - lwflx_sn_up(i)) + lhflx_sn(i) + shflx_sn(i)

         ELSE

            ! -------------------------
            ! + Without absorption of short wave radiation
            ! -------------------------

            for_sn(i)  = swflx_sn_net(i) + (lwflx_sn_dn(i) - lwflx_sn_up(i)) + lhflx_sn(i) + shflx_sn(i)

         ENDIF

      ENDDO




! =============================================================================
! - END subroutine: calculate_atmospheric_forcing
! ============================================================================

   END SUBROUTINE calculate_atmospheric_forcing






!------------------------------------------------------------------------------
! End of module mo_nix_blanc
!------------------------------------------------------------------------------

END MODULE mo_nix_meteo_util




