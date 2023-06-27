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
! Begin of module mo_nix_snow_util
! ------------------------------------------------------------------------------

MODULE mo_nix_snow_util

   USE mo_kind,                    ONLY: wp

   USE mo_nix_constants,           ONLY: eps_div

   USE mo_physical_constants,      ONLY: t0_melt => tmelt   , &  ! absolute zero for temperature
      rho_w   => rhoh2o       ! density if liquid water

   USE mo_nix_constants,           ONLY: rho_i               , & ! density of ice
      rho_a               , & ! density of air
      specific_heat_air   , &
      specific_heat_ice   , &
      specific_heat_water


! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: calc_hn_density
   PUBLIC :: calc_hn_amounts
   PUBLIC :: update_nix_state

CONTAINS


! =============================================================================
! + Begin subroutine: calc_hn_density
! ============================================================================

   SUBROUTINE calc_hn_density(nvec, ivstart, ivend, rho_hn, uv, t_a)

      ! Subroutine arguments
      INTEGER                          :: &
      !
         nvec              , & ! < array dimensions
         ivstart           , & ! < start index for computations in the parallel program
         ivend                 ! < end index for computations in the parallel program

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
      !
         rho_hn                ! < new snow density

      REAL (KIND = wp), DIMENSION(nvec) :: &
      !
         uv                 , & ! < wind speed
         t_a                    ! < air temperature first (from surface) atmospheric level


      ! Local variables
      INTEGER :: &
         i                    ! loop index in x-direction


      REAL    (KIND = wp), PARAMETER ::  &

         rh     = 0.8        , & ! relative humidity                          (-)

         beta01 = 3.28       , & ! coefficients
         beta1  = 0.03       , & !
         beta02 = -0.36      , &
         beta2  = -0.75      , &
         beta3  = 0.3

      REAL    (KIND = wp) ::  &

         vw                  , & ! wind speed limited to 2 m s-1                  (m/s)
         t_c                 , & ! air temperature converted to degrees Celsius   (Celsius)
         arg                 ! argument for power low

      REAL    (KIND = wp), PARAMETER ::  &

         rho_hn_min = 100.0_wp     , &
         rho_hn_max = 150.0_wp

      INTEGER ::                              &
         itype_snow_density = 2


      ! ------------------------------------------------------------------------------
      ! Calculate new snow density -  FIXME: Switch hardcoded could become a namelist switch itype_snow_density
      ! ------------------------------------------------------------------------------

      DO i = ivstart, ivend

         ! ---------------------
         ! SNOWPACK
         ! ---------------------

         IF(itype_snow_density == 1) THEN



         ENDIF


         ! ---------------------
         ! TERRA
         ! ---------------------

         IF(itype_snow_density == 2) THEN

            rho_hn(i) = rho_hn_min + (rho_hn_max - rho_hn_min) * ((t_a(i)-258.15_wp) / (t0_melt-258.15_wp))**2

         ENDIF

         ! ---------------------
         ! Fixed new snow density - Note: Switch hardcoded could become a namelist switch itype_snow_density
         ! ---------------------

         IF(itype_snow_density == 3) THEN

            rho_hn(i) = 100.0_wp

         ENDIF


         ! ---------------------
         ! Apply limits
         ! ---------------------

         rho_hn(i) = MAX(rho_hn_min, MIN(rho_hn(i), rho_hn_max))



      END DO ! end of i



! =============================================================================
! - END subroutine: calc_hn_density
! ============================================================================

   END SUBROUTINE calc_hn_density


! =============================================================================
! + Begin subroutine: calc_hn_amounts
! ============================================================================

   SUBROUTINE calc_hn_amounts(nvec, ivstart, ivend, hn_sn, zsnow_rate, dt)

      ! Subroutine arguments
      INTEGER                          :: &
      !
         nvec              , & ! < array dimensions
         ivstart           , & ! < start index for computations in the parallel program
         ivend                 ! < end index for computations in the parallel program

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
      !
         hn_sn                 ! < new snow density


      REAL (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
      !
         zsnow_rate             ! < rate of snow fall

      REAL (KIND = wp), INTENT(IN) :: &
         dt                     ! < time step

      ! Local variables
      INTEGER :: &
         i                    ! loop index in x-direction


      ! ------------------------------------------------------------------------------
      ! Calculate new snow amounts - NOTE: Don't allow negative values of zsnow_rate
      ! ------------------------------------------------------------------------------

      DO i = ivstart, ivend

         hn_sn(i) = hn_sn(i) + ( MAX(0.0_wp, zsnow_rate(i)) * dt )

      END DO

! =============================================================================
! - END subroutine: calc_hn_amounts
! ============================================================================

   END SUBROUTINE calc_hn_amounts


! =============================================================================
! + Begin subroutine: update_nix_state
! ============================================================================

   SUBROUTINE update_nix_state(nvec, ivstart, ivend, top, ke_snow , &
   &                   dzm_sn, rho_sn                     , &
   &                   theta_i, theta_w, theta_a          , &
   &                   t_sn, t_sn_n, hm_sn, zm_sn                 , &
   &                   hcap_sn, hcon_sn,  mass_sn         , &
   &                   h_snow, t_sn_sfc )

      ! Subroutine Arguments
      INTEGER, INTENT(IN)   :: &
         nvec    , &
         ivstart , &
         ivend

      INTEGER, DIMENSION(nvec), INTENT(IN) :: &
         top                 ! index of the first (top) layer index       (-)

      INTEGER, INTENT(IN)  ::  &
         ke_snow             ! number of snow layers

      REAL (KIND = wp)    , DIMENSION(nvec,ke_snow),  INTENT(INOUT) :: &
         dzm_sn        , &   ! layer thickness between main levels        (m)
         rho_sn        , &   ! snow layer density                         (kg/m3)
         theta_i       , &   ! volumetric ice content                     (-)
         theta_w       , &   ! water ice content                          (-)
         theta_a       , &   ! air ice content                            (-)
         t_sn          , &      ! snow temperature (main level)              (K)
         hm_sn         , &   ! height of snow main levels                 (m)
         zm_sn         , &   ! depth of snow main levels                  (m)
         hcap_sn       , &   ! heat capacity                              ( )
         hcon_sn       , &   ! heat conductivity                          ( )
         mass_sn             ! snow layer mass                            (kg)

      REAL (KIND = wp)    , DIMENSION(nvec,ke_snow+1),  INTENT(INOUT) :: &
         t_sn_n ! snow temperatures nodal

      REAL (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
         h_snow        , &      ! snow height                                (m)
         t_sn_sfc            ! snow surface temperature                   (K)

      ! Local variables
      INTEGER  :: &
         i       , &
         ksn                 ! loop index over snow layers (vertical)

      ! ------------------------------------------------------------------------------
      ! Update nix dependent state variables
      ! ------------------------------------------------------------------------------


      DO i = ivstart, ivend

         IF(top(i) .GE. 1) THEN

            ! -------------------------
            ! Height of snow (main) levels
            ! -------------------------

            DO ksn=1,top(i)

               IF(ksn .EQ. 1) THEN
                  hm_sn(i,ksn) = dzm_sn(i,ksn)
               ELSE
                  hm_sn(i,ksn) = hm_sn(i,ksn-1) + dzm_sn(i,ksn)
               ENDIF

            ENDDO ! ksn


            ! -------------------------
            ! Snow Height
            ! -------------------------

            h_snow(i) = hm_sn(i,top(i))

            ! --------------------------
            ! Depth of snow layer (main) levels
            ! -------------------------

            DO ksn = top(i),1,-1
               zm_sn(i,(top(i)+1) - ksn ) = hm_sn(i,ksn) !invert height vector
            ENDDO

            ! --------------------------
            ! Volumetric air content
            ! -------------------------

            DO ksn = top(i),1,-1
               theta_a(i,ksn) = max(0.0_wp,1.0_wp - theta_i(i,ksn) - theta_w(i,ksn))
            ENDDO

            ! --------------------------
            ! Snow layer density
            ! -------------------------

            DO ksn = 1, top(i), 1

               IF(theta_i(i,ksn) .EQ. 0.0_wp) THEN
                  rho_sn(i,ksn) = 0.0_wp
               ELSE
                  rho_sn(i,ksn) = theta_i(i,ksn)*rho_i + theta_w(i,ksn)*rho_w
               ENDIF

            ENDDO

            ! FIXME: Why are we doing this again, but the other way round???
            DO ksn = top(i), 1, -1

               IF(theta_i(i,ksn) .EQ. 0.0_wp) THEN
                  rho_sn(i,ksn) = 0.0_wp
               ELSE
                  rho_sn(i,ksn) = theta_i(i,ksn)*rho_i + theta_w(i,ksn)*rho_w
               ENDIF

            ENDDO

            ! --------------------------
            ! Heat capacity
            ! --------------------------

            DO ksn = 1, top(i), 1

               IF(rho_sn(i,ksn) .LT. eps_div) THEN
                  hcap_sn(i,ksn) = 0.0_wp
               ELSE
                  hcap_sn(i,ksn) = (  rho_a   * theta_a(i,ksn) * specific_heat_air     &
                     + rho_i   * theta_i(i,ksn) * specific_heat_ice     &
                     + rho_w   * theta_w(i,ksn) * specific_heat_water)  &
                     / rho_sn(i,ksn)
               ENDIF

            ENDDO

            ! --------------------------
            ! Heat conductivity
            ! --------------------------

            DO ksn = 1, top(i), 1

               IF(rho_sn(i,ksn) .LT. eps_div) THEN
                  hcon_sn(i,ksn) = 0.0_wp
               ELSE
                  hcon_sn(i,ksn) = 2.22_wp * EXP(1.88_wp * LOG(rho_sn(i,ksn)/rho_i))
               ENDIF

            ENDDO

            ! --------------------------
            ! Snow layer mass
            ! --------------------------

            DO ksn = 1, top(i), 1

               IF(ksn .EQ. 1) THEN
                  mass_sn(i,ksn) = zm_sn(i,ksn) * ((theta_i(i,ksn) * rho_i) + (theta_w(i,ksn) * rho_w))
               ELSE
                  mass_sn(i,ksn) = ABS(zm_sn(i,ksn) - zm_sn(i,ksn-1))  * ((theta_i(i,ksn) * rho_i) + (theta_w(i,ksn) * rho_w))
               ENDIF

            ENDDO


            ! --------------------------
            ! Snow surface temperature
            ! --------------------------

            ! Update snow surface temperature ! Now with nodal temperatures
            t_sn_sfc(i) = t_sn_n(i,top(i)+1)


         ENDIF

      ENDDO



! =============================================================================
! - END subroutine: update_nix_state
! ============================================================================

   END SUBROUTINE update_nix_state





!------------------------------------------------------------------------------
! End of module mo_nix_snow_util
!------------------------------------------------------------------------------

END MODULE mo_nix_snow_util




