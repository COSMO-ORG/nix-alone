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
! Begin of module mo_nix_misc_util
! ------------------------------------------------------------------------------

MODULE mo_nix_misc_util

   USE mo_kind,                    ONLY: wp

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: save_nix_state
   PUBLIC :: update_lnd_state


CONTAINS

! =============================================================================
! + Begin subroutine: save_nix_state
! ============================================================================

   SUBROUTINE save_nix_state(nvec, ivstart, ivend, ke_soil, ke_snow  , &
   &                     zm_sn_old, theta_i_old, theta_w_old     , &
   &                     zm_sn    , theta_i    , theta_w )


      ! Subroutine arguments
      INTEGER, INTENT(IN)                          :: &
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow            , & ! < number of snow layers
         ke_soil

      REAL (KIND=wp), DIMENSION(nvec,ke_snow), INTENT(OUT) :: &
         zm_sn_old          , &
         theta_i_old        , &
         theta_w_old

      REAL (KIND=wp), DIMENSION(nvec,ke_snow), INTENT(IN) :: &
         zm_sn              , &
         theta_i            , &
         theta_w

      ! Local variables
      INTEGER      :: &

         i                  , & ! loop index in x-direction
         ksn                    ! loop index in z-direction (snow layers)


      ! ------------------------------------------------------------------------------
      ! + Save NIX state
      ! ------------------------------------------------------------------------------

      DO ksn = 1, ke_snow, 1
         DO i = ivstart, ivend

            zm_sn_old(i,ksn)   = zm_sn(i,ksn)

            theta_i_old(i,ksn) = theta_i(i,ksn)
            theta_w_old(i,ksn) = theta_w(i,ksn)

         ENDDO
      ENDDO



! =============================================================================
! - END subroutine: save_nix_state
! ============================================================================

   END SUBROUTINE save_nix_state



! =============================================================================
! + Begin subroutine: update_ln_state
! ============================================================================

   SUBROUTINE update_lnd_state(               &
   &                      nvec        , &
   &                      ivstart     , &
   &                      ivend       , &
   &                      ke_soil     , &
   &                      ke_snow     , &
   &                      top         , & ! < NIX  : index of first (top) snow layer            (  -  )
!
   &                      t_snow_new  , & ! < TERRA: temperature of the snow-surface            (  K  )
   &                      t_sn_sfc    , & ! < NIX  : temperature of the snow-surface            (  K  )
!
   &                      qv_s        , & ! < TERRA: specific humidity at the surface           (kg/kg)
   &                      qv0_sn      , & ! < NIX  : specific humidity at the surface           (kg/kg)
!
   &                      w_snow_new  , & ! < TERRA: water content of snow                       (m H2O)
   &                      swe_sn      , & ! < NIX  : water content of snow                       (m H20)
!
   &                      h_snow      , & ! < TERRA: snow height                                 (  m  )
   &                      hm_sn       , & ! < NIX  : height of snow layers (main level)          (  m  )
!
   &                      tch         , & ! < TERRA: turbulent transfer coefficient for heat     ( -- )
   &                      tcm         , & ! < TERRA: turbulent transfer coefficient for momentum ( -- )
   &                      tch_sn      , & ! < NIX  : turbulent transfer coefficient              ( -- )
!
   &                      runoff_s    , & ! < TERRA: surface water runoff                        (kg/m2)
   &                      runoff_sn   , & ! < NIX  : runoff (melt-water) through bottom layer    (kg/m2)
!
   &                      zshfl_snow  , & ! < TERRA: sensible heat flux snow/air interface       (W/m2)
   &                      shflx_sn    , & ! < NIX  : sensible heat flux snow/air interface       (W/m2)
!
   &                      zlhfl_snow  , & ! < TERRA: latent heat flux snow/air interface         (W/m2)
   &                      lhflx_sn      ) ! < NIX  : latent heat flux snow/air interface         (W/m2)


      ! Subroutine arguments
      INTEGER, INTENT(IN)                          :: &
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow            , & ! < number of snow layers
         ke_soil

      INTEGER, DIMENSION(nvec), INTENT(IN) :: &
         top                    ! < index of first (top) snow layer             - NIX


      REAL (KIND = wp), DIMENSION(nvec), INTENT(OUT) :: &
         t_snow_new         , & ! < temperature of the snow surface             - TERRA
         t_sn_sfc           , & ! < temperature of the snow surface             - NIX
         qv_s               , & ! < specific humidity at the surface            - TERRA
         qv0_sn             , & ! < specific humidity at the surface            - NIX
         w_snow_new         , & ! < water content of snow                       - TERRA
         swe_sn             , & ! < water content of snow                       - NIX
         h_snow             , & ! < snow height                                 - TERRA
         tch                , & ! < turbulent transfer coefficient for heat     - TERRA
         tcm                , & ! < turbulent transfer coefficient for momentum - TERRA
         tch_sn             , & ! < turbulent transfer coefficient              - NIX
         runoff_s           , & ! < surface water runoff                        - TERRA
         runoff_sn          , & ! < runoff (melt-water) through bottom layer    - NIX
         zshfl_snow         , & ! < sensible heat flux snow/air interface       - TERRA
         shflx_sn           , & ! < sensible heat flux snow/air interface       - NIX
         zlhfl_snow         , & ! < latent heat flux snow/air interface         - TERRA
         lhflx_sn               ! < latent heat flux snow/air interface         - NIX

      REAL (KIND = wp), DIMENSION(nvec, ke_snow), INTENT(OUT) :: &
         hm_sn                  ! < height of snow layers (main level) - NIX

      ! Local variables
      INTEGER :: &
         i, ksn, kso

      ! ------------------------------------------------------------------------------
      ! + Update prognostic land state from TERRA with NIX state variables.
      ! ------------------------------------------------------------------------------

      ! -----------------------
      ! 2D fields - one to one i.e. no processing required
      ! -----------------------

      DO i = ivstart, ivend

         IF(top(i) .GE. 1) THEN

            t_snow_new(i)  = t_sn_sfc(i)                 ! temperature of snow surface
            w_snow_new(i)  = swe_sn(i)/1000.0_wp         ! snow water equivialent FIXME: Add hn_sn (storage) to swe
            h_snow(i)      = hm_sn(i,top(i))             ! snow height
            zshfl_snow(i)  = shflx_sn(i)                 ! sensible heat flux
            zlhfl_snow(i)  = lhflx_sn(i)                 ! latent heat flux

         ENDIF


!      tch(i)         = tch_sn(i)                  ! transfercoeffieceint for heat     FIXME: Let's not do this yet not sure what will happen!
!      tcm(i)         = tch_sn(i)                  ! transfer coefficient for momentum FIXME: equals heat! ok?

!      qv_s(i)       = qv0_sn(i)                   ! specific humidity at surface FIXME: not sure if that is correct to just overwrite it
         !                                     compare calculation in sfc_terra.90
      ENDDO

      ! -----------------------
      ! 2D fields -  processing required
      ! -----------------------

      DO i = ivstart, ivend

         ! Surface runoff FIXME: also here not sure if we can just add or overwrite it
!      runoff_s(i)   = runoff_s(i) + runoff_sn(i)  ! surface runoff

      ENDDO




! =============================================================================
! - END subroutine: update_lnd_state
! ============================================================================

   END SUBROUTINE update_lnd_state





!------------------------------------------------------------------------------
! End of module mo_nix_misc_util
!------------------------------------------------------------------------------

END MODULE mo_nix_misc_util





