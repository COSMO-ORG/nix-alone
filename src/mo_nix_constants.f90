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
! Begin of module mo_nix_constants
! ------------------------------------------------------------------------------

MODULE mo_nix_constants

   USE mo_kind,                    ONLY: wp

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC      ! All constants and variables in this module are public





   ! 1. Physical constants  NOTE: a bunch of these already exists globally, but
   !                              might have slithly different values so double check and replace at some point
   ! -------------------------------------------------

   REAL  (KIND=wp), PARAMETER           ::  &

      rho_i                       = 917.0_wp               , & ! density of ice                 (kg/m**3)
      rho_a                       = 1.1_wp                 , & ! density of air                 (kg/m**3)

      specific_heat_air           = 1004.67_wp             , & !
      specific_heat_ice           = 2100.0_wp              , & !
      specific_heat_water         = 4190.0_wp              , & !

      theta_s                     = 1.0_wp                 , & ! Saturated Water Content, for now we say 1.0
      theta_r                     = 0.0_wp                 , & ! Minimum amount of liquid water that will remain
      ctalb                       = 0.004_wp               , & ! VS : added snow physical constants

      e_snow                      = 0.98_wp                    ! lw emissivity for snow

   REAL  (KIND=wp), PARAMETER :: eps_div  = 1.0E-6_wp

   INTEGER , PARAMETER :: itype_heatcond = 1
!------------------------------------------------------------------------------
! End of module mo_nix_blanc
!------------------------------------------------------------------------------

END MODULE mo_nix_constants




