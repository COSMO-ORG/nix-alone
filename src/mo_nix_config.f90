!>
!! @brief configuration setup for Swiss multi-layer snow scheme NIX
!!
!! @author Sascha Bellaire, MCH (2022-07-21)
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nix_config

   USE mo_kind,               ONLY: wp
!   USE mo_io_units,           ONLY: filename_max
!   USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
!   USE mo_time_config,        ONLY: time_config
!   USE mtime,                 ONLY: timedelta, event, getPTStringFromMS, MAX_TIMEDELTA_STR_LEN, MAX_DATETIME_STR_LEN, &
!   &                              newTimedelta, newEvent

   IMPLICIT NONE

   PUBLIC                          ! All switches and variables in this module are public

   !<--------------------------------------------------------------------------
   !
   ! Basic configuration setup for the Swiss multi-layer snow cover scheme (NIX)
   !
   !>--------------------------------------------------------------------------

   ! namelist variables

   LOGICAL :: lnix        ! > run with the Swiss multi-layer snow model NIX

   INTEGER, PARAMETER ::  &
      ke_snow = 10      , & ! > run with this amount of snow layers
      ke_soil = 1       , & ! > run with this amount of soil layers
      nvec = 1

   REAL (KIND = wp), PARAMETER ::   &
      zdt = 900.0_wp              , & ! model timestep : note this should correspond to the meteotimestep     (s)
      min_height_layer = 0.002_wp , & ! minimum layer thickness
      max_height_layer = 0.01_wp      ! maximum layer thickness

   REAL (KIND = wp) :: &

      z0_sn  = 0.002_wp                ! roughness length

CONTAINS

   !>-----------------------------------------------------------------------------------
   !
   ! Further configurations of the Swiss multi-layer snow scheme NIX
   !
   !         - Note: Currently non required
   !
   !<-----------------------------------------------------------------------------------

   SUBROUTINE configure_nix()

      ! Fill configure subroutine if needed

   END SUBROUTINE configure_nix


END MODULE mo_nix_config
