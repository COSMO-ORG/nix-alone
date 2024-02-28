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
! Begin of module mo_nix_init
! ------------------------------------------------------------------------------

MODULE mo_nix_init

   USE mo_kind,                    ONLY: wp

   USE mo_nix_config,              ONLY: ke_snow, nvec

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: nix_init


CONTAINS

! =============================================================================
! + Begin subroutine: init_nix
! ============================================================================

   SUBROUTINE nix_init(nvec      , &
   &         ke_snow   , &
   &         ivstart   , &
   &         ivend     , &
   &         t_sn      , &
   &         theta_i   , & ! volumetric ice content
   &         theta_w   , & ! volumetric water content
   &         theta_a   , & ! volumetric air content
   &         dzm_sn    , & ! snow layer thickness
   &         hn_sn     , & ! new snow amount
   &         top_sn    , & ! index of first (top) snow layer
   &         h_snow      ) ! snow depth                        (m H2O)



      ! Subroutine Arguments
      INTEGER, INTENT(IN)    :: &
         nvec        , & ! array dimensions
         ke_snow     , & ! number of snow layers
         ivstart     , & ! start index for computations in the parallel program
         ivend           ! end index for computations in the parallel program

      REAL (KIND = wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
         t_sn                , & ! snow layer temperature (main level)
         theta_i             , & ! volumetric ice content
         theta_w             , & ! volumetric water content
         theta_a             , & ! volumetric air content
         dzm_sn                  ! snow layer thickness

      REAL(KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
         hn_sn                   ! new snow amount

      INTEGER, DIMENSION(nvec), INTENT(INOUT) :: &
         top_sn                  ! index of first (top) snow layer

      REAL    (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
         h_snow          ! snow depth

      ! Local variables
      INTEGER :: &
         i            , &     ! loop index in x-direction
         ksn                  ! loop index in y-direction

      ! ------------------------------------------------------------------------------
      ! Section 1 - Hard cold start - All snow is wiped out.
      ! ------------------------------------------------------------------------------

      !Set snow height and top index to zero
      DO i = ivstart, ivend

         h_snow(i) = 0.0_wp      ! clear all snow
         top_sn(i) = 0           ! rest top level index

         hn_sn(i)  = 0.0_wp      ! reset new snow amounts - storage

      END DO

      ! Reset snow profiles
      DO ksn = 1, ke_snow
         DO i = ivstart, ivend

            dzm_sn(i,ksn)  = 0.0_wp
            t_sn(i,ksn)    = 0.0_wp
            theta_i(i,ksn) = 0.0_wp
            theta_w(i,ksn) = 0.0_wp
            theta_a(i,ksn) = 0.0_wp

         ENDDO
      ENDDO


! =============================================================================
! - END subroutine: init_nix
! ============================================================================

   END SUBROUTINE nix_init




!------------------------------------------------------------------------------
! End of module mo_nix_blanc
!------------------------------------------------------------------------------

END MODULE mo_nix_init




