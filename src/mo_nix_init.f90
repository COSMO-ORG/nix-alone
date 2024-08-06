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

   USE mo_nix_config,              ONLY: ke_snow          ,  &
                                         nvec             ,  &
                                         itype_nix_start  ,  &
                                         min_newsnow_layer,  &
                                         min_height_layer ,  &
                                         max_height_layer

   USE mo_physical_constants,      ONLY: t0_melt => tmelt    ! absolute zero for temperature

   USE mo_nix_constants,           ONLY: rho_i ! density of ice (kg/m^3)

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE :: nix_snow_analysis

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: nix_init


CONTAINS

! =============================================================================
! + Begin subroutine: nix_snow_analysis
! ============================================================================

SUBROUTINE nix_snow_analysis(ke_snow   , &
   &         t_sn      , &
   &         theta_i   , & ! volumetric ice content
   &         theta_w   , & ! volumetric water content
   &         theta_a   , & ! volumetric air content
   &         dzm_sn    , & ! snow layer thickness
   &         hn_sn     , & ! new snow amount
   &         top_sn    , & ! index of first (top) snow layer
   &         h_snow    , & ! snow depth in NIX (m)
   &         h_snow_in , & ! target snow depth (m), coming from snow analysis
   &         ta          ) ! reference air temperature to determine melting conditions


      ! Subroutine Arguments
      INTEGER, INTENT(IN)    :: &
         ke_snow         ! number of snow layers

      REAL (KIND = wp), DIMENSION(ke_snow), INTENT(INOUT) :: &
         t_sn                , & ! snow layer temperature (main level)
         theta_i             , & ! volumetric ice content
         theta_w             , & ! volumetric water content
         theta_a             , & ! volumetric air content
         dzm_sn                  ! snow layer thickness

      REAL(KIND = wp), INTENT(INOUT) :: &
         hn_sn                   ! new snow amount

      INTEGER, INTENT(INOUT) :: &
         top_sn                  ! index of first (top) snow layer

      REAL    (KIND = wp), INTENT(INOUT) :: &
         h_snow          ! snow depth

      REAL    (KIND = wp), INTENT(IN) :: &
         h_snow_in        , &   ! snow depth
         ta                     ! reference temperature to determine melting conditions

      ! Local variables
      INTEGER :: &
         i            , &     ! loop index in x-direction
         ksn          , &     ! loop index in y-direction
         top_sn_new

      REAL (KIND = wp), PARAMETER :: &
         corr_max_hs = 0.2_wp , & ! maximum snow depth to be corrected.
         rho_wet = 250.0_wp   , & ! typical wet snow density
         rho_dry = 100.0_wp       ! typical dry snow density

      REAL (KIND = wp)       :: &
         corr                     ! correction factor

      ! Case 1: enough snow or no snow in both analysis and model --> nothing to do, simply return
      IF (     (h_snow_in .GT. corr_max_hs .AND. h_snow .GT. corr_max_hs)   &
          .OR. (h_snow_in .EQ. 0.0_wp      .AND. h_snow .EQ. 0.0_wp)) THEN
         RETURN
      ENDIF

      ! Case 2: no snow in snow analysis, snow in NIX --> remove all snow
      IF (h_snow_in .LT. min_newsnow_layer) THEN
         h_snow = 0.0_wp      ! clear all snow
         top_sn = 0           ! reset top level index
         hn_sn  = 0.0_wp      ! reset new snow amounts - storage

         ! Reset snow profiles
         DO ksn = 1, ke_snow
            dzm_sn(ksn)  = 0.0_wp
            t_sn(ksn)    = 0.0_wp
            theta_i(ksn) = 0.0_wp
            theta_w(ksn) = 0.0_wp
            theta_a(ksn) = 0.0_wp
         ENDDO

         RETURN
      ENDIF

      ! Case 3: no snow in NIX, snow in snow analysis --> build snowcover
      IF (h_snow .EQ. 0.0_wp .AND. h_snow_in .LT. min_newsnow_layer) THEN
         top_sn = MIN(ke_snow, FLOOR(h_snow_in / min_height_layer))
         hn_sn  = 0.0_wp      ! reset new snow amounts - storage

         ! Assume snow profiles
         h_snow = 0.0_wp
         DO ksn = 1, top_sn
            dzm_sn(ksn)  = h_snow / top_sn
            h_snow       = h_snow + dzm_sn(ksn)
            t_sn(ksn)    = MIN(t0_melt, ta)     ! Assume snow temperature = air temperature, while not exceeding melting point
            IF (ta .GE. t0_melt) THEN
               ! Assume melting conditions
               theta_i(ksn) = rho_wet / rho_i
               theta_w(ksn) = 0.04_wp
            ELSE
               ! Assume dry conditions
               theta_i(ksn) = rho_dry / rho_i
               theta_w(ksn) = 0.0_wp
            ENDIF
            theta_a(ksn) = 1.0_wp - theta_i(ksn) - theta_w(ksn)
         ENDDO

         RETURN

      ENDIF


      ! Case 4: more snow in analysis than in NIX
      IF (h_snow_in .GT. h_snow) THEN
         ! Determine new number of layers: not less than already present, not more than allowed
         top_sn_new = MAX(top_sn, MIN(ke_snow, FLOOR(h_snow_in / max_height_layer)))
         DO ksn = top_sn_new, 1, -1                                ! Need to loop in reverse, to avoid overwriting elements we still need
            i            = NINT((real(top_sn) * (real(ksn-1) / real(top_sn_new)) + .5_wp))   ! Determine layer mapping
            dzm_sn(ksn)  = dzm_sn(i) * (h_snow_in / h_snow) * (real(top_sn) / real(top_sn_new))
            t_sn(ksn)    = t_sn(i)
            theta_i(ksn) = theta_i(i)
            theta_w(ksn) = theta_w(i)
            theta_a(ksn) = 1.0_wp - theta_i(ksn) - theta_w(ksn)
         ENDDO
         top_sn = top_sn_new  ! assign new top level index
         hn_sn  = 0.0_wp      ! reset new snow amounts - storage
         RETURN
      ENDIF


      ! Case 5: less snow in analysis than in NIX, simply scale all layers
      corr = (h_snow_in / h_snow)
      IF (h_snow_in .LT. h_snow) THEN
         h_snow = 0.0_wp
         DO ksn = 1, top_sn
            dzm_sn(ksn)  = dzm_sn(ksn) * corr
            h_snow       = h_snow + dzm_sn(ksn)
         ENDDO
         hn_sn  = 0.0_wp      ! reset new snow amounts - storage
         RETURN
      ENDIF

! =============================================================================
! - END subroutine: nix_snow_analysis
! ============================================================================

   END SUBROUTINE nix_snow_analysis



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
   &         h_snow    , &
   &         ta          ) ! reference temperature to determine melting conditions or not



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
         h_snow               ! snow depth

      REAL    (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
         ta               ! reference temperature to determine melting conditions or not

         ! Local variables
      INTEGER :: &
         i            , &     ! loop index in x-direction
         ksn                  ! loop index in y-direction

      REAL(KIND = wp) :: &
         h_snow_in

      IF (itype_nix_start .EQ. 1) THEN

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

      ELSE IF (itype_nix_start .EQ. 3 .OR. itype_nix_start .EQ. 4) THEN
         ! ------------------------------------------------------------------------------
         ! Section 3 - Warm start - required nix fields are taken from initial condition
         !         4 - Warm start + snow analysis
         ! ------------------------------------------------------------------------------

         !Set snow height and top index to zero
         DO i = ivstart, ivend

            top_sn(i) = 0
            h_snow_in = h_snow(i)   ! Store input h_snow to apply snow analysis
            h_snow(i) = 0.0_wp
            hn_sn(i)  = 0.0_wp      ! reset new snow amounts - storage

            DO ksn = 1, ke_snow

               IF (dzm_sn(i,ksn) .GT. 0) THEN
                  top_sn(i) = ksn                            ! update top level index
                  h_snow(i) = h_snow(i) + dzm_sn(i,ksn)      ! update snow depth

               ELSE
                  ! Reset remainder of snow profile
                  dzm_sn(i,ksn)  = 0.0_wp
                  t_sn(i,ksn)    = 0.0_wp
                  theta_i(i,ksn) = 0.0_wp
                  theta_w(i,ksn) = 0.0_wp
                  theta_a(i,ksn) = 0.0_wp

               ENDIF

            ENDDO

            IF (itype_nix_start .EQ. 4) THEN

               CALL nix_snow_analysis(ke_snow   , &
               &         t_sn(i,:)      , &
               &         theta_i(i,:)   , & ! volumetric ice content
               &         theta_w(i,:)   , & ! volumetric water content
               &         theta_a(i,:)   , & ! volumetric air content
               &         dzm_sn(i,:)    , & ! snow layer thickness
               &         hn_sn(i)       , & ! new snow amount
               &         top_sn(i)      , & ! index of first (top) snow layer
               &         h_snow(i)      , &
               &         h_snow_in      , &
               &         ta(i)            ) ! snow depth                        (m H2O))
            ELSE
               h_snow_in = h_snow(i)
            ENDIF

         ENDDO

      ELSE

         WRITE (0,*) "ERROR: unknown itype_nix_start: ", itype_nix_start
         CALL EXIT (1)

      ENDIF


! =============================================================================
! - END subroutine: init_nix
! ============================================================================

   END SUBROUTINE nix_init




!------------------------------------------------------------------------------
! End of module mo_nix_init
!------------------------------------------------------------------------------

END MODULE mo_nix_init




