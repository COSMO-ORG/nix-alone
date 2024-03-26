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
! Begin of module mo_nix_stratigraphy
! ------------------------------------------------------------------------------

MODULE mo_nix_stratigraphy

   USE mo_kind,                    ONLY: wp

   USE mo_physical_constants,      ONLY: t0_melt => tmelt   , &  ! absolute zero for temperature
      rho_w   => rhoh2o       ! density if liquid water

   USE mo_nix_constants,             ONLY: eps_div

   USE mo_nix_constants,           ONLY: rho_i
   USE mo_nix_config,              ONLY: min_height_layer, max_height_layer


! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: nix_stratigraphy
   PUBLIC :: split_bottom_layer

CONTAINS

! =============================================================================
! + Begin subroutine: nix_stratigraphy
! ============================================================================

   SUBROUTINE nix_stratigraphy(nvec, ivstart, ivend           , &
   &                 top, ke_snow, dzm_sn, rho_sn   , &
   &                 theta_i, theta_w, theta_a      , &
   &                 t_sn, t_sn_n, hn_sn, rho_hn, t , alpha_sn)


      ! Subroutine Arguments
      INTEGER, INTENT(IN)   :: &
         nvec          , & ! < array dimensions
         ivstart       , & ! < start index for computations in the parallel program
         ivend             ! < end index for computations in the parallel program

      INTEGER, DIMENSION(nvec), INTENT(INOUT) :: &
         top                 ! index of the first (top) layer index       (-)

      INTEGER, INTENT(IN)  ::  &
         ke_snow             ! number of snow layers

      REAL    (KIND = wp), DIMENSION(nvec,1:ke_snow),  INTENT(INOUT) :: &
         dzm_sn        , &   ! layer thickness between main levels        (m)
         rho_sn        , &   ! snow layer density                         (kg/m3)
         theta_i       , &   ! volumetric ice content                     (-)
         theta_w       , &   ! water ice content                          (-)
         theta_a       , &   ! air ice content                            (-)
         t_sn                ! snow temperature (main level)              (K)

      REAL    (KIND = wp), DIMENSION(nvec,1:ke_snow+1),  INTENT(INOUT) :: &
         t_sn_n

      REAL    (KIND = wp), DIMENSION(nvec),  INTENT(INOUT) :: &
         hn_sn       , &        ! new snow amount                            (m)
         alpha_sn


      REAL    (KIND = wp), DIMENSION(nvec),  INTENT(IN) :: &
         rho_hn        , &   ! new snow density                           (kg/m3)
         t                ! air temperature first level                (K)

      ! Local variables

      INTEGER ::           &

         i            , &   ! loop index in x-direction
         ksn                ! loop index in z-direction (snow layers)

      REAL    (KIND = wp), DIMENSION(1:ke_snow) :: &
         tmp_sn             ! temporary utility array used for copying

      REAL    (KIND = wp), DIMENSION(1:ke_snow+1) :: &
         tmp_sn_n             ! temporary utility array used for copying
      ! ------------------------------------------------------------------------------
      ! Begin subroutine
      ! ------------------------------------------------------------------------------

      DO i = ivstart, ivend

         IF( (hn_sn(i)+eps_div)/rho_hn(i) .GE. max_height_layer) THEN  ! There is enough new snow for a full layer ...

            IF (top(i) .LT. ke_snow) THEN ! ... and maximum number of layers is not reached ...

               ! -----------------------
               ! Add layers
               ! -----------------------

               ! Update layer index
               top(i) = top(i) + 1  ! new index of top layer

               ! Assign properties
               dzm_sn(i,top(i))  = hn_sn(i)/rho_hn(i)
               rho_sn(i,top(i))  = rho_hn(i)
               theta_i(i,top(i)) = rho_hn(i) / rho_i
               theta_w(i,top(i)) = 0.0_wp
               theta_a(i,top(i)) = 1.0_wp - theta_i(i,top(i)) - theta_w(i,top(i))
               IF (top(i) .EQ. 1) THEN
                  ! First snow element also needs to define the first nodal temperature
                  t_sn_n(i,top(i))   = min(t(i),t0_melt)
               ENDIF
               t_sn_n(i,top(i) + 1) = t_sn_n(i,top(i))
               t_sn(i,top(i))       = 0.5_wp * (t_sn_n(i,top(i)) + t_sn_n(i,top(i)+1))

               ! Reset new snow storage
               hn_sn(i) = 0.0_wp

               alpha_sn(i) = 0.85_wp ! FIXME: link the hard-coded value to some config store

            ELSE ! ... otherwise merge bottom two layers to make room for new layer

               ! -----------------------
               ! Merge bottom two layers
               ! -----------------------

               theta_i(i,1) =   (theta_i(i,2)*dzm_sn(i,2)                                        &
                  +  theta_i(i,1)*dzm_sn(i,1)) / (dzm_sn(i,2) + dzm_sn(i,1))    ! volumetric ice content

               theta_w(i,1) =   (theta_w(i,2)*dzm_sn(i,2)                                        &
                  +  theta_w(i,1)*dzm_sn(i,1)) / (dzm_sn(i,2) + dzm_sn(i,1))    ! volumetric water content

               t_sn(i,1)    =   (t_sn(i,2)*dzm_sn(i,2)                                           &
                  +  t_sn(i,1)*dzm_sn(i,1))    / (dzm_sn(i,2) + dzm_sn(i,1))    ! snow layer temperature


               dzm_sn(i,1)  = dzm_sn(i,2) + dzm_sn(i,1)                                     ! layer thickness


               ! Update top layer index
               top(i) = top(i) - 1


               ! -----------------------------
               ! Move layers down to free top layer
               ! -----------------------------

               ! Layer thickness
               tmp_sn(1:ke_snow) = dzm_sn(i,1:ke_snow)

               DO ksn = ke_snow, 3, -1
                  dzm_sn(i,ksn-1)  = tmp_sn(ksn)
               ENDDO

               ! Volumetric ice content
               tmp_sn(1:ke_snow) = theta_i(i,1:ke_snow)

               DO ksn = ke_snow, 3, -1
                  theta_i(i,ksn-1)  = tmp_sn(ksn)
               ENDDO

               ! Volumetric water content
               tmp_sn(1:ke_snow) = theta_w(i,1:ke_snow)

               DO ksn = ke_snow, 3, -1
                  theta_w(i,ksn-1)  = tmp_sn(ksn)
               ENDDO

               ! Layer temperature
               tmp_sn(1:ke_snow) = t_sn(i,1:ke_snow)

               DO ksn = ke_snow, 3, -1
                  t_sn(i,ksn-1)  = tmp_sn(ksn)
               ENDDO

               ! nodal temperature
               tmp_sn_n(1:ke_snow+1) = t_sn_n(i,1:ke_snow+1)
               do ksn = ke_snow+1,3,-1
                  t_sn_n(i,ksn-1) = tmp_sn_n(ksn)
               enddo

               ! -------------------
               ! Add new snow layer and assign properties
               ! -------------------

               ! Limit top layer index it can/should only be ke_snow here
               top(i) = MIN(top(i) + 1 , ke_snow)

               ! Assign properties
               dzm_sn(i,top(i))   = hn_sn(i)/rho_hn(i)

               rho_sn(i,top(i))   = rho_hn(i)

               theta_i(i,top(i))  = rho_hn(i) / rho_i
               theta_w(i,top(i))  = 0.0_wp           ! new snow is always dry
               theta_a(i,top(i))  = 1.0_wp - theta_i(i,top(i)) - theta_w(i,top(i))

               t_sn_n(i,top(i) + 1) = t_sn_n(i,top(i))
               t_sn(i,top(i))       = 0.5_wp * (t_sn_n(i,top(i)) + t_sn_n(i,top(i)+1))

               ! Reset new snow storage
               hn_sn(i) = 0.0_wp

               alpha_sn(i) = 0.85_wp ! FIXME: link the hard-coded value to some config store

            ENDIF

         ELSE ! ... not enough snow and check if bottom layer can be splitted

            IF(top(i) .LT. ke_snow .AND. dzm_sn(i,1) .GT. (2.0_wp*max_height_layer)) THEN ! bottom layer can be split

               ! -----------------------
               ! Move layers up
               ! -----------------------

               ! Layer thickness
               tmp_sn(1:ke_snow) = dzm_sn(i,1:ke_snow)

               DO ksn = 2, ke_snow-1, 1
                  dzm_sn(i,ksn+1) = tmp_sn(ksn)
               ENDDO

               ! Volumetric ice content
               tmp_sn(1:ke_snow) = theta_i(i,1:ke_snow)

               DO ksn = 2, ke_snow-1, 1
                  theta_i(i,ksn+1) = tmp_sn(ksn)
               ENDDO

               ! Volumetric water content
               tmp_sn(1:ke_snow) = theta_w(i,1:ke_snow)


               DO ksn = 2, ke_snow-1, 1
                  theta_w(i,ksn+1) = tmp_sn(ksn)
               ENDDO

               ! Layer temperature
               tmp_sn(1:ke_snow) = t_sn(i,1:ke_snow)
               DO ksn = 2, ke_snow-1, 1
                  t_sn(i,ksn+1) = tmp_sn(ksn)
               ENDDO

               ! nodal temperatures
               tmp_sn_n(1:ke_snow+1) = t_sn_n(i,1:ke_snow+1)
               do ksn = 2,ke_snow,1
                  t_sn_n(i,ksn+1) = tmp_sn_n(ksn)
               enddo

               ! ---------------------
               ! Split bottom layer
               ! --------------------

               dzm_sn(i,1)  = dzm_sn(i,1) - max_height_layer
               dzm_sn(i,2)  = max_height_layer

               theta_i(i,2) = theta_i(i,1)
               theta_w(i,2) = theta_w(i,1)

               t_sn(i,2)    = t_sn(i,1)

               ! adjusting nodal temperatures
               t_sn_n(i,2) = t_sn_n(i,1) + ((t_sn_n(i,3) - t_sn_n(i,1)) / (dzm_sn(i,1) + dzm_sn(i,2))) * dzm_sn(i,1)

               ! Update to layer index
               top(i) = top(i) + 1


            ENDIF ! bottom layer can be split

         ENDIF  ! end not enough snow for a new layer


      END DO ! end of i



! =============================================================================
! - END subroutine: nix_stratigraphy
! ============================================================================

   END SUBROUTINE nix_stratigraphy



! =============================================================================
! + Begin subroutine: split_bottom_layer
! ============================================================================

   SUBROUTINE split_bottom_layer(nvec, ivstart, ivend         , &
   &                   top, ke_snow, dzm_sn         , &
   &                   theta_i, theta_w, theta_a    , &
   &                   t_sn, t_sn_n )


      ! Subroutine Arguments
      INTEGER, INTENT(IN)   :: &
         nvec          , & ! < array dimensions
         ivstart       , & ! < start index for computations in the parallel program
         ivend             ! < end index for computations in the parallel program

      INTEGER, DIMENSION(nvec), INTENT(INOUT) :: &
         top                 ! index of the first (top) layer index       (-)

      INTEGER, INTENT(IN)  ::  &
         ke_snow             ! number of snow layers

      REAL    (KIND = wp), DIMENSION(nvec,ke_snow),  INTENT(INOUT) :: &
         dzm_sn        , &   ! layer thickness between main levels        (m)
         theta_i       , &   ! volumetric ice content                     (-)
         theta_w       , &   ! water ice content                          (-)
         theta_a       , &   ! air ice content                            (-)
         t_sn                ! snow temperature (main level)              (K)

      REAL    (KIND = wp), DIMENSION(nvec,ke_snow+1),  INTENT(INOUT) :: t_sn_n

      ! Local variables
      INTEGER ::           &

         i            , &   ! loop index in x-direction
         ksn                ! loop index in z-direction (snow layers)

      REAL    (KIND = wp), DIMENSION(nvec,ke_snow) :: &
         tmp_sn             ! temporary utility array used for copying

      REAL    (KIND = wp), DIMENSION(ke_snow+1) :: &
         tmp_sn_n             ! temporary utility array used for copying

      ! ------------------------------------------------------------------------------
      ! Begin subroutine
      ! ------------------------------------------------------------------------------

      DO i = ivstart, ivend

         IF (top(i) .EQ. 1) THEN ! ... and maximum number of layers is not reached ...

            ! -----------------------
            ! Move layers up
            ! -----------------------

            ! Layer thickness
            tmp_sn(i,1:ke_snow) = dzm_sn(i,1:ke_snow)

            DO ksn = 2, ke_snow-1, 1
               dzm_sn(i,ksn+1) = tmp_sn(i,ksn)
            ENDDO

            ! Volumetric ice content
            tmp_sn(i,1:ke_snow) = theta_i(i,1:ke_snow)

            DO ksn = 2, ke_snow-1, 1
               theta_i(i,ksn+1) = tmp_sn(i,ksn)
            ENDDO

            ! Volumetric water content
            tmp_sn(i,1:ke_snow) = theta_w(i,1:ke_snow)

            DO ksn = 2, ke_snow-1, 1
               theta_w(i,ksn+1) = tmp_sn(i,ksn)
            ENDDO

            ! Layer temperature
            tmp_sn(i,1:ke_snow) = t_sn(i,1:ke_snow)

            DO ksn = 2, ke_snow-1, 1
               t_sn(i,ksn+1) = tmp_sn(i,ksn)
            ENDDO

            ! nodal temperatures
            tmp_sn_n(1:ke_snow+1) = t_sn_n(i,1:ke_snow+1)
            do ksn = 2,ke_snow,1
               t_sn_n(i,ksn+1) = tmp_sn_n(ksn)
            enddo

            ! -----------------------
            ! Split bottom layer
            ! -----------------------

            dzm_sn(i,1) = dzm_sn(i,1) * 0.5_wp
            dzm_sn(i,2) = dzm_sn(i,1)

            theta_i(i,2) = theta_i(i,1)
            theta_w(i,2) = theta_w(i,1)

            t_sn(i,2)    = t_sn(i,1)

            t_sn_n(i,2) = t_sn(i,2)


            ! Update dependent variables
            top(i) = top(i) + 1

         ENDIF

      ENDDO ! end of i


! =============================================================================
! - END subroutine: nix_stratigraphy
! ============================================================================

   END SUBROUTINE split_bottom_layer






!------------------------------------------------------------------------------
! End of module mo_nix_blanc
!------------------------------------------------------------------------------

END MODULE mo_nix_stratigraphy




