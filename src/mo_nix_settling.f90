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
! Begin of module mo_nix_settling
! ------------------------------------------------------------------------------

MODULE mo_nix_settling

   USE mo_kind,                    ONLY: wp

   USE mo_nix_config,              ONLY: min_height_layer

   USE mo_physical_constants,      ONLY: t0_melt => tmelt   , & ! absolute zero for temperature
      rho_w   => rhoh2o  , & ! density of liquid water (kg/m^3)
      g       => grav        ! [m/s2] av. gravitational acceleration

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: aggregate_layers
   PUBLIC :: settling

CONTAINS

! =============================================================================
! + Begin subroutine: aggregate_layers
! ============================================================================

   SUBROUTINE aggregate_layers(nvec, ivstart, ivend, ke_soil, ke_snow     , &
   &                      top, dzm_sn, hm_sn, zm_sn                  , &
   &                      theta_i, theta_i_old, theta_w, theta_a     , &
   &                      hcap_sn, hcon_sn, hdif_sn                  , &
   &                      t_sn, t_sn_n, rho_sn, mass_sn, dt)

      ! Subroutine arguments
      INTEGER, INTENT(IN)                          :: &
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow            , & ! < number of snow layers
         ke_soil

      INTEGER, DIMENSION(nvec), INTENT(INOUT) :: &
         top                    ! top layer index

      REAL (KIND=wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
         dzm_sn             , & ! < snow layer thickness
         zm_sn              , & ! <            depth main level
         hm_sn              , & ! <            height main level
         theta_i            , & ! <            volumetric ice content
         theta_i_old        , & ! <            volumetric ice content - old values before phase change
         theta_w            , & ! <            volumetric water content
         theta_a            , & ! <            volumetric air content
         hcap_sn            , & ! <            heat capacity
         hcon_sn            , & ! <            heat conductivity
         hdif_sn            , & ! <            heat diffusivity
         t_sn               , & ! <            temperature
         rho_sn             , & ! <            density
         mass_sn                ! <            mass

      REAL (KIND=wp), DIMENSION(nvec,ke_snow+1), INTENT(INOUT) :: t_sn_n

      REAL (KIND=wp), INTENT(IN)  :: &
         dt


      ! Local variables
      INTEGER ::           &
         i,j,k                  , &
         ksn


      ! ------------------------------------------------------------------------------
      ! + Aggregate Layers
      ! ------------------------------------------------------------------------------

      DO i = ivstart, ivend

         IF(top(i) .GE. 2) THEN ! snow on the ground

            ! Aggregate profile if feasible - layers might have melted during phase change
            DO ksn = top(i), 2, -1


               IF(dzm_sn(i,ksn) .LT. min_height_layer .OR. theta_i(i,ksn) .LT. 0.01_wp) THEN ! layer is quite thin

                  ! Aggregate with lower layer - Adjust properties ...
                  theta_i(i,ksn-1)     = (theta_i(i,ksn)*dzm_sn(i,ksn)                                      &
                     + theta_i(i,ksn-1)*dzm_sn(i,ksn-1)) / (dzm_sn(i,ksn) + dzm_sn(i,ksn-1))      ! volumetric ice content

                  theta_i_old(i,ksn-1) = (theta_i_old(i,ksn)*dzm_sn(i,ksn)                                      &
                     + theta_i_old(i,ksn-1)*dzm_sn(i,ksn-1)) / (dzm_sn(i,ksn) + dzm_sn(i,ksn-1))  ! volumetric ice content

                  theta_w(i,ksn-1)     = (theta_w(i,ksn)*dzm_sn(i,ksn)                                      &
                     + theta_w(i,ksn-1)*dzm_sn(i,ksn-1)) / (dzm_sn(i,ksn) + dzm_sn(i,ksn-1))      ! volumetric water content

                  t_sn(i,ksn-1)        = ( t_sn(i,ksn)*dzm_sn(i,ksn) * rho_sn(i,ksn)                         &
                     + t_sn(i,ksn-1)*dzm_sn(i,ksn-1) * rho_sn(i,ksn-1) ) /  &
                  & (dzm_sn(i,ksn)*rho_sn(i,ksn) + dzm_sn(i,ksn-1) * rho_sn(i,ksn-1) )      ! volumetric ice content

                  ! CHECK THE ABOVE REFORMULATION for t_sn !!!!!

                  dzm_sn(i,ksn-1)      = dzm_sn(i,ksn) + dzm_sn(i,ksn-1)                          ! layer thickness


                  t_sn_n(i,ksn) = 2.0_wp * t_sn(i,ksn-1) - t_sn_n(i,ksn-1) ! adjusting nodal temperature for the upper node of the merged cell

                  ! ... and reset properties - FIXME: Doing this here means we have to give additional fields we don't need here to the
                  !                                   subroutine. Doing this reset in update_nix_state() might make more sense but would
                  !                                   require some adaptations like a new field a marker that identifies
                  !                                   melted/aggregated layers.
                  hm_sn(i,ksn)       = 0.0_wp
                  zm_sn(i,ksn)       = 0.0_wp
                  dzm_sn(i,ksn)      = 0.0_wp

                  theta_i(i,ksn)     = 0.0_wp
                  theta_w(i,ksn)     = 0.0_wp
                  theta_a(i,ksn)     = 0.0_wp

                  rho_sn(i,ksn)      = 0.0_wp
                  t_sn(i,ksn)        = 0.0_wp
                  mass_sn(i,ksn)     = 0.0_wp

                  hcap_sn(i,ksn)     = 0.0_wp
                  hcon_sn(i,ksn)     = 0.0_wp
                  hdif_sn(i,ksn)     = 0.0_wp

                  theta_i_old(i,ksn) = 0.0_wp

               ENDIF

            ENDDO ! end of ksn

            ! Move layers accordingly and find new top index
            j = 0
            DO k = 1,top(i),1

               IF( (dzm_sn(i,k) .GT. min_height_layer) .AND. (theta_i(i,k) .GT. 0.01_wp) ) THEN
                  !IF( theta_i(i,lay) .GT. 0.01_wp ) then
                  ! FIXME: Technially we want the second line as if-statement.
                  !        However, code freezes on GPU for COSMO when we use this line!
                  !        Investigate if that still holds for ICON

                  j=j+1

                  dzm_sn(i,j)      = dzm_sn(i,k)
                  theta_i(i,j)     = theta_i(i,k)
                  theta_w(i,j)     = theta_w(i,k)
                  rho_sn(i,j)      = rho_sn(i,k)
                  t_sn(i,j)        = t_sn(i,k)

                  t_sn_n(i,j) = t_sn_n(i,k)
                  t_sn_n(i,j+1) = t_sn_n(i,k+1)

                  theta_i_old(i,j) = theta_i_old(i,k)

               ENDIF

            ENDDO
            top(i) = j

            ! Reset layer properties
            DO ksn = top(i)+1, ke_snow, 1

               hm_sn(i,ksn)       = 0.0_wp
               zm_sn(i,ksn)       = 0.0_wp
               dzm_sn(i,ksn)      = 0.0_wp

               theta_i(i,ksn)     = 0.0_wp
               theta_w(i,ksn)     = 0.0_wp
               theta_a(i,ksn)     = 0.0_wp
               theta_i_old(i,ksn) = 0.0_wp

               rho_sn(i,ksn)      = 0.0_wp
               t_sn(i,ksn)        = 0.0_wp
               t_sn_n(i,ksn+1)    = 0.0_wp
               mass_sn(i,ksn)     = 0.0_wp

               hcap_sn(i,ksn)     = 0.0_wp
               hcon_sn(i,ksn)     = 0.0_wp
               hdif_sn(i,ksn)     = 0.0_wp

            ENDDO


         ENDIF ! snow on the ground

      ENDDO ! end of i




! =============================================================================
! - END subroutine: aggregate_layers
! ============================================================================

   END SUBROUTINE aggregate_layers



! =============================================================================
! + Begin subroutine: settling
! ============================================================================


   SUBROUTINE settling(nvec, ivstart, ivend, ke_soil, ke_snow     , &
   &             top, dzm_sn, theta_i, theta_i_old, theta_w  , &
   &             t_sn, rho_sn, mass_sn, dt, zm_sn, h_snow)

      ! Subroutine arguments
      INTEGER, INTENT(IN)                          :: &
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow            , & ! < number of snow layers
         ke_soil

      INTEGER, DIMENSION(nvec), INTENT(INOUT) :: &
         top                    ! top layer index

      REAL (KIND=wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
         dzm_sn             , & ! < snow layer thickness
         theta_i            , & ! <            volumetric ice content
         theta_i_old        , & ! <            volumetric ice content - old values before phase change
         theta_w            , & ! <            volumetric water content
         t_sn               , & ! <            temperature
         rho_sn             , & ! <            density
         mass_sn            , & ! <            mass
         zm_sn                  ! <            depth - main level

      REAL (KIND=wp), INTENT(IN)  :: &
         dt

      REAL (KIND=wp), DIMENSION(nvec), INTENT(INOUT) :: &
         h_snow                 ! < snow height



      ! Local variables
      INTEGER ::           &
         i                  , &
         ksn

      REAL (KIND = wp) :: &
         rate_1             , & ! settling rate for ice loss (s**-1)
         rate_2             , & ! overburden stress (s**-1)
         overburden         , & ! overburden load
         tot_rate           , & ! total settling rates (s**-1)
         ddz                    ! change of layer thickness

      REAL (KIND = wp)                          :: &
         dz_old             , & ! old layer thickness before settling     (m)
         dz_new                 ! new layer thickness after settling      (m)

      ! After Vionnet (2012)
      REAL(KIND=wp), PARAMETER ::  &

         a_eta = 0.1_wp     , &   !  default  0.1_wp
         b_eta = 0.023_wp   , &   !           0.023_wp
         c_eta = 250.0_wp   , &   !           250.0_wp
         eta_0 = 7.62237E6_wp     !           7.62237E6_wp

      REAL (KIND = wp) :: &

         f1                  , &
         f2                  , &
         eta



      ! ------------------------------------------------------------------------------
      ! + Start settling routine
      ! ------------------------------------------------------------------------------

      ! ----------------------
      ! Calculate settling rate
      ! ---------------------

      DO i = ivstart, ivend

         IF(top(i) .GE. 1) THEN

            ! ----------------------
            ! Initalizations
            ! ----------------------

            ddz          = 0.0_wp   ! total change in layer thickness
            overburden   = 0.0_wp   ! overburden stress (weight of the overlying layers times g)

            DO ksn = top(i), 1, -1

               ! Reset values for settling rates
               rate_1    = 0.0_wp      ! due to ice loss
               rate_2    = 0.0_wp      ! due to overburden stress

               tot_rate  = 0.0_wp      ! total settliing rate

               dz_old    = 0.0_wp      ! old layer thickness (before settling) ...
               dz_new    = 0.0_wp      ! ... and after settling

               ! ----------------------
               ! Calculate settling rate
               ! ---------------------

               !  Settling due to loss of ice
               IF(theta_i(i,ksn) <  theta_i_old(i,ksn)) THEN ! melting ocured in that layers
                  rate_1 = -1.0_wp/dt * MAX(0.0_wp,  (theta_i_old(i,ksn) - theta_i(i,ksn)) / theta_i_old(i,ksn) )
               ENDIF

               ! Settling due to overburden stress - Vionnet (2012)
               f1   = 1.0_wp / (1.0_wp + 60.0_wp * ((theta_w(i,ksn)*rho_w*dzm_sn(i,ksn)) / (rho_w * dzm_sn(i,ksn))))
               f2   = 4.0_wp

               eta = f1 * f2 * eta_0 * (rho_sn(i,ksn)/c_eta) * exp(a_eta*(t0_melt - t_sn(i,ksn)) + b_eta*rho_sn(i,ksn))
               rate_2 = -1.0_wp * (overburden + (mass_sn(i,ksn)*g/2.0_wp)) / eta

               ! Increase overburden stress FIXME: How to deal with slope angle should be overburden = m*g*cos(alpha)
               overburden = overburden + (mass_sn(i,ksn) * g)

               ! ----------------------
               ! Calculate change
               ! ---------------------

               ! ... of all (ice loss, overburden, destructive) settling rates (1/s)
               tot_rate = (rate_1*dt) + (rate_2*dt)

               ! ... of layer thickness, i.e. sum over all layers (m)
               ddz      = ddz + MAX(-1.0_wp * dzm_sn(i,ksn), dzm_sn(i,ksn) * tot_rate)

               dz_old = dzm_sn(i,ksn)
               dz_new = dz_old + MAX(-1.0_wp * dzm_sn(i,ksn), dzm_sn(i,ksn) * tot_rate)

               ! ... volumetric contents
               theta_i(i,ksn) = MAX(0.0_wp, theta_i(i,ksn) * (dz_old / dz_new))    ! ice content
               theta_w(i,ksn) = MAX(0.0_wp, theta_w(i,ksn) * (dz_old / dz_new))    ! water content

               ! ... of layer thickness (m)
               dzm_sn(i,ksn) = dz_new

            ENDDO ! ksn

            ! -------------------
            ! + Re-calculate layer depth from layer thickness
            ! -------------------

            DO ksn = top(i), 1, -1
               IF(ksn .EQ. top(i)) THEN
                  zm_sn(i,ksn) = dzm_sn(i,ksn)
               ELSE
                  zm_sn(i,ksn) = zm_sn(i,ksn+1) + dzm_sn(i,ksn)
               ENDIF
            ENDDO

            ! Calculate change in snow depth
            h_snow(i) = zm_sn(i,1)

         ENDIF

      ENDDO ! end of i







! =============================================================================
! - END subroutine: settling
! ============================================================================

   END SUBROUTINE settling





!------------------------------------------------------------------------------
! End of module mo_nix_settling
!------------------------------------------------------------------------------

END MODULE mo_nix_settling






