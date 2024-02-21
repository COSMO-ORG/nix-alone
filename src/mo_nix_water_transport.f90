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
! Begin of module mo_nix_water_transport
! ------------------------------------------------------------------------------

MODULE mo_nix_water_transport

   USE mo_kind,                    ONLY: wp

   USE mo_physical_constants,      ONLY: t0_melt => tmelt   , & ! absolute zero for temperature
      lh_v    => alv     , & ! latent heat of vapourization
      lh_s    => als     , & ! latent heat of sublimation
      lh_f    => alf     , & ! latent heat of fusion
      rho_w   => rhoh2o      ! density of liquid water (kg/m^3)

   USE mo_nix_constants,           ONLY: rho_i

   USE mo_nix_constants,             ONLY: eps_div

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: surface_mass_flux
   PUBLIC :: water_transport



CONTAINS

! =============================================================================
! + Begin subroutine: surface_mass_flux
! ============================================================================

   SUBROUTINE surface_mass_flux(nvec, ivstart, ivend, ke_soil, ke_snow  , &
   &                       top, lhflx_sn, dt, dzm_sn, rho_sn       , &
   &                       theta_i, theta_w, t_sn_sfc)

      ! Subroutine arguments
      INTEGER, INTENT(IN)                          :: &
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow            , & ! < number of snow layers
         ke_soil

      INTEGER, DIMENSION(nvec), INTENT(IN) :: &
         top                    ! top layer index

      INTEGER      :: &

         i                  , & ! loop index in x-direction
         ksn                , & ! loop index in z-direction (snow layers)
         kso

      REAL (KIND=wp), DIMENSION(nvec), INTENT(INOUT) :: &
         lhflx_sn               ! latent heat flux

      REAL (KIND=wp), INTENT(IN) :: &
         dt                     ! time step

      REAL (KIND=wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
         dzm_sn             , & ! snow layer thickness
         rho_sn             , & !            density
         theta_i            , & !            volumetric ice content
         theta_w                !            volumetric water content

      REAL (KIND=wp), DIMENSION(nvec), INTENT(IN) :: &
         t_sn_sfc               !            snow surface temperature



      ! Local variables
      REAL (KIND = wp)        :: &

         lhflx_sn_store   , & ! variable to keep track of used lhflx_sn
         dL               , & ! change of layer thickness
         dM               , & ! change of layer mass
         M                , & ! initial mass and volmetric content (water of ice)
         hoar             , & ! hoar mass
         dzm_sn_old           ! save old value of layer thickness

      REAL (KIND = wp), PARAMETER  :: &
         eps_div2 = 1.0E-12   ! very small number FIXME: add it to constants and make it availabel via USE
      !                          same for phase_change


      ! ------------------------------------------------------------------------------
      ! + Water Transport
      ! ------------------------------------------------------------------------------
      DO i = ivstart, ivend

         lhflx_sn_store = lhflx_sn(i)

         IF(top(i) .GE. 1) THEN ! snow on the ground

            ! -----------------------------------------------------------------------
            ! + Compute Sublimation/Deposition/Condensation/Evaporation
            ! -----------------------------------------------------------------------

            ! Initiate some values
            dL   = 0.0_wp
            dM   = 0.0_wp
            M    = 0.0_wp
            hoar = 0.0_wp

            ! --------------------------
            ! Latent heat flux towards the surface - mass gain
            ! ---------------------------

            IF(lhflx_sn_store .GT. eps_div2) THEN ! add mass

               IF(t_sn_sfc(i) .LT. t0_melt) THEN  ! add ice

                  dM =lhflx_sn_store * (dt)/lh_s     ! calculate mass change
                  lhflx_sn_store = 0.0_wp            ! reset latent heat flux, i.e. energy was used
                  hoar = dM

                  ! Adjust layer properties accordingly, keep snow density constant

                  dzm_sn_old = dzm_sn(i,top(i))
                  dL = dM/rho_sn(i,top(i))

                  dzm_sn(i,top(i))  = dzm_sn(i,top(i)) + dL
                  theta_i(i,top(i)) = theta_i(i,top(i)) * (dzm_sn_old/dzm_sn(i,top(i)))
                  theta_i(i,top(i)) = theta_i(i,top(i)) + (dM/(rho_i*dzm_sn(i,top(i))))
                  theta_w(i,top(i)) = theta_w(i,top(i)) * (dzm_sn_old/dzm_sn(i,top(i)))

               ELSE ! add water

                  dM = lhflx_sn_store * (dt)/lh_v                                          ! calculate mass change
                  lhflx_sn_store = 0.0_wp                                                  ! reset latent heat, i.e. energy was used
                  theta_w(i,top(i)) = theta_w(i,top(i)) + dM/(rho_w*dzm_sn(i,top(i)))   ! update volumetric water content

               ENDIF

            ELSE

               ! --------------------------
               ! Latent heat flux away from the surface - mass loss
               ! --------------------------

               IF(lhflx_sn_store .LT. (-1.0_wp*eps_div)) THEN ! additional check in case lh ist super small, but sligtly positive

                  ksn_loop: DO ksn = top(i), 1, -1 ! loop through snow layers

                     IF(theta_w(i,ksn) .GT. eps_div) THEN ! there is water, i.e. evaporate first

                        ! Calculate mass change
                        dM = lhflx_sn_store * (dt)/lh_v
                        M = theta_w(i,ksn) * rho_w * dzm_sn(i,ksn)

                        ! Check that you only take the available amount of water
                        IF(-dM .GE. M) THEN

                           dM = -M
                           theta_w(i,ksn) = theta_w(i,ksn) + dM/(rho_w*dzm_sn(i,ksn))

                        ELSE

                           theta_w(i,ksn) = theta_w(i,ksn) + dM/(rho_w*dzm_sn(i,ksn))

                        ENDIF

                        lhflx_sn_store = lhflx_sn_store - dM*lh_v/(dt) ! update energy used

                     ELSEIF (theta_i(i,ksn) .GT. eps_div) THEN ! there is no water then sublimate ice matrix

                        dM = lhflx_sn_store * (dt)/lh_s
                        M = theta_i(i,ksn) * rho_i * dzm_sn(i,ksn)

                        IF(-dM .GT. M) THEN ! all ice can be sublimated

                           dM = -M
                           theta_i(i,ksn) = 0.0_wp
                           dzm_sn(i,ksn)  = 0.0_wp

                        ELSE

                           dzm_sn_old = dzm_sn(i,ksn)
                           dL = dM/rho_sn(i,ksn)

                           dzm_sn(i,ksn)  = dzm_sn(i,ksn) + dL
                           theta_i(i,ksn) = theta_i(i,ksn) * (dzm_sn_old/dzm_sn(i,ksn))
                           theta_i(i,ksn) = theta_i(i,ksn) + (dM/(rho_i*dzm_sn(i,ksn)))
                           theta_w(i,ksn) = theta_w(i,ksn) * (dzm_sn_old/dzm_sn(i,ksn))

                        ENDIF

                        lhflx_sn_store = lhflx_sn_store - dM*lh_s/(dt) ! update energy used

                     ENDIF

                     ! there is still energy left, which should technically be used
                     ! by the soil layer for now let's erase it

                     IF(abs(lhflx_sn_store) .LT. (1.0_wp*eps_div)) THEN
                        lhflx_sn_store = 0.0_wp
                        exit ksn_loop
                     ENDIF

                  ENDDO ksn_loop ! end of ksn

               ENDIF ! additional check for lhflx_sn

            ENDIF ! surface temperature

         ENDIF ! snow on the ground

      ENDDO ! end of i


! =============================================================================
! - END subroutine: surface_mass_flux
! ============================================================================

   END SUBROUTINE surface_mass_flux



! =============================================================================
! + Begin subroutine: water_transport
! ============================================================================

   SUBROUTINE water_transport(nvec, ivstart, ivend, ke_soil, ke_snow  , &
   &                     top, lhflx_sn, dt, dzm_sn, rho_sn       , &
   &                     theta_i, theta_w, theta_a, hcap_sn      , &
   &                     runoff_sn, t_sn)

      ! Subroutine arguments
      INTEGER, INTENT(IN)                          :: &
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow            , & ! < number of snow layers
         ke_soil

      INTEGER, DIMENSION(nvec), INTENT(IN) :: &
         top                    ! top layer index

      REAL (KIND=wp), DIMENSION(nvec), INTENT(INOUT) :: &
         lhflx_sn           , & ! latent heat flux
         runoff_sn

      REAL (KIND=wp), INTENT(IN) :: &
         dt                     ! time step

      REAL (KIND=wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
         dzm_sn             , & ! snow layer thickness
         rho_sn             , & !            density
         theta_i            , & !            volumetric ice content
         theta_w            , & !            volumetric water content
         theta_a            , & !            volumetric air content
         hcap_sn            , & !            heat capacity
         t_sn                   !            temperature

      ! Local variables
      INTEGER      :: &

         i                  , & ! loop index in x-direction
         ksn                    ! loop index in z-direction (snow layers)

      REAL    (KIND = wp) :: &

         frac_rho          , &    ! fraction of density water to ice
         limit_theta_i     , &    ! abc-formula
         w_up              , &    ! water content of upper layer
         w_low             , &    ! water content of lower layer
         dtheta_w_up       , &    ! available water in upper layer
         dtheta_w_low      , &    ! available water in lower layer
         dtheta_w_low_x    , &    ! backup variable for dtheta_w_low
         theta_w_bot       , &    ! volumetric water content of bottom snow layer
         excess_water      , &    ! excess water
         res_wat_cont      , &    ! potential residual water content
         dtheta_w_sub_wtr  , &    ! additional storage capacity due to refreezing
         dz_up             , &    ! thickness above the layer of interest                                (m)
         dz_low                   ! thickness below the layer of interest                                (m)

      REAL (KIND = wp), DIMENSION(ke_snow)  :: &

         w_res               ! effective residual water content

      ! ------------------------------------------------------------------------------
      ! + Water Transport
      ! ------------------------------------------------------------------------------

      DO i = ivstart, ivend

         ! --------------------------
         ! + Initializations
         ! ---------------------------

         w_res(:)         = 0.0_wp
         dtheta_w_sub_wtr = 0.0_wp
         excess_water     = 0.0_wp
         dtheta_w_up      = 0.0_wp
         dtheta_w_low     = 0.0_wp


         frac_rho      = rho_w/rho_i
         limit_theta_i =   1.0_wp - frac_rho * ((1.0_wp + 0.0165_wp * frac_rho)                   &
            - SQRT((1.0_wp + 0.0165_wp * frac_rho)*(1.0_wp + 0.0165_wp * frac_rho)  &
            - 4.0_wp * frac_rho * 0.0264_wp)) / (2.0_wp * frac_rho)

         DO ksn = 1, top(i), 1

            ! Determine the additional storage capacity (of water) due to refreezing
            dtheta_w_sub_wtr = (hcap_sn(i,ksn) * rho_sn(i,ksn)) * (1.0_wp / lh_f ) * &
            & (1.0_wp / rho_w) * MAX(0.0_wp, (t0_melt - t_sn(i,ksn)))

            ! Estimate resiudal water (RWC) content by volume; Coleou and Lesaffre, 1998, Ann. Glaciol., 26, 64-68
            IF(theta_i(i,ksn) .GT. 0.23_wp) THEN
               res_wat_cont = 0.0264_wp + 0.0099_wp * (1.0_wp - theta_i(i,ksn)) / theta_i(i,ksn)
            ELSE
               res_wat_cont = 0.08_wp - 0.1023_wp * (theta_i(i,ksn) - 0.03_wp)
            ENDIF

            ! Limit residual water content
            res_wat_cont = MIN(res_wat_cont, 0.08_wp)  ! NOTE: only needed in case of theta_i < 0.03
            ! Effective residual water content
            w_res(ksn) = MIN((1.0_wp - theta_i(i,ksn)) * rho_i/rho_w, res_wat_cont + dtheta_w_sub_wtr)
            w_res(ksn) = MAX(0.0_wp, w_res(ksn))

         ENDDO ! end of ksn


         ! --------------------------
         ! Now start moving the water and adjust properties accordingly
         ! --------------------------

         DO ksn = top(i), 2, -1

            ! Reset excess water to zero
            IF(theta_i(i,ksn) .LT. eps_div) THEN ! no more ice in this layer only residual water which needs to be moved
               excess_water = theta_w(i,ksn) ! add residual water to excess water and ...
               theta_w(i,ksn) = 0.0_wp       ! reset water content of said layer
            ELSE
               excess_water = 0.0_wp
            ENDIF

            ! water content of the upper layer
            w_up = theta_w(i,ksn)

            IF(ksn .EQ. top(i) .AND. w_up .GT. 0.0_wp .AND. w_up .LE. w_res(ksn)) THEN

               ! In that case you need to update the volumetric air content and the
               ! density of the top element as it may have caught some rain! Only top element should be considered,
               ! as when rain would have  infiltrated lower elements as well, w_up > w_res.
               theta_a(i,ksn) = MAX(0.0_wp, 1.0_wp - theta_w(i,ksn) - theta_i(i,ksn))
               rho_sn(i,ksn)  = (theta_i(i,ksn) * rho_i) + (theta_w(i,ksn) * rho_w)

               ! We want positive densities
               IF(rho_sn(i,ksn) .LT. 0.0_wp) THEN
                  rho_sn(i,ksn) = 0.0_wp
               ENDIF

            ENDIF

            IF(w_up .GT. w_res(ksn) .OR. excess_water .GT. 0.0_wp) THEN

               ! ...water is being transfered
               dz_up       = dzm_sn(i,ksn)
               dz_low      = dzm_sn(i,ksn-1)
               w_low       = theta_w(i,ksn-1)
               dtheta_w_up = MAX(0.0_wp, w_up - w_res(ksn))

               IF(dtheta_w_up .GT. 0.0_wp .OR. excess_water .GT. 0.0_wp) THEN

                  ! ... dtheta_w_low is determined by also taking excess_water into
                  ! account. Maybe excess_water can be stored in this layer.

                  dtheta_w_low = dtheta_w_up * (dzm_sn(i,ksn)/dzm_sn(i,ksn-1)) + (excess_water/dzm_sn(i,ksn-1))

                  ! now check whether there is enough air left - you might not be able
                  ! to move the water
                  ! or/and water may refreeze and expand specifically, you might want
                  ! to create a water table over ice

                  IF( (dtheta_w_low + w_low) .GT. (rho_i/rho_w * (1.0_wp - theta_i(i,ksn-1))) ) THEN

                     ! Deal with excess water ... Look how much you can leave in the lower layer (ksn-1).
                     ! If you have too much water even for the lower element (more
                     ! melt or rain per time step than can be kept in this element), water is transferred
                     ! to excess_water. Excess_water moves the water downward, trying to insert the
                     ! water in lower layer.
                     dtheta_w_low_x = dtheta_w_low    ! make backup
                     dtheta_w_low   = MAX(0.0_wp, (rho_i/rho_w * (1.0_wp - theta_i(i,ksn-1)) - w_low))
                     ! All the water that could not be stored in lower layer is considered excess_water.
                     excess_water = (dtheta_w_low_x - dtheta_w_low) * dzm_sn(i,ksn-1)

                  ELSE

                     excess_water = 0.0_wp

                  ENDIF

                  ! update volumetric contents, masses and density
                  theta_w(i,ksn)   = w_up  - dtheta_w_up
                  theta_w(i,ksn-1) = w_low + dtheta_w_low

               ENDIF ! end positive water movement

            ENDIF ! end if( W_upper > Wres )

         ENDDO ! ksn

         ! -------------------------------------
         ! Special treatment for lowermost snow layer, i.e runoff at bottom layer
         ! Note: If we want ponding on surface layer, i.e. glacier ice, sea ice, rock
         ! etc. we need to do it here, e.g. if itype = ...
         ! -------------------------------------

         theta_w_bot = theta_w(i,1)

         IF (theta_w_bot .GT. w_res(1)) THEN

            ! Adjust dependent values accordingly
            theta_w(i,1) = w_res(1)
            theta_a(i,1) = 1.0_wp - theta_w(i,1) - theta_i(i,1)
            ! Put all excess water of bottom layer  in runoff, i.e. move it out of the snow cover
            ! Note: if one comments this out you get ponding
            runoff_sn(i) = runoff_sn(i) + excess_water +  dzm_sn(i,1) * (theta_w_bot - w_res(1))

         ENDIF

      ENDDO ! end of i








! =============================================================================
! - END subroutine: water_transport
! ============================================================================

   END SUBROUTINE water_transport



!------------------------------------------------------------------------------
! End of module mo_nix_water_transport
!------------------------------------------------------------------------------

END MODULE mo_nix_water_transport




