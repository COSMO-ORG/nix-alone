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
! Begin of module mo_nix_phase_change
! ------------------------------------------------------------------------------

MODULE mo_nix_phase_change

   USE mo_kind,                    ONLY: wp

   USE mo_physical_constants,      ONLY: t0_melt => tmelt   , & ! absolute zero for temperature
      lh_f    => alf     , & ! latent heat of fusion
      rho_w   => rhoh2o      ! density of liquid water (kg/m^3)

   USE mo_nix_constants,           ONLY: theta_r, theta_s, rho_i

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: phase_change


CONTAINS

! =============================================================================
! + Begin subroutine: phase_change
! ============================================================================

   SUBROUTINE phase_change(nvec, ivstart, ivend, ke_soil, ke_snow  , &
   &                  top, dzm_sn, theta_i, theta_w, theta_a  , &
   &                  t_sn, t_sn_n, hcap_sn, rho_sn, dt)

      ! Subroutine arguments
      INTEGER, INTENT(IN)                          :: &
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow            , & ! < number of snow layers
         ke_soil

      INTEGER, DIMENSION(nvec), INTENT(IN) :: &
         top                    ! top layer index

      REAL (KIND=wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
         dzm_sn             , & ! < snow layer depth
         theta_i            , & ! <            volumetric ice content
         theta_w            , & ! <            volumetric water content
         theta_a            , & ! <            volumetric air content
         t_sn               , & ! <            temperature
         hcap_sn            , & ! <            capacity
         rho_sn                 ! <            density

      REAL (KIND=wp), DIMENSION(nvec,ke_snow+1), INTENT(INOUT) :: t_sn_n

      REAL (KIND=wp), INTENT(IN)  :: &
         dt


      ! Local variables
      INTEGER ::           &
         i                  , &
         ksn

      LOGICAL, DIMENSION(nvec)  :: &

         melt_flag     , & ! flag indicating melting; IF THEN TRUE
         freeze_flag       !      indicating freezing

      REAL (KIND=wp), PARAMETER :: &
         eps_div2 = 1.0E-12_wp

      REAL    (KIND = wp) ::  &
         dT_sub_melt        , & ! difference between current temperature and melting temperature
         A_sub_melt         , & ! coefficient A_sub_melt (see notes below)
         dtheta_i           , & ! change in volumetric ice content
         dtheta_w           , & ! change in volumetric water content
         q_mf               , & !
         q_rest




      ! ------------------------------------------------------------------------------
      ! + Calculate phase change
      ! ------------------------------------------------------------------------------

      ! ---------------------------
      ! + Initializations
      ! ---------------------------

      DO i = ivstart, ivend

         IF(top(i) .GE. 1) THEN ! snow on the ground

            ! Initiate
            melt_flag(i)   = .FALSE.
            freeze_flag(i) = .FALSE.

            DO ksn = top(i), 1, -1
               ! Metling
               IF(t_sn(i,ksn) .GT. t0_melt .AND. theta_i(i,ksn) .GT. eps_div2) THEN
                  melt_flag(i) = .TRUE.
               ENDIF
               ! Freezing
               IF(t_sn(i,ksn) .LT. t0_melt .AND. theta_w(i,ksn) .GT. (theta_r + eps_div2)) THEN
                  freeze_flag(i) = .TRUE.
               ENDIF
            ENDDO  !end of snow layers

         ENDIF

      ENDDO ! end of i

      ! ---------------------------
      ! + Melting
      ! ---------------------------

      DO i = ivstart, ivend

         IF(top(i) .GE. 1) THEN ! snow on the ground

            IF(melt_flag(i)) THEN ! melting

               q_rest = 0.0_wp

               DO ksn = top(i), 1, -1

                  ! --------------------
                  ! Now see if any melting is going on -- this implies that
                  ! (1) the temperature of the element is above/equal the melting temperature
                  ! (2) there is something to melt and
                  ! (3) there is enough room to place the meltwater
                  ! --------------------

                  IF(t_sn(i,ksn) .GE. t0_melt .AND. theta_i(i,ksn) .GT. 0.0_wp .AND. theta_w(i,ksn) .LT. theta_s) THEN

                     ! difference dT between actual melting temperature and layer temperature
                     dT_sub_melt = 0.0_wp
                     dT_sub_melt = t0_melt - t_sn(i,ksn)

                     ! --------------------
                     ! Now we take into account that there might be some extra energy that could
                     ! not be used by the element above because of complete melting
                     ! --------------------

                     dT_sub_melt = dT_sub_melt - (q_rest / (hcap_sn(i,ksn) * rho_sn(i,ksn) * dzm_sn(i,ksn)))

                     ! Only do it when there is real potential to melt
                     IF(dT_sub_melt .LT. 0.0_wp) THEN

                        ! --------------------
                        ! Determine the DECREASE in ice content and the INCREASE of water
                        ! content. Adapt A_sub_melt to compute mass changes
                        ! --------------------

                        A_sub_melt = (hcap_sn(i,ksn) * rho_sn(i,ksn)) / (rho_i * lh_f)
                        dtheta_i   = A_sub_melt * dT_sub_melt
                        dtheta_w   = - (rho_i / rho_w) * dtheta_i

                        ! --------------------
                        ! It could happen that there is enough energy available to melt more ice
                        ! than is present. You can only melt so much ice as is there ....
                        ! --------------------
                        IF( (theta_i(i,ksn) + dtheta_i) .LT. 0.0_wp) THEN
                           dtheta_i = - theta_i(i,ksn)
                           dtheta_w = - (rho_i / rho_w) * dtheta_i
                           dT_sub_melt = dtheta_i / A_sub_melt
                        ENDIF


                        ! --------------------
                        ! It could also be that you are trying to produce more water than is allowed.
                        ! --------------------
                        IF( (theta_w(i,ksn) + dtheta_w) .GT. theta_s) THEN
                           dtheta_w = theta_s - theta_w(i,ksn)
                           dtheta_i = - (rho_w/rho_i) *dtheta_w
                           dT_sub_melt = dtheta_i / A_sub_melt
                        ENDIF

                        ! --------------------
                        ! Reset/Recalculate properties
                        ! --------------------

                        ! Layer temperature and transfered energy
                        t_sn(i,ksn) = t_sn(i,ksn) + dT_sub_melt
                        IF(t_sn(i,ksn) .LE. t0_melt) THEN ! if melting occured it can only be at melt point
                           q_rest     = 0.0_wp
                           t_sn(i,ksn)  = t0_melt

                           t_sn_n(i,ksn) = t0_melt
                           t_sn_n(i,ksn+1) = t0_melt
                        ELSE
                           q_rest = hcap_sn(i,ksn) * rho_sn(i,ksn) * dzm_sn(i,ksn) * (t_sn(i,ksn) - t0_melt)
                           t_sn(i,ksn)   = t0_melt
                           t_sn_n(i,ksn) = t0_melt
                           t_sn_n(i,ksn+1) = t0_melt
                        ENDIF

                        ! Volumetric freezing power
                        q_mf     = q_mf + ((dtheta_i * rho_i * lh_f) / (dt) )
                        dtheta_w = dtheta_w

                        ! Contents of ice, water and air
                        theta_i(i,ksn) = theta_i(i,ksn) + dtheta_i
                        theta_w(i,ksn) = theta_w(i,ksn) + dtheta_w

                     ENDIF ! deltaT check

                  ENDIF ! melt check

               ENDDO ! ksn

               !FIXME: At this point in COSMO we do an update of hm_sn, h_snow, zm_sn, theta_a,
               !       rho_sn, hcap_sn, hcon_sn and mass_sn. Why Do we really need this?
               !       Suggest to do it in main after phase change has been called.

            ENDIF ! check for melt flag

         ENDIF ! snow on the ground

      ENDDO ! end of i

      ! ---------------------------
      ! + Freezing
      ! ---------------------------

      DO i = ivstart, ivend

         IF(top(i) .GE. 1) THEN ! snow on the ground

            IF(freeze_flag(i)) THEN ! freeze

               DO  ksn = top(i), 1, -1

                  ! -------------------------------------------
                  ! Freezing within the snowpack can occur if
                  ! (1) the temperature of the element is below freezing and
                  ! (2) if water is present to be refrozen
                  ! -------------------------------------------

                  IF(t_sn(i,ksn) .LT. t0_melt .AND. theta_w(i,ksn) .GT. theta_r) THEN

                     ! Difference dT between actual layer temperature and freezing temperature
                     dT_sub_melt = 0.0_wp
                     dT_sub_melt = t0_melt - t_sn(i,ksn)

                     ! Adapt A_sub_melt to compute mass change
                     A_sub_melt = (hcap_sn(i,ksn) * rho_sn(i,ksn)) / (rho_i * lh_f)

                     ! Compute change in volumetric contenst
                     dtheta_i = A_sub_melt * dT_sub_melt
                     dtheta_w = -(rho_i/rho_w) * dtheta_i

                     ! Make sure that there is enough water to refreeze
                     IF( (theta_w(i,ksn) + dtheta_w) .LT. theta_r) THEN
                        dtheta_w = -ABS(theta_w(i,ksn) - theta_r)
                        dtheta_i = -(rho_w / rho_i) * dtheta_w
                        dT_sub_melt       = dtheta_i / A_sub_melt
                     ENDIF

                     ! See if the layer is pure ice
                     IF( (theta_i(i,ksn) + theta_r + dtheta_i) .GE. 1.0_wp) THEN
                        dtheta_w = - ABS(theta_w(i,ksn) - theta_r)
                        dtheta_i = - (rho_w/rho_i) * dtheta_w
                        theta_i(i,ksn)  = 1.0_wp - theta_r
                        theta_w(i,ksn)  = theta_r
                        theta_a(i,ksn)  = 0.0_wp
                     ELSE
                        theta_i(i,ksn) = theta_i(i,ksn) + dtheta_i
                        theta_w(i,ksn) = theta_w(i,ksn) + dtheta_w
                        theta_a(i,ksn) = MAX(0.0_wp, 1.0_wp - theta_i(i,ksn) - theta_w(i,ksn))
                     ENDIF

                     ! --------------------
                     ! Reset/Recalculate properties
                     ! --------------------

                     ! Set some limits
                     IF(theta_w(i,ksn) .GE. 1.0_wp) THEN
                        theta_w(i,ksn) = 1.0_wp
                     ENDIF

                     !Compute the volumetric refreezing power
                     q_mf     = q_mf + ((dtheta_i * rho_i * lh_f) / (dt) )
                     dtheta_w = dtheta_w
                     t_sn(i,ksn)     = t_sn(i,ksn) + dT_sub_melt

                     t_sn_n(i,ksn)   = t_sn_n(i,ksn)   + dT_sub_melt
                     t_sn_n(i,ksn+1) = t_sn_n(i,ksn+1) + dT_sub_melt


                  ENDIF ! freeze check

               ENDDO ! ksn

               !FIXME: Same here. Why the update?

            ENDIF ! freeze

         ENDIF ! there is snow or not

      ENDDO ! end loop over horizontal pixels


! =============================================================================
! - END subroutine: phase_change
! ============================================================================

   END SUBROUTINE phase_change


!------------------------------------------------------------------------------
! End of module mo_nix_phase_change
!------------------------------------------------------------------------------

END MODULE mo_nix_phase_change

