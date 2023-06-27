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
! Begin of module mo_nix_substrate
! ------------------------------------------------------------------------------

MODULE mo_nix_substrate

   USE mo_kind,                    ONLY: wp

   USE mo_nix_constants,           ONLY: itype_heatcond

   USE sfc_terra_data  ! All variables from this data module are used by
   ! this module. These variables start with letter "c"

   USE mo_physical_constants,      ONLY: rho_w => rhoh2o  ! density of liquid water (kg/m^3)


! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: calculate_soil_properties

CONTAINS


! =============================================================================
! + Begin subroutine: calculate_soil_properties
! ============================================================================

   SUBROUTINE calculate_soil_properties(nvec, ivstart, ivend, ke_soil,        &
   &                                   soiltyp_subs, w_so, w_so_ice, t_so,   &
   &                                   plcov, rootdp, zmls,                  &
   &                                   hcap_so, hcon_so )

      ! Subroutine arguments
      INTEGER                          :: &
      !
         nvec              , & ! < array dimensions
         ivstart           , & ! < start index for computations in the parallel program
         ivend             , & ! < end index for computations in the parallel program
         ke_soil               ! < number of soil layers minus climatological layer

      INTEGER, DIMENSION(nvec), INTENT(IN) :: &
         soiltyp_subs         ! type of the soil (keys 0-9)

      REAL    (KIND = wp), DIMENSION(nvec,0:ke_soil+1), INTENT(IN) :: &
         t_so                 ! soil temperature (main level)                 (  K  )

      REAL    (KIND = wp), DIMENSION(nvec,ke_soil+1), INTENT(IN) :: &
         w_so            , & ! total water conent (ice + liquid water)       (m H20)
         w_so_ice

      REAL    (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
         plcov            , & ! fraction of plant cover                         --
         rootdp               ! depth of the roots

      REAL    (KIND = wp), DIMENSION(ke_soil+1), INTENT(IN) :: &
         zmls                 ! processing soil level structure

      REAL (KIND = wp), DIMENSION(nvec,ke_soil+1), INTENT(OUT) :: &
         hcap_so          , & ! heat capacity of soil layers             [    ]
         hcon_so              ! heat conductivity of soil laters



      ! Local variables
      INTEGER :: &
         i                , & ! < loop index in x-direction
         kso              , & ! < loop index in z-direction (soil layers)
         mstyp                ! < soil type index


      REAL(KIND=wp)  ::  &

         zwqg             , & ! < mean of fcap and pwp
         z4wdpv           , & ! 4*zwqg/porv
         zzz                  ! utility variable

      REAL(KIND=wp) ::  &
         zfcap      (nvec,ke_soil+1)   , & ! field capacity of soil
         zpwp       (nvec,ke_soil+1)   , & ! plant wilting point  (fraction of volume)
         zporv      (nvec,ke_soil+1)   , & ! pore volume (fraction of volume)
         zalamtmp   (nvec,ke_soil)     , & ! heat conductivity
         zdlam      (nvec)             , & ! heat conductivity parameter
         zalam      (nvec,ke_soil)     , & ! heat conductivity
         hzalam     (nvec,ke_soil+1)   , & ! heat conductivity
         zw_fr      (nvec,ke_soil+1)   , & ! fractional total water content of soil layers
         zrocg      (nvec,ke_soil+1)   , & ! total volumetric heat capacity of soil
         zrocg_soil (nvec,ke_soil+1)   , & ! volumetric heat capacity of bare soil
         ziw_fr     (nvec,ke_soil+1)   , & ! fractional ice content of soil layer
         zlw_fr     (nvec,ke_soil+1)   , & ! fractional liqu. water content of soil layer
         zroc       (nvec,ke_soil+1)       ! heat capacity

      ! ------------------------------------------------------------------------------
      ! Calculate heat conductivity
      ! ------------------------------------------------------------------------------

      IF(itype_heatcond == 1) THEN

         DO i = ivstart, ivend

            mstyp     = soiltyp_subs(i)                  ! soil type

            zfcap     (i,:) = cfcap (mstyp)              ! field capacity
            zpwp      (i,:) = cpwp  (mstyp)              ! plant wilting point
            zporv     (i,:) = cporv (mstyp)              ! pore volume
            zdlam     (i)   = cala1 (mstyp)-cala0(mstyp) ! heat conductivity parameter

         ENDDO

         DO kso = 1, ke_soil
            DO i = ivstart, ivend
               zwqg         = 0.5_wp*(zfcap(i,kso) + zpwp(i,kso))
               z4wdpv       = 4._wp*zwqg/zporv(i,kso)

               zalamtmp(i,kso) =              zdlam(i)                         &
                  * (0.25_wp + 0.30_wp*zdlam(i)           &
                  / (1._wp+0.75_wp*zdlam(i)))             &
                  * MIN (z4wdpv, 1.0_wp + (z4wdpv-1.0_wp)   &
                  *(1.0_wp+0.35_wp*zdlam(i))              &
                  /(1.0_wp+1.95_wp*zdlam(i)))
            ENDDO
         ENDDO


         DO kso = 1, ke_soil ! FIXME: This excessive copy between arrays is not neccessary!
            ! However keep it for now to be 'identical' in terms of naming convention to TERRA
            DO i = ivstart, ivend
               zalam(i,kso)   = zalam(i,kso) + zalamtmp(i,kso)
               hzalam(i,kso)  = zalam(i,kso)
               hcon_so(i,kso) = hzalam(i,kso) ! final field for output
            ENDDO
         ENDDO

      ENDIF

      ! ------------------------------------------------------------------------------
      ! Calculate heat capacity
      ! ------------------------------------------------------------------------------

      DO i = ivstart, ivend
         zw_fr  (i,ke_soil+1)  = w_so(i,ke_soil+1)/zdzhs(ke_soil+1)
      ENDDO

      ! REORDER
      DO kso   = 1, ke_soil
         DO i = ivstart, ivend
            zw_fr   (i,kso)     = w_so(i,kso)/zdzhs(kso)
         ENDDO
      ENDDO

      DO kso = 1, ke_soil
         DO i = ivstart, ivend

            ! Scale soil heat capacity with organic fraction -> Chadburn et al., 2015
            IF (zmls(kso) < rootdp(i)) THEN
               zzz = plcov(i)*(rootdp(i)-zmls(kso))/rootdp(i)
               zrocg(i,kso)=(1.0_wp-zzz)*zrocg_soil(i,kso)+zzz*0.58E+06_wp
            END IF

         ENDDO
      ENDDO


      DO i = ivstart, ivend !FIXME: See comment above
         DO   kso = 1,ke_soil+1

            ziw_fr(i,kso) = w_so_ice(i,kso)/zdzhs(kso)                             ! ice frac.

            zlw_fr(i,kso) = zw_fr(i,kso) - ziw_fr(i,kso)                           ! liquid water frac.

            zroc(i,kso)   = zrocg(i,kso) + rho_w*zlw_fr(i,kso)*chc_w +          &  ! soil  heat capacity
               rho_w*ziw_fr(i,kso)*chc_i

            hcap_so(i,kso)   = zroc(i,kso)

         END DO
      END DO      !soil layers




! =============================================================================
! - END subroutine: calculate_soil_properties
! ============================================================================

   END SUBROUTINE calculate_soil_properties



!------------------------------------------------------------------------------
! End of module mo_nix_substrate.f90
!------------------------------------------------------------------------------

END MODULE mo_nix_substrate




