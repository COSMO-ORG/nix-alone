! Source: icon-nwp/src/atm_phy_schemes/mo_lookup_tables_constants.f90
! Constants used for the computation of lookup tables of the saturation
! mixing ratio over liquid water (*c_les*) or ice(*c_ies*)
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------
MODULE mo_lookup_tables_constants

    USE, INTRINSIC :: iso_fortran_env, ONLY: wp => real64
    USE mo_physical_constants, ONLY: alv, als, cpd, rd, rv, tmelt
  
  
    IMPLICIT NONE
  
    PRIVATE
  
    PUBLIC :: c1es, c2es, c3les, c3ies, c4les, c4ies, c5les, c5ies, &
      &       c5alvcp, c5alscp, alvdcp, alsdcp
  
    REAL (wp), PARAMETER :: c1es  = 610.78_wp              !
    REAL (wp), PARAMETER :: c2es  = c1es*rd/rv             !
    REAL (wp), PARAMETER :: c3les = 17.269_wp              !
    REAL (wp), PARAMETER :: c3ies = 21.875_wp              !
    REAL (wp), PARAMETER :: c4les = 35.86_wp               !
    REAL (wp), PARAMETER :: c4ies = 7.66_wp                !
    REAL (wp), PARAMETER :: c5les = c3les*(tmelt-c4les)    !
    REAL (wp), PARAMETER :: c5ies = c3ies*(tmelt-c4ies)    !
    REAL (wp), PARAMETER :: c5alvcp = c5les*alv/cpd        !
    REAL (wp), PARAMETER :: c5alscp = c5ies*als/cpd        !
    REAL (wp), PARAMETER :: alvdcp  = alv/cpd              !
    REAL (wp), PARAMETER :: alsdcp  = als/cpd              !
    !$ACC DECLARE COPYIN(c1es, c2es, c3les, c3ies, c4les, c4ies, c5les, c5ies) &
    !$ACC   COPYIN(c5alvcp, c5alscp, alvdcp, alsdcp)
  
END MODULE mo_lookup_tables_constants




!  Source: icon-nwp/src/atm_phy_schemes/mo_satad.f90
!
!  Description:
!  This module provides service utilities for meteorological calculations.
!
! Routines (module procedure)
!
!     - pres_sat_ice
!       Saturation water vapour pressure
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------
MODULE mo_icon_functions

    USE, INTRINSIC :: iso_fortran_env, ONLY: ireals => real64, iintegers =>  int32

    USE mo_physical_constants, ONLY: b3    => tmelt !!

    USE mo_lookup_tables_constants, ONLY:  &
    b1    => c1es  , & !! constants for computing the sat. vapour
    b2w   => c3les , & !! pressure over water (l) and ice (i)
    b2i   => c3ies , & !!               -- " --
    b4w   => c4les , & !!               -- " --
    b4i   => c4ies , & !!               -- " --
    b234w => c5les , & !!               -- " --
    b234i => c5ies     !!               -- " --

    PUBLIC  :: sat_pres_ice

    INTEGER, PARAMETER :: ipsat = 1    ! (1) Tetens (1930)
                                       ! (2) Murphy-Koop for liq and ice 

    real(ireals), parameter  ::      &
        &  c1i_mk = 9.550426,    &  ! coefficients in Murphy and Koop saturation vapor pressure
        &  c2i_mk = 5723.265,    &  ! over ice and over  liquid water
        &  c3i_mk = 3.53068,     &
        &  c4i_mk = 0.00728332,  &
        &  c1w_mk = 54.842763,   &
        &  c2w_mk = 6763.22,     &
        &  c3w_mk = 4.210,       &
        &  c4w_mk = 0.000367,    &
        &  c5w_mk = 53.878,      &
        &  c6w_mk = 1331.22,     &
        &  c7w_mk = 9.44523,     &
        &  c8w_mk = 0.014025,    &
        &  xi_mk  = 0.0415,      &
        &  t0_mk  = 218.8

    CONTAINS

    ELEMENTAL FUNCTION sat_pres_ice(temp)
    IMPLICIT NONE
    REAL (KIND=ireals)              :: sat_pres_ice
    REAL (KIND=ireals), INTENT(IN)  :: temp

    !$ACC ROUTINE SEQ

    IF (ipsat <= 1) THEN
        sat_pres_ice = b1*EXP( b2i*(temp-b3)/(temp-b4i) )
    ELSEIF (ipsat==2 .OR. ipsat==3) THEN
        sat_pres_ice = psati_murphykoop(temp)
    ENDIF

    END FUNCTION sat_pres_ice

    ELEMENTAL FUNCTION psati_murphykoop(tk)
    IMPLICIT NONE
    REAL(KIND=ireals)             :: psati_murphykoop
    REAL(KIND=ireals), intent(IN) :: tk
    
    ! Eq. (7) of Murphy and Koop (2005)
    psati_murphykoop = exp(c1i_mk - c2i_mk/tk + c3i_mk*log(tk) - c4i_mk*tk )
        
    END FUNCTION psati_murphykoop

END MODULE mo_icon_functions