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
! Begin of module mo_nix_heat_equation
! ------------------------------------------------------------------------------

MODULE mo_nix_heat_equation

  USE mo_kind,                    ONLY: wp
  USE mo_physical_constants,      ONLY: stbo             ! Stefan Boltzman Konstante

  USE mo_nix_constants,           ONLY: eps_div, e_snow, specific_heat_water

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

PUBLIC :: heat_equation_implicit

CONTAINS

! =============================================================================
! + Begin subroutine: heat_equation_implicit
! ============================================================================

   SUBROUTINE heat_equation_implicit(nvec, ivstart, ivend                   , &
       &                             ke_snow, ke_soil, top                  , &
       &                             dzm_sn, hcon_sn, hcap_sn, hdif_sn      , &
       &                             t_sn, t_sn_n                           , &
       &                             swflx_sn_abs,lwflx_sn_dn               , &
       &                             lwflx_sn_up,lhflx_sn, shflx_sn         , &
       &                             zrain_rate                             , &
       &                             hcon_so, t_so, dt                      , &
       &                             t_sn_sfc, tch_sn, rho_sn, t, theta_w)


      ! Subroutine arguments
      INTEGER, INTENT(IN)                          :: &
        nvec               , & ! < array dimensions
        ivstart            , & ! < start index for computations in the parallel program
        ivend              , & ! < end index for computations in the parallel program
        ke_snow            , & ! < number of snow layers
        ke_soil                ! < number of soil layers

      INTEGER, DIMENSION(nvec), INTENT(IN) :: &
        top                    ! top layer index

      REAL (KIND=wp), DIMENSION(nvec,ke_snow), INTENT(IN) :: &
        dzm_sn              , & ! < snow layer depth
        hcon_sn             , & ! <            conductivity
        hcap_sn             , & ! <            capacity
        hdif_sn             , & ! <            difusivity
        swflx_sn_abs        , & ! <            absorbed short-wave radiation
        theta_w                 ! <            liquid water content

      REAL (KIND=wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
        t_sn                    ! <            temperature

      REAL (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
        lwflx_sn_dn         , &
        lwflx_sn_up         , &
        lhflx_sn, shflx_sn


      REAL (KIND=wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) ::  &    ! TODO/HACK: can intent be IN, as in ICON?
        rho_sn

      REAL (KIND=wp), DIMENSION(nvec,ke_snow+1), INTENT(INOUT) :: &
        t_sn_n

      REAL (KIND=wp), DIMENSION(nvec,ke_soil+1), INTENT(IN) :: &
        hcon_so                 ! < soil layer conductivity

      REAL (KIND=wp), DIMENSION(nvec,0:ke_soil+1), INTENT(IN) :: &
        t_so                    ! <            temperature

      REAL (KIND=wp), INTENT(IN)  ::  &
        dt                      ! time step

      REAL (KIND=wp), DIMENSION(nvec), INTENT(INOUT)  ::  &
        t_sn_sfc               ! < snow surface temperature

      REAL (KIND=wp), DIMENSION(nvec), INTENT(IN)  ::  &
        tch_sn

      REAL (KIND=wp), DIMENSION(nvec), INTENT(IN) :: &
         t

      REAL (KIND=wp), DIMENSION(nvec), INTENT(IN) :: &
         zrain_rate     

      ! Local variables
      INTEGER :: ksn, i

      REAL (KIND=wp), DIMENSION(ke_snow+1) :: &
         a_matrix  , &
         b_matrix  , &
         c_matrix  , &
         d_matrix  , &
         tmp_t

       REAL (KIND=wp) :: &
         emiss        , &
         t_emiss      , &
         coeff        , &
         delta_rad

       REAL (KIND=wp) :: &
         gamma_r       , &      
         gamma_soil

       REAL (KIND=wp), PARAMETER :: &
         eps_div = 1.0E-6_wp

       REAL (KIND=wp), DIMENSION(ke_snow+1) :: &
          tmp_t_sn_n   
         
       REAL (KIND=wp), dimension(ke_snow+1) :: &
          U,      &
          dU,     &
          ddU

       REAL (KIND=wp) :: tol = 0.0001_wp    ! Solver tolerance for convergence in heat equation solver

       INTEGER        :: maxiter = 200      ! Maximum number of iterations in heat equation solver
       INTEGER        :: iteration

       REAL (KIND=wp)  :: c, k, maxddU
       REAL (KIND=wp)  :: alpha_shf

    ! ------------------------------------------------------------------------------
    ! + Solve 1D heat equation - Implicit
    ! ------------------------------------------------------------------------------

    DO i = ivstart, ivend

      IF(top(i) .GT. 1) THEN  !!!snow on the ground

        U(:) = t_sn_n(i,:)
        dU   = 0.0_wp

        iterloop: DO iteration = 1,maxiter

        ! ---------------------------
        ! + Preparations
        ! ---------------------------

          a_matrix = 0.0_wp
          b_matrix = 0.0_wp
          c_matrix = 0.0_wp
          d_matrix = 0.0_wp
          ddU = dU
          dU = 0.0_wp

          ! Net longwave radiation coefficient: non-linear dependence on snow surface temperature
          emiss     = lwflx_sn_dn(i)/(stbo*t(i)*t(i)*t(i)*t(i))
          t_emiss   = sqrt(sqrt(emiss)) * t(i)
          delta_rad = stbo * (t_emiss + U(top(i)+1)) * (t_emiss * t_emiss + U(top(i)+1) * U(top(i)+1))

          ! Exchange coefficient for sensible heat flux
          IF ( abs( t(i) - t_sn_sfc(i) ) < eps_div  ) THEN
            alpha_shf = shflx_sn(i) / ( eps_div ) !1.0
          ELSE
            alpha_shf = shflx_sn(i) / ( t(i) - t_sn_sfc(i) )
          ENDIF

          gamma_r    = zrain_rate(i) * specific_heat_water;    ! exchange coefficient for rain energy
          gamma_soil = (hcon_so(i,1)/0.005_wp)                 ! FIXME: hardcoded length of the top soil layer

          ! Set up the solver for the tridiagonal form: a_matrix(i)*x(i-1) + b_matrix(i)*x(i) + c_matrix(i)*x(i+1) = d_matrix(i)

          DO ksn = 1,top(i),1
          
            c = (dzm_sn(i,ksn) * rho_sn(i,ksn) * hcap_sn(i,ksn)) / ( 6.0_wp * dt )
            k = hcon_sn(i,ksn) / dzm_sn(i,ksn)

            ! SNOWPACK: Se[0][0] = Se[1][1] = k;
            b_matrix(ksn)   = b_matrix(ksn)   + k
            b_matrix(ksn+1) = b_matrix(ksn+1) + k

            ! SNOWPACK: Se[0][1] = Se[1][0] = -k;
            a_matrix(ksn+1)   = a_matrix(ksn+1) - k
            c_matrix(ksn)     = c_matrix(ksn)   - k

            ! SNOWPACK: Add the implicit time integration term to the right hand side
            d_matrix(ksn)   = d_matrix(ksn)   - (k * t_sn_n(i,ksn)   - k * t_sn_n(i,ksn+1))
            d_matrix(ksn+1) = d_matrix(ksn+1) - (k * t_sn_n(i,ksn+1) - k * t_sn_n(i,ksn)  )

            ! SNOWPACK: Now add the heat capacitity matrix
            b_matrix(ksn)   = b_matrix(ksn)   + 2.0_wp * c
            b_matrix(ksn+1) = b_matrix(ksn+1) + 2.0_wp * c

            a_matrix(ksn+1) = a_matrix(ksn+1) + c
            c_matrix(ksn)   = c_matrix(ksn)   + c

            ! SNOWPACK: Heat the element via short-wave radiation
            d_matrix(ksn+1) = d_matrix(ksn+1) + swflx_sn_abs(i,ksn)

            ! Add upper boundary conditions (the surface energy balance)
            IF (ksn .EQ. top(i)) THEN
              IF(theta_w(i,ksn) .GT. 0.0_wp) then
                ! Explicit
                d_matrix(ksn+1) = d_matrix(ksn+1) +                                                                &
                  &                 + (lwflx_sn_dn(i)                                                              & ! Net longwave radiation
                  &                 - e_snow*stbo*t_sn_n(i,ksn+1)*t_sn_n(i,ksn+1)*t_sn_n(i,ksn+1)*t_sn_n(i,ksn+1)) &
                  &                 + lhflx_sn(i)                                                                  & ! Latent heat flux
                  &                 + shflx_sn(i)                                                                  & ! Sensible heat flux
                  &                 + gamma_r * (t(i) - t_sn_n(i,ksn+1))                                             ! Rain energy
               ELSE
                ! Implicit
                d_matrix(ksn+1) = d_matrix(ksn+1) + lhflx_sn(i)                                                      ! Latent heat flux
                d_matrix(ksn+1) = d_matrix(ksn+1) + alpha_shf * t(i)                                                 ! Sensible heat flux linearization
                d_matrix(ksn+1) = d_matrix(ksn+1) + delta_rad * t_emiss                                              ! Longwave radiation linearization
                d_matrix(ksn+1) = d_matrix(ksn+1) + gamma_r * t(i)                                                   ! Rain energy linearization
                b_matrix(ksn+1) = b_matrix(ksn+1) + alpha_shf + delta_rad + gamma_r
                d_matrix(ksn+1) = d_matrix(ksn+1) - (alpha_shf + delta_rad + gamma_r) * t_sn_n(i,ksn+1)
               ENDIF
             ENDIF

             ! Add lower boundary condition (Dirichlet BC)
             IF (ksn .EQ. 1) THEN
               b_matrix(ksn) = 1E12_wp
             ENDIF

           ENDDO !ksn

           ! ------------------------------------------------------------
           ! Solve the system - Thomas Algorithm
           ! ------------------------------------------------------------

           ! Step 1: forward elimination
           DO ksn=2,top(i)+1
             coeff  = a_matrix(ksn)/b_matrix(ksn-1)
             b_matrix(ksn) = b_matrix(ksn) - coeff * c_matrix(ksn-1)
             d_matrix(ksn) = d_matrix(ksn) - coeff * d_matrix(ksn-1)
           ENDDO

           dU = 0.0_wp
           ! Step 2: back substitution
           dU(top(i)+1) = d_matrix(top(i)+1)/b_matrix(top(i)+1)

           DO ksn=top(i),1,-1
             dU(ksn) = (d_matrix(ksn) - c_matrix(ksn) * dU(ksn+1))/b_matrix(ksn)
           ENDDO

           DO ksn = 1,top(i)+1
             ddU(ksn) = dU(ksn) - ddU(ksn)
             IF (maxddU .lt. abs(ddU(ksn)) .OR. ksn .EQ. 1) THEN
               maxddU = abs(ddU(ksn))
             ENDIF
                  U(ksn) = U(ksn) + ddU(ksn)
           ENDDO

           IF (maxddU < tol) THEN
             EXIT iterloop
           ENDIF

        ENDDO iterloop ! iteration

        ! Update nodal temperature
        DO ksn=1,top(i)+1
          t_sn_n(i,ksn) = U(ksn)
        ENDDO

        ! Update element temperature
        DO ksn=1,top(i)
          t_sn(i,ksn) = 0.5_wp * (t_sn_n(i,ksn) + t_sn_n(i,ksn+1))
        ENDDO

        t_sn_sfc(i) = t_sn_n(i,top(i)+1)
        ! shflx_sn(i) = alpha_shf * ( t(i) - t_sn_sfc(i) )


      ENDIF !snow on the ground

    ENDDO !end of i

! =============================================================================
! - END subroutine: heat_equation_implicit
! ============================================================================

  END SUBROUTINE heat_equation_implicit


!------------------------------------------------------------------------------
! End of module mo_nix_heat_equation
!------------------------------------------------------------------------------

END MODULE mo_nix_heat_equation




