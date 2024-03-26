module mo_nix_heat_equation

   use mo_kind,                    only: wp
   use mo_physical_constants,      only: stbo
   use mo_nix_constants,           only: eps_div, e_snow
! ------------------------------------------------------------------------------
! declarations
! ------------------------------------------------------------------------------

   implicit none

   private

!------------------------------------------------------------------------------
! anything public?
!------------------------------------------------------------------------------

   public :: heat_equation_wrapper

contains

   subroutine heat_equation_semi_implicit(nvec, ivstart, ivend, ke_soil, ke_snow  , &
   &                                 top, zm_sn, hcon_sn, hcap_sn, hdif_sn   , &
   &                                 t_sn, swflx_sn_abs, hcon_so, hcap_so    , &
   &                                 t_so, dt, for_sn, t_sn_sfc)


      ! subroutine arguments
      integer, intent(in)                          :: &
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow            , & ! < number of snow layers
         ke_soil

      integer, dimension(nvec), intent(in) :: &
         top                    ! top layer index

      ! vs: intent for t_sn !?
      real (kind=wp), dimension(nvec,ke_snow), intent(inout) :: &
         zm_sn              , & ! < snow layer depth
         hcon_sn            , & ! <            conductivity
         hcap_sn            , & ! <            capacity
         hdif_sn            , & ! <            difusivity
         t_sn               , & ! <            temperature
         swflx_sn_abs           ! <            absorbed short-wave radiation

      real (kind=wp), dimension(nvec,ke_soil+1), intent(in) :: &
         hcon_so            , & ! < soil layer conductivity
         hcap_so            , & ! <            capacity
         t_so                   ! <            temperature

      real (kind = wp), intent(in)  ::  &
         dt                     ! time step

      ! vs: intent for t_sn_sfc should be in-out !?
      real (kind = wp), dimension(nvec), intent(inout)  ::  &
         for_sn             , & ! < total atmospheric forcing
         t_sn_sfc               ! < snow surface temperature

      ! local variables
      real (kind=wp), dimension(-ke_snow+1:1) :: &
         zm_col             , & ! depth of main levels

         hcon_col           , & ! heat conductivity
         hcap_col           , & ! heat capacity
         hdif_col           , & ! heat diffusivity

         rho_col            , &  ! density

         t_col              , &  ! temperature of the snow/substrate column

         swflx_sn_abs_col   , &  ! absorbed short-wave radiaton in each layer

         alpha_col          , &  ! utility variables for building and solving the tri-diagonal matrix
         gamma_col          , &  !

         a                  , &  !
         b                  , &  !
         c                  , &  !
         d                  , &  !
         e                       ! final snow layer temperature

      real (kind=wp)       ::  &

         dlw_u_sn           , &  ! derivative of upwelling longwave radiation for snow on the ground    (w m-2)
         dz_up              , &  ! thickness above the layer of interest                                (m)
         dz_low             , &  ! thickness below the layer of interest                                (m)
         beta_col                ! utility variable

      integer      :: &

         i                  , & ! loop index in x-direction
         ksn                , & ! loop index in z-direction (snow layers)
         kso                , & ! loop index in z-direction (soil layers)
         counter

      real  (kind=wp), parameter           ::  &

         cn = 0.5_wp            ! cranck-nicholson factor



      ! ------------------------------------------------------------------------------
      ! + solve 1d heat equation
      ! ------------------------------------------------------------------------------

      ! ---------------------------
      ! + initializations
      ! ---------------------------

      do i = ivstart, ivend

         if(top(i) .gt. 1) then  !!!snow on the ground

            zm_col            = 0.0_wp
            hcon_col          = 0.0_wp
            hcap_col          = 0.0_wp
            rho_col           = 0.0_wp
            t_col             = 0.0_wp
            swflx_sn_abs_col  = 0.0_wp

            counter = 1

            do ksn = -top(i)+1, 1, 1

               if(ksn .le. 0) then  ! snow layers

                  if(ksn .eq. -top(i)+1) then

                     zm_col(ksn)           = zm_sn(i,top(i))
                     hcon_col(ksn)         = hcon_sn(i,top(i))
                     hcap_col(ksn)         = hcap_sn(i,top(i))
                     t_col(ksn)            = t_sn(i,top(i))
                     swflx_sn_abs_col(ksn) = swflx_sn_abs(i,top(i))

                  else

                     zm_col(ksn)            = zm_sn(i,top(i)-counter)
                     hcon_col(ksn)          = hcon_sn(i,top(i)-counter)
                     hcap_col(ksn)          = hcap_sn(i,top(i)-counter)
                     t_col(ksn)             = t_sn(i,top(i)-counter)
                     swflx_sn_abs_col(ksn)  = swflx_sn_abs(i,top(i)-counter)

                     counter = counter + 1

                  endif

               else  ! soil layers !fixme: we are using or used capacity and conductivity of teh first soil layer for all!!!
                  !       loop index stops at 1 so first layer only.

                  zm_col(ksn)            = zm_sn(i,1) + 0.005_wp
                  hcon_col(ksn)          = hcon_so(i,1)                  ! hcon(0)
                  hcap_col(ksn)          = hcap_so(i,1)                   ! hcap(0)
                  t_col(ksn)             = t_so(i,ksn)
                  swflx_sn_abs_col(ksn)  = 0.0_wp

               endif

            enddo ! end of column (snow+soil)


            ! ---------------------------
            ! + pre-calculations
            ! ---------------------------

            ! derivative of emitted long wave radiation
            dlw_u_sn = 0.0_wp

            ! calculate factors (diffusion) for the linear equations for ...
            do ksn = -top(i)+1, 1, 1

               if(ksn .le. 0) then ! ... except
                  alpha_col(ksn)   = dt / hcap_col(ksn)
                  hdif_col(ksn) = hcon_col(ksn) * (t_col(ksn+1) - t_col(ksn)) /  (zm_col(ksn+1) - zm_col(ksn))
               else !(ksn == substrate) ! ... bottom layer
                  alpha_col(ksn) = dt / hcap_col(ksn)
                  hdif_col(ksn) = 0.0_wp
               endif

            end do

            ! ---------------------------
            ! + setup tridiagonal matrix for set of linear equations for each layer ...
            ! ---------------------------

            do ksn = -top(i)+1, 1

               if(ksn .eq. -top(i)+1) then ! ... top layer

                  dz_low = zm_col(ksn+1) - zm_col(ksn)

                  a(ksn)     = 0.0_wp
                  b(ksn)     = 1 + (1.0_wp - cn) * alpha_col(ksn) * hcon_col(ksn)/dz_low - alpha_col(ksn) * dlw_u_sn
                  c(ksn)     = -   (1.0_wp - cn) * alpha_col(ksn) * hcon_col(ksn)/dz_low

                  d(ksn)     = t_col(ksn) + alpha_col(ksn) * (for_sn(i) - dlw_u_sn*t_col(ksn) + cn*hdif_col(ksn))

               elseif (ksn .le. 0) then ! ... inner layers

                  dz_up  = zm_col(ksn)   - zm_col(ksn-1)
                  dz_low = zm_col(ksn+1) - zm_col(ksn)

                  a(ksn) = -         (1.0_wp - cn)   * alpha_col(ksn) *  hcon_col(ksn-1)/dz_up
                  b(ksn) = 1.0_wp +  (1.0_wp - cn)   * alpha_col(ksn) * (hcon_col(ksn)  /dz_low + hcon_col(ksn-1)/dz_up)
                  c(ksn) = -         (1.0_wp - cn)   * alpha_col(ksn) *  hcon_col(ksn)  /dz_low

                  d(ksn) = t_col(ksn) + cn*alpha_col(ksn) * (hdif_col(ksn) - hdif_col(ksn-1)) + alpha_col(ksn)*swflx_sn_abs_col(ksn)

               elseif (ksn .eq. 1) then ! bottom layers

                  dz_up = zm_col(ksn)   - zm_col(ksn-1)

                  a(ksn) = -         (1.0_wp - cn) * alpha_col(ksn) * hcon_col(ksn-1)/dz_up
                  b(ksn) = 1.0_wp +  (1.0_wp - cn) * alpha_col(ksn) * hcon_col(ksn-1)/dz_up
                  c(ksn) = 0.0_wp

                  d(ksn)     = t_col(ksn) - cn*alpha_col(ksn)*hdif_col(ksn-1) + alpha_col(ksn)*hdif_col(ksn)

               endif

            enddo

            ! ---------------------------
            ! + solve the system - thomas algorithm
            ! ---------------------------

            beta_col = b(-top(i)+1)

            ! forward substitution
            do ksn = -top(i)+1, 1, 1

               if(ksn .ge. -top(i)+1) then

                  if(ksn .eq. -top(i)+1) then
                     e(ksn) = d(ksn) / beta_col
                  else
                     gamma_col(ksn) = c(ksn-1) / beta_col
                     beta_col       = b(ksn) - a(ksn) * gamma_col(ksn)
                     e(ksn)     = (d(ksn) - a(ksn) * e(ksn-1)) / beta_col
                  endif
               endif

            enddo

            ! backward substitution
            do ksn = 0, -top(i)+1, -1

               if(ksn .ge. -top(i)+1) then
                  e(ksn) = e(ksn) - gamma_col(ksn+1) * e(ksn+1)
               endif

            enddo

            ! ---------------------------
            ! + do some updating required for the next sections
            ! ---------------------------

            ! snow surface temperature
            t_sn_sfc(i) = e(-top(i)+1)

            ! snow layer temperature
            counter = 1
            do ksn = top(i),1,-1

               if(ksn .eq. top(i)) then
                  t_sn(i,ksn)    = e(-top(i)+1)
               else
                  t_sn(i,ksn)    = e(-top(i)+1+counter)
                  counter = counter + 1
               endif

            enddo

         endif ! end snow on the ground

      enddo ! end of i


   end subroutine heat_equation_semi_implicit

   subroutine heat_equation_wrapper(nvec, ivstart, ivend, &
   & ke_snow, top,  &
   & dzm_sn, hcon_sn, hcap_sn, hdif_sn   , &
   & t_sn, t_sn_n, &
   & swflx_sn_abs,lwflx_sn_dn, &
   & lwflx_sn_up,lhflx_sn, &
   & shflx_sn,  &
   & hcon_so, &
   & t_so,dt,t_sn_sfc, tch_sn, rho_sn, t, theta_w)

      ! subroutine arguments
      integer, intent(in)                          :: &
         nvec               , & ! < array dimensions
         ivstart            , & ! < start index for computations in the parallel program
         ivend              , & ! < end index for computations in the parallel program
         ke_snow                ! < number of snow layers

      integer, dimension(nvec), intent(in) :: &
         top                    ! top layer index

      ! vs: intent for t_sn !?
      real (kind=wp), dimension(nvec,ke_snow), intent(inout) :: &
         dzm_sn             , & ! < snow layer depth
         hcon_sn            , & ! <            conductivity
         hcap_sn            , & ! <            capacity
         hdif_sn            , & ! <            diffusivity
         t_sn               , & ! <            temperature
         swflx_sn_abs       , & ! <            absorbed short-wave radiation
         theta_w                ! <            liquid water content

      real (kind=wp), dimension(nvec), intent(inout) :: &
         lwflx_sn_dn        , &
         lwflx_sn_up        , &
         lhflx_sn, shflx_sn


      real (kind=wp), dimension(nvec,ke_snow), intent(inout) :: rho_sn

      real (kind=wp), dimension(nvec,ke_snow+1), intent(inout) :: &
         t_sn_n

      real (kind=wp), dimension(ke_snow+1) :: &
         U,      &
         dU,     &
         ddU

      real (kind=wp), dimension(nvec), intent(in) :: &
         hcon_so            , & ! < soil layer conductivity
         t_so                   ! <            temperature

      real (kind = wp), intent(in)  ::  &
         dt                     ! time step

      ! vs: intent for t_sn_sfc should be in-out !?
      real (kind = wp), dimension(nvec), intent(inout)  ::  &
      &   t_sn_sfc              ! < snow surface temperature

      real (kind = wp), dimension(nvec), intent(in)  ::  &
      &   tch_sn

      real (kind=wp), dimension(nvec), intent(in) :: t

      !!! local variables
      real (kind=wp) :: tol = 0.0001_wp    ! Solver tolerance for convergence in heat equation solver
      integer :: maxiter = 200             ! Maximum number of iterations in heat equation solver
      integer :: ksn, i, iteration

      real(kind=wp), dimension(ke_snow+1) :: a_matrix,b_matrix,c_matrix,d_matrix
      real(kind=wp) :: emiss, t_emiss, coeff, delta_rad
      real(kind=wp) :: c, k, maxddU
      real(kind=wp) :: alpha_shf, gamma_soil


      do i = ivstart, ivend
         if(top(i) .gt. 1) then            ! snow on the ground
            U = t_sn_n(i,:)
            ddU = 0.0_wp

            iterloop: do iteration = 1,maxiter

               a_matrix = 0.0_wp
               b_matrix = 0.0_wp
               c_matrix = 0.0_wp
               d_matrix = 0.0_wp
               ddU = dU
               dU = 0.0_wp

               ! Net longwave radiation coefficient: non-linear dependence on snow surface temperature
               emiss   = lwflx_sn_dn(i)/(stbo*t(i)*t(i)*t(i)*t(i))
               t_emiss = sqrt(sqrt(emiss)) * t(i)
               delta_rad = stbo * (t_emiss + U(top(i)+1)) * (t_emiss * t_emiss + U(top(i)+1) * U(top(i)+1))

               ! Exchange coefficient for sensible heat flux
               if ( abs( t(i) - t_sn_sfc(i) ) < eps_div  ) then
                  alpha_shf = shflx_sn(i) / ( eps_div ) !1.0
               else
                  alpha_shf = shflx_sn(i) / ( t(i) - t_sn_sfc(i) )
               endif

               gamma_soil = 0.0_wp !(hcon_so(i)/0.005_wp) ! VS: length of the top soil layer

               ! Set up the solver for the tridiagonal form:
               ! a_matrix(i)*x(i-1) + b_matrix(i)*x(i) + c_matrix(i)*x(i+1) = d_matrix(i)
               do ksn = 1,top(i),1
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
                  if (ksn .EQ. top(i)) THEN
                     ! TODO add rain energy
                     if(theta_w(i,ksn) .gt. 0.0_wp) then
                        ! Explicit
                        d_matrix(ksn+1) = d_matrix(ksn+1) + &
                        &                 + (lwflx_sn_dn(i) & ! Net longwave radiation
                        &                 - e_snow*stbo*t_sn_n(i,ksn+1)*t_sn_n(i,ksn+1)*t_sn_n(i,ksn+1)*t_sn_n(i,ksn+1)) &
                        &                 + lhflx_sn(i)     & ! Latent heat flux
                        &                 + shflx_sn(i)       ! Sensible heat flux
                     else
                        ! Implicit
                        d_matrix(ksn+1) = d_matrix(ksn+1) + lhflx_sn(i)                               ! Latent heat flux
                        d_matrix(ksn+1) = d_matrix(ksn+1) + alpha_shf * t(i)                          ! Sensible heat flux linearization
                        d_matrix(ksn+1) = d_matrix(ksn+1) + delta_rad * t_emiss                       ! Longwave radiation linearization
                        b_matrix(ksn+1) = b_matrix(ksn+1) + alpha_shf + delta_rad
                        d_matrix(ksn+1) = d_matrix(ksn+1) - (alpha_shf + delta_rad) * t_sn_n(i,ksn+1)
                     endif
                  endif

                  ! Add lower boundary condition (Dirichlet BC)
                  if (ksn .EQ. 1) THEN
                     b_matrix(ksn) = 1E12_wp
                  endif

               enddo


               ! ------------------------------------------------------------
               ! Solve the system - Thomas Algorithm
               ! ------------------------------------------------------------

               ! step 1: forward elimination
               do ksn=2,top(i)+1
                  coeff  = a_matrix(ksn)/b_matrix(ksn-1)
                  b_matrix(ksn) = b_matrix(ksn) - coeff * c_matrix(ksn-1)
                  d_matrix(ksn) = d_matrix(ksn) - coeff * d_matrix(ksn-1)
               enddo

               dU = 0.0_wp
               ! step 2: back substitution
               dU(top(i)+1) = d_matrix(top(i)+1)/b_matrix(top(i)+1)

               do ksn=top(i),1,-1
                  dU(ksn) = (d_matrix(ksn) - c_matrix(ksn) * dU(ksn+1))/b_matrix(ksn)
               enddo

               do ksn = 1,top(i)+1
                  ddU(ksn) = dU(ksn) - ddU(ksn)
                  if (maxddU .lt. abs(ddU(ksn)) .or. ksn .eq. 1) then
                     maxddU = abs(ddU(ksn))
                  endif
                  U(ksn) = U(ksn) + ddU(ksn)
               enddo

               if (maxddU < tol) then
                  exit iterloop
               endif

            enddo iterloop ! iteration

            do ksn=1,top(i)+1
               t_sn_n(i,ksn) = U(ksn)
            enddo


            do ksn=1,top(i)
               t_sn(i,ksn) = 0.5_wp * (t_sn_n(i,ksn) + t_sn_n(i,ksn+1))
            enddo

            t_sn_sfc(i) = t_sn_n(i,top(i)+1)
            ! shflx_sn(i) = alpha_shf * ( t(i) - t_sn_sfc(i) )


         endif ! end snow on the ground

      enddo ! end of i

   end subroutine heat_equation_wrapper
!------------------------------------------------------------------------------
! end of module mo_nix_heat_equation
!------------------------------------------------------------------------------

end module mo_nix_heat_equation
