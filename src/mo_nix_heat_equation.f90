module mo_nix_heat_equation

   use mo_kind,                    only: wp
   use mo_physical_constants,      only: stbo
   use mo_nix_constants, only: eps_div
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
         ke_snow             ! < number of snow layers

      integer, dimension(nvec), intent(in) :: &
         top                    ! top layer index

      ! vs: intent for t_sn !?
      real (kind=wp), dimension(nvec,ke_snow), intent(inout) :: &
         dzm_sn              , & ! < snow layer depth
         hcon_sn            , & ! <            conductivity
         hcap_sn            , & ! <            capacity
         hdif_sn            , & ! <            difusivity
         t_sn               , & ! <            temperature
         swflx_sn_abs       , &
         theta_w ! <            absorbed short-wave radiation

      real (kind=wp), dimension(nvec), intent(inout) :: &
         lwflx_sn_dn        , &
         lwflx_sn_up        , &
         lhflx_sn, shflx_sn


      real (kind=wp), dimension(nvec,ke_snow), intent(inout) :: rho_sn

      real (kind=wp), dimension(nvec,ke_snow+1), intent(inout) :: &
         t_sn_n

      real (kind=wp), dimension(nvec), intent(in) :: &
         hcon_so            , & ! < soil layer conductivity
         t_so                   ! <            temperature

      real (kind = wp), intent(in)  ::  &
         dt                     ! time step

      ! vs: intent for t_sn_sfc should be in-out !?
      real (kind = wp), dimension(nvec), intent(inout)  ::  &
      &   t_sn_sfc               ! < snow surface temperature

      real (kind = wp), dimension(nvec), intent(in)  ::  &
      &   tch_sn

      real (kind=wp), dimension(nvec), intent(in) :: t

      !!! local variables
      integer :: l_top, ksn, counter, i
      real(kind=wp), dimension(ke_snow+1) :: a_matrix,b_matrix,c_matrix,d_matrix,tmp_t
      real(kind=wp) :: emiss, t_emiss, coeff, delta_rad
      real(kind=wp) :: alpha_solver_down,alpha_solver_up
      real(kind=wp) :: beta_solver_down,beta_solver_up
      real(kind=wp) :: gamma_sh, gamma_soil


      do i = ivstart, ivend
         l_top = top(i)
         if(l_top .gt. 1) then  !!!snow on the ground

            ! PREP
            emiss   = lwflx_sn_dn(i)/(stbo*t(i)*t(i)*t(i)*t(i))
            t_emiss = sqrt(sqrt(emiss)) * t(i)
            delta_rad = stbo * (t_emiss + t_sn_sfc(i)) * (t_emiss * t_emiss + t_sn_sfc(i) * t_sn_sfc(i))

            ! exchange coefficient for sensible heat flux
            if ( abs( t(i) - t_sn_sfc(i) ) < eps_div  ) then
               gamma_sh = 1.0
            else
               gamma_sh = shflx_sn(i) / ( t(i) - t_sn_sfc(i) )
            endif

            gamma_soil = (hcon_so(i)/0.005_wp) ! VS: length of the top soil layer

            counter = 0
            do ksn = l_top+1,1,-1
               counter = counter+1

               if(counter .eq. 1) then ! ... top node ! neumann bc ! KSN = l_top + 1

                  alpha_solver_up   = (dzm_sn(i,ksn-1) * rho_sn(i,ksn-1) * hcap_sn(i,ksn-1)) / ( 6.0_wp * dt )
                  beta_solver_up    = ( hcon_sn(i,ksn-1) ) * (1.0_wp/dzm_sn(i,ksn-1))

                  a_matrix(counter)     =  0.0_wp
                  b_matrix(counter)     =  2.0_wp * alpha_solver_up + beta_solver_up + delta_rad + gamma_sh ! delta_rad for LW ! Gamma_sh for Sensible
                  c_matrix(counter)     =  alpha_solver_up - beta_solver_up

                  d_matrix(counter)     =  2.0_wp * alpha_solver_up * t_sn_n(i,ksn) + alpha_solver_up * t_sn_n(i,ksn-1) &
                  &                        + delta_rad * t_emiss &
                  &                        + lhflx_sn(i) & ! for latent heat
                  &                        + gamma_sh * t(i)  & ! for sensible
                  &                        + swflx_sn_abs(i,ksn-1)


               elseif (counter .eq. l_top+1) then ! ... bottom node ! neumann BC ! KSN = 1

                  alpha_solver_up   = (dzm_sn(i,ksn) * rho_sn(i,ksn) * hcap_sn(i,ksn)) / ( 6.0_wp * dt )
                  beta_solver_up    = ( hcon_sn(i,ksn) ) * (1.0_wp/dzm_sn(i,ksn))


                  a_matrix(counter) = alpha_solver_up - beta_solver_up
                  b_matrix(counter) = 2.0_wp * alpha_solver_up + beta_solver_up + gamma_soil
                  c_matrix(counter) = 0.0_wp

                  d_matrix(counter) = alpha_solver_up * t_sn_n(i,ksn+1) + 2.0_wp * alpha_solver_up * t_sn_n(i,ksn) &
                  &                   + gamma_soil * t_so(i)

               else  ! middle nodes

                  alpha_solver_up  = ( dzm_sn(i,ksn) * rho_sn(i,ksn) * hcap_sn(i,ksn) ) / ( 6.0_wp * dt )
                  beta_solver_up   = ( hcon_sn(i,ksn) ) * (1.0_wp/dzm_sn(i,ksn))

                  alpha_solver_down = ( dzm_sn(i,ksn-1) * rho_sn(i,ksn-1) * hcap_sn(i,ksn-1) )  / ( 6.0_wp * dt )
                  beta_solver_down  = ( hcon_sn(i,ksn-1) ) * (1.0_wp/dzm_sn(i,ksn-1))

                  a_matrix(counter) = alpha_solver_up - beta_solver_up
                  b_matrix(counter) = 2.0_wp * (alpha_solver_up + alpha_solver_down) + beta_solver_up + beta_solver_down
                  c_matrix(counter) = alpha_solver_down - beta_solver_down

                  d_matrix(counter) = alpha_solver_up * t_sn_n(i,ksn+1) +  &
                  &                                2.0_wp * (alpha_solver_up + alpha_solver_down) * t_sn_n(i,ksn) + &
                  &                                alpha_solver_down * t_sn_n(i,ksn-1)

               endif ! if block for splitting between top, middle and bottom nodes
            enddo


            ! ------------------------------------------------------------
            ! Solve the system - Thomas Algorithm
            ! ------------------------------------------------------------

            ! step 1: forward elimination
            do ksn=2,l_top+1
               coeff  = a_matrix(ksn)/b_matrix(ksn-1)
               b_matrix(ksn) = b_matrix(ksn) - coeff * c_matrix(ksn-1)
               d_matrix(ksn) = d_matrix(ksn) - coeff * d_matrix(ksn-1)
            enddo

            tmp_t = 0.0_wp
            ! step 2: back substitution
            tmp_t(l_top+1) = d_matrix(l_top+1)/b_matrix(l_top+1)

            do ksn=l_top,1,-1
               !write(*,*) ksn,d_matrix(ksn),c_matrix(ksn), t_sn_now(i,ksn),b_matrix(ksn)
               tmp_t(ksn) = (d_matrix(ksn) - c_matrix(ksn) * tmp_t(ksn+1))/b_matrix(ksn)
            enddo


            counter=0
            do ksn=l_top+1,1,-1
               counter=counter+1
               t_sn_n(i,ksn) = tmp_t(counter)
            enddo


            do ksn=1,l_top
               t_sn(i,ksn) = 0.5_wp * (t_sn_n(i,ksn) + t_sn_n(i,ksn+1))
            enddo

            t_sn_sfc(i) = t_sn_n(i,l_top+1)
            ! shflx_sn(i) = gamma_sh * ( t(i) - t_sn_sfc(i) )




         endif ! end snow on the ground

      enddo ! end of i










   end subroutine heat_equation_wrapper
!------------------------------------------------------------------------------
! end of module mo_nix_heat_equation
!------------------------------------------------------------------------------

end module mo_nix_heat_equation
