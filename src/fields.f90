module fields

! ------------------------------------------------------------------------------
!
! description: this module declares all fields, which are used by more than one
!              module, i.e. non local arrays. all fields are declared as
!              allocatble arrays, then allocated in the main setup and
!              deallocted during clean-up.
!
!
! current code owner: sascha bellaire
!
! e-mmail:  sascha.bellaire@gmail.com
!
! history:
! version    date       name
! ---------- ---------- ----
!
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! + begin declarations
! ------------------------------------------------------------------------------

! -------------------
! modules used
! ------------------

   use mo_kind, only :  wp           ! kind-type parameter for real variables

! ------------------------------------------------------------------------------
! - end declarations
! ------------------------------------------------------------------------------

   implicit none

! ------------------------------------------------------------------------------
! + begin module procedures
! ------------------------------------------------------------------------------

   ! -------------------
   ! input
   ! ------------------
   real (kind = wp), dimension(:,:), allocatable :: &
      prr_con    , & ! precipitation rate of rain, convective        (kg/m2*s)
      prs_con    , & ! precipitation rate of snow, convective        (kg/m2*s)
      prr_gsp    , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
      prs_gsp    , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
      prg_gsp    , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
      u          , & ! zonal wind speed                              (m/s)
      v          , & ! meridional wind speed                         (m/s)
      t          , & ! temperature                                   (k)
      qv         , & ! specific water vapor content                  (kg/kg)
      ps         , & ! surface pressure                              ( pa  )
      t_so       , &
      iswr, ilwr

   real (kind = wp), dimension(:, :), allocatable :: &
      t_sn                           , & ! snow layer temperature (main level)
      theta_i                          , & ! volumetric ice content
      theta_w                          , & ! volumetric water content
      theta_a                          , & ! volumetric air content
      dzm_sn                               ! snow layer thickness

   real (kind=wp), dimension(:,:), allocatable :: t_sn_n

   integer, dimension(:), allocatable :: &
      top_sn

   real (kind = wp), dimension(:), allocatable :: &
      hn_sn                            , & ! new snow amount
      t_sn_sfc                         , & ! snow surface temperature                    [ k ]
      swflx_sn_net                     , & ! short wave radiation - net                  [w/m**2]
      swflx_sn_dn                      , & ! short wave radiation - downward             [w/m**2]
      swflx_sn_up                          ! short wave radiation - upward               [w/m**2]

   real (kind = wp), dimension(:,:), allocatable:: &
      swflx_sn_abs                         ! inout short wave radiation - absorbed             [w/m**2]

   real (kind = wp), dimension(:), allocatable :: &
      lwflx_sn_up                       , & ! long wave radiation  - upward               [w/m**2]
      lwflx_sn_dn                       , & ! long wave radiation  - downward             [w/m**2]
      alpha_sn                         , & ! surface albedo (snow)                       [-]
      shflx_sn                         , & ! turbulent sensible heat flux                [w/m**2]
      lhflx_sn                         , & ! turbulent latent heat flux                  [w/m**2]
      tch_sn                           , & ! transfer coefficient                        [    ]
      for_sn                           , & ! total atmospheric forcing at snow surface   [w/m**2]
      qv0_sn                               ! specific humidity ar snow surface

   real (kind = wp), dimension(:,:), allocatable :: &
      hcap_sn                          , & ! snow layer heat capacity                    [   ]
      hcon_sn                          , & ! snow layer heat conductivity                [   ]
      hdif_sn                          , & ! snow layer heat diffusivity                 [   ]
      rho_sn                           , & ! snow layer density                          [kg/m**3]
      mass_sn                          , & ! snow layer mass                             [kg]
      hm_sn                            , & ! height (from bottom) of main snow level     [ m ]
      zm_sn                                ! depth (from top) of main snow level         [ m ]

   real (kind = wp), dimension(:), allocatable :: &
!
      swe_sn                           , & ! snow water equivialent                      [   ]
      runoff_sn                            ! total melt runoff                           [   ]

   real (kind = wp), dimension(:), allocatable :: &
      h_snow                               ! snow height                                 [ m ]

   !  real  (kind=wp),     allocatable ::           &

   !     t            (:)   , &   ! air temperature                                       (k)
   !     ps           (:)   , &   ! surface pressure                                      (k)
   !     qv           (:)   , &   ! specific humidity                                     (p)
   !     rh           (:)   , &   ! relative humidity                                     (%)
   !     uv           (:)   , &   ! wind speed                                            (m/s)
   !     t_so         (:)   , &   ! soil temperature                                      (k)

   !     tot_prec     (:)   , &   ! precipitations rate (water equivialent)               (mm/s)

   !     swdir_s      (:)   , &   ! incoming direct shortwave radiation                          (w m-2)
   !     swdifd_s     (:)   , &   ! incoming fiffuse shortwave radiation
   !     lwd_s        (:)         ! long wave radiation (down)                            (w m-2)


   !  ! -------------------
   !  ! output
   !  ! ------------------

   !  real  (kind=wp),     allocatable ::           &

   !  ! snow cover properties
   !     t_sn_now         (:,:)   , &   ! multi-layer snow temperature                          (k)

   !     theta_i_now      (:,:)   , &   ! volumetric ice content                                (-)
   !     theta_w_now      (:,:)   , &   ! volumetric ice content                                (-)
   !     theta_a_now      (:,:)   , &   ! volumetric ice content                                (-)

   !     rho_sn       (:,:)   , &   ! layer density                                         (kg m-3)

   !     hm_sn        (:,:)   , &   ! main level height                                    (m)
   !     zm_sn        (:,:)   , &   ! main level depth                                      (m)

   !     dzm_sn_now       (:,:)   , &  ! thickness between main levels                          (m)

   !     m_sn         (:,:)   , &   ! layer mass                                            (kg)

   !     hcap_sn      (:,:)   , &   ! heat capacity                                         ( )
   !     hcon_sn      (:,:)   , &   ! heat conductivity                                     ( )
   !     hdif_sn      (:,:)   , &   ! heat diffusivity                                      ( )

   !     h_snow       (:)      , &   ! snow cover height                                     (m)

   !     t_sn_sfc     (:)      , &   ! snow surface temperature                              (k)

   !     swabs_sn       (:) , &   ! absorbed short wave radiation (each layer)            (w m-2)

   !     lwu_s           , &
   !     lh_sn           , &
   !     sh_sn           , &

   !  ! meteorological

   !     sw_net    (:)         , &   ! net. short wave radiation                             (w m-2)
   !     sw_u      (:)         , &   ! upward shortwave radiation                            (w m-2)
   !     sw_d      (:)         , &
   !     lw_u      (:)         , &   ! long wave radiation (upwelling)                       (w m-2)
   !     lw_d      (:)         , &

   !     lh        (:)         , &   ! latent heat                                           (w m-2)
   !     sh        (:)         , &   ! sensible heat                                         (w m-2)

   !     albedo    (:)         , &   ! albedo                                                ( - )
   !     tch       (:)         , &   ! transfer coefficient

   !     hn_sn_now (:)         , &   ! new snow amounts                                      (m)

   !     zsnow_rate         , &

   !     zuv                , &

   !  ! snow cover properties
   !     runoff             , &   ! melt water runoff                                     (mm)
   !     rho_hn             , &   ! density of new snow                                   (kg m-3)


   !  ! old values of certain key paramters
   !     h_snow_old         , &

   !     zm_sn_old(:)       , &

   !     theta_i_old(:)     , &
   !     theta_a_old(:)     , &
   !     theta_w_old(:)

   !  integer,              allocatable ::           &

   !     top                , &   ! index of top snow layer
   !     bot                , &   ! index of bottom layer, i.e. max_n_snow

   !     n_snow_now               ! number of snow layer                                  (-)


   !  real (kind=wp),       allocatable :: &

   !     zm          (:)         , &  ! depth of main levels

   !     hcon        (:)         , & ! heat conductivity
   !     hcap        (:)         , & ! heat capacity
   !     hdif        (:)         , & ! heat diffusivity

   !     rho         (:)         , &  ! density

   !     t_sol       (:)         , &  ! temperature of the snow/substrate column

   !     sw_abs      (:)         , &  ! absorbed short-wave radiaton in each layer

   !     alpha       (:)         , &  ! utility variables for building and solving the tri-diagonal matrix
   !     gamma_sol   (:)         , &  !

   !     a           (:)         , &  !
   !     b           (:)         , &  !
   !     c           (:)         , &  !
   !     d           (:)         , &  !
   !     e           (:)              ! final snow layer temperature



   !  real (kind= wp), allocatable  :: &

   !     w_res      (:)              ! effective residual water content


! ------------------------------------------------------------------------------
! - end module procedures
! ------------------------------------------------------------------------------

! =============================================================================
! - end module for ...
! ==============================================================================

end module fields


