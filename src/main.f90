program  main

   use mo_kind,  only : wp  ! definitions of kind parameter double vs. single precision
   use fields                       ! contains all required global fields
   use allocation                   ! allocates all global fields
   use input                        ! handles the input data
   use config                       ! model configuration
   use output                       ! handles the output
   use mo_nix_constants             ! physical constants
   use mo_nix_main
   use mo_nix_init

   integer :: n
! -------------------
! read meteorological input
! -------------------
   call read_input()
! -------------------
! allocate all global arrays.
!-------------------
   call allocate_fields()
! ---------------------------------------------------------------------------------
! + section i: initializations
! ---------------------------------------------------------------------------------

   call nix_init(nvec      , &
   &         ke_snow   , &
   &         1   , &
   &         1   , &
   &         t_sn      , &
   &         theta_i   , &
   &         theta_w   , &
   &         theta_a   , &
   &         dzm_sn    , &
   &         hn_sn     , &
   &         top_sn    , &
   &         h_snow      )

! begin loop over number of timesteps
   do n = 1,nsteps,1

      !write(*,*) 'n: ', n
      call nix_core(             &
      ! Utility variables (IN)
      &     nvec                             , &
      &     1                                , &
      &     1                                , &
      &     1                                , &
      &     ke_soil                          , &
      &     ke_snow                          , &
      &     2                                , &
      &     zdt                               , &
      ! NIX input variables (IN)
      &     prr_con (:,n)                    , &
      &     prs_con (:,n)                    , &
      &     prr_gsp (:,n)                    , &
      &     prs_gsp (:,n)                    , &
      &     prg_gsp (:,n)                    , &
      &     u (:,n)                          , &
      &     v (:,n)                          , &
      &     t (:,n)                          , &
      &     qv (:,n)                         , &
      &     ps (:,n)                         , &
      ! NIX prognostic variables (INOUT)
      &     t_sn                             , &
      &     t_sn_n                           , &
      &     theta_i                          , &
      &     theta_w                          , &
      &     theta_a                          , &
      &     dzm_sn                           , &
      &     hn_sn                            , &
      &     top_sn                           , &
      ! NIX specific fluxes and related variables
      &     t_sn_sfc                         , &
      &     swflx_sn_net                     , &
      &     iswr (:,n)                       , &
      &     swflx_sn_up                      , &
      &     swflx_sn_abs                     , &
      &     lwflx_sn_up                      , &
      &     ilwr (:,n)                       , &
      &     alpha_sn                         , &
      &     shflx_sn                         , &
      &     lhflx_sn                         , &
      &     tch_sn                           , &
      &     for_sn                           , &
      &     qv0_sn                           , &
      ! Heat equation
      &     hcap_sn                          , &
      &     hcon_sn                          , &
      &     hdif_sn                          , &
      ! Additional fields
      &     rho_sn                           , &
      &     mass_sn                          , &
      &     hm_sn                            , &
      &     zm_sn                            , &
      &     swe_sn                           , &
      &     runoff_sn                        , &
      &     h_snow                           , &
      ! Soil properties
      &     t_so (:,n)                       )

      call write_output(nvec, 1, 1, n, 48                  , &
        &                 top_sn, ke_snow, dzm_sn, rho_sn  , &
        &                 theta_i, theta_w, theta_a        , &
        &                 t_sn, t_sn_n, hn_sn, t , alpha_sn)

   enddo ! nsteps

   ! deallocate arrays
   call deallocate_fields()

   ! print *, 'done, done!'

end program  main

