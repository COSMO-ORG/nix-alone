module config

! ------------------------------------------------------------------------------
! description: this one contains all required model configurations.
! ------------------------------------------------------------------------------

   use mo_kind, only : wp           ! kind-type parameter for real variables

   implicit none

   public ! all constants and variables in this module are public


   integer ::  &

      nsteps    = 1270080                 , & ! total number of tiemsteps
      ke_snow   = 10                      , & ! total number of snow layers
      ke_soil   = 1                       , & ! total number of soil layers
      nvec      = 1
   real (kind=wp)                   ::    &

      zdt      = 900.0_wp         , &    ! model timestep : note this should correspond to the meteotimestep     (s)
      sdt      = 20.0_wp                  ! internal timestep for heat equation

   real (kind=wp), parameter        ::    &

      min_height_layer = 0.002_wp       , &  ! minimum layer thickness
      max_height_layer = 0.01_wp              ! maximum layer thickness

   real (kind=wp), parameter        ::    &

      height_meteo_values = 10.0_wp       , & ! height of meteo values above ground soil                     (m)
      altitude            = 2540.0_wp     , & ! elevation above sea level                                    (m)

      roughness_length    = 0.002_wp          ! roughness length                                             (m)

! ------------------------------------------------------------------------------
! - end declarations
! ------------------------------------------------------------------------------


! =============================================================================
! - end module
! ==============================================================================

end module config


