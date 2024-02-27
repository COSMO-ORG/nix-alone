module allocation

! ------------------------------------------------------------------------------
!
! description: this module contains subroutine:
!
!           - allocate_fields for allocation of all global field
!           - deallocate_fields for de-allocation of all global fields
!
! ------------------------------------------------------------------------------

   use mo_kind, only : wp   ! kind-type parameter for real variables

   use mo_nix_config, only :   &
      ke_snow                , &     ! maximum number of snow layers
      ke_soil                , &     ! number of soil layers
      nvec

   use fields

! ------------------------------------------------------------------------------
! - end declarations
! ------------------------------------------------------------------------------

   implicit none

contains

! ------------------------------------------------------------------------------
! + begin subroutine alloc_fields
! ------------------------------------------------------------------------------

   subroutine allocate_fields()

! ------------------------------------------------------------------------------
!
! description:
!   this routine allocates space for all global fields and initializes them
!
! method:
!   all fields are allocated with the allocate statemenet
!
! -----------------------------------------------------------------------------

      ! -------------------
      ! output
      ! -------------------
      ! write(*,*) 'nvec: ', nvec

      allocate ( t_sn     (nvec,1:ke_snow) ) ; t_sn      = 0.0_wp
      allocate ( theta_i  (nvec,1:ke_snow) ) ; theta_i   = 0.0_wp
      allocate ( theta_w  (nvec,1:ke_snow) ) ; theta_w   = 0.0_wp
      allocate ( theta_a  (nvec,1:ke_snow) ) ; theta_a   = 0.0_wp
      allocate ( dzm_sn   (nvec,1:ke_snow) ) ; dzm_sn    = 0.0_wp

      allocate ( t_sn_n   (nvec,1:ke_snow) ) ; t_sn_n    = 0.0_wp

      allocate ( top_sn   (nvec) )           ; top_sn = 0

      allocate ( hn_sn        (nvec) ) ; hn_sn = 0.0_wp
      allocate ( t_sn_sfc     (nvec) ) ; t_sn_sfc  = 273.15_wp
      allocate ( swflx_sn_net (nvec) ) ; swflx_sn_net = 0.0_wp
      allocate ( swflx_sn_dn  (nvec) ) ; swflx_sn_dn = 0.0_wp
      allocate ( swflx_sn_up  (nvec) ) ; swflx_sn_up = 0.0_wp

      allocate ( swflx_sn_abs (nvec,1:ke_snow) ) ; swflx_sn_abs = 0.0_wp

      allocate ( lwflx_sn_up (nvec) ) ; lwflx_sn_up = 0.0_wp
      allocate ( lwflx_sn_dn (nvec) ) ; lwflx_sn_dn = 0.0_wp
      allocate ( alpha_sn    (nvec) ) ; alpha_sn    = 0.3_wp
      allocate ( shflx_sn    (nvec) ) ; shflx_sn    = 0.0_wp
      allocate ( lhflx_sn    (nvec) ) ; lhflx_sn    = 0.0_wp
      allocate ( tch_sn      (nvec) ) ; tch_sn    = 0.0_wp
      allocate ( for_sn      (nvec) ) ; for_sn    = 0.0_wp
      allocate ( qv0_sn      (nvec) ) ; qv0_sn    = 0.0_wp


      allocate ( hcap_sn (nvec,1:ke_snow) ) ; hcap_sn   = 0.0_wp
      allocate ( hcon_sn (nvec,1:ke_snow) ) ; hcon_sn   = 0.0_wp
      allocate ( hdif_sn (nvec,1:ke_snow) ) ; hdif_sn   = 0.0_wp
      allocate ( rho_sn  (nvec,1:ke_snow) ) ; rho_sn    = 0.0_wp
      allocate ( mass_sn (nvec,1:ke_snow) ) ; mass_sn   = 0.0_wp
      allocate ( hm_sn   (nvec,1:ke_snow) ) ; hm_sn     = 0.0_wp
      allocate ( zm_sn   (nvec,1:ke_snow) ) ; zm_sn     = 0.0_wp


      allocate ( swe_sn    (nvec) ) ; swe_sn    = 0.0_wp
      allocate ( runoff_sn (nvec) ) ; runoff_sn = 0.0_wp
      allocate ( h_snow    (nvec) ) ; h_snow    = 0.0_wp

! ------------------------------------------------------------------------------
! - end subroutine alloc_fields
! ------------------------------------------------------------------------------

   end subroutine allocate_fields

! ------------------------------------------------------------------------------
! + begin subroutine alloc_fields
! ------------------------------------------------------------------------------

   subroutine deallocate_fields()

! ------------------------------------------------------------------------------
!
! description:
!   this routine deallocates all global fields.
!
! method:
!   all fields are allocated with the deallocate statemenet
!
! -----------------------------------------------------------------------------

      ! -------------------
      ! output
      ! -------------------
      deallocate ( t_sn    )
      deallocate ( theta_i )
      deallocate ( theta_w )
      deallocate ( theta_a )
      deallocate ( dzm_sn  )

      deallocate ( t_sn_n )

      deallocate ( top_sn )

      deallocate ( hn_sn        )
      deallocate ( t_sn_sfc     )
      deallocate ( swflx_sn_net )
      deallocate ( swflx_sn_dn  )
      deallocate ( swflx_sn_up  )

      deallocate ( swflx_sn_abs )

      deallocate ( lwflx_sn_up)
      deallocate ( lwflx_sn_dn)
      deallocate ( alpha_sn   )
      deallocate ( shflx_sn   )
      deallocate ( lhflx_sn   )
      deallocate ( tch_sn     )
      deallocate ( for_sn     )
      deallocate ( qv0_sn     )


      deallocate ( hcap_sn )
      deallocate ( hcon_sn )
      deallocate ( hdif_sn )
      deallocate ( rho_sn  )
      deallocate ( mass_sn )
      deallocate ( hm_sn   )
      deallocate ( zm_sn   )


      deallocate ( swe_sn    )
      deallocate ( runoff_sn )
      deallocate ( h_snow    )



! ------------------------------------------------------------------------------
! - end subroutine alloc_fields
! ------------------------------------------------------------------------------

   end subroutine deallocate_fields





! =============================================================================
! - end module for ...
! ==============================================================================

end module allocation


