module output
! ------------------------------------------------------------------------------
!
! description: this module handles the output of the snow cover scheme
!
! ------------------------------------------------------------------------------

   use mo_kind, only : wp           ! kind-type parameter for real variables
   use fields
   use config

   implicit none

contains


   subroutine write_output(n)

! ------------------------------------------------------------------------------

      integer,    intent(in)   :: n

!       ! ----------------------
!       ! profile information
!       ! ----------------------

!       ! open file
!       open(22, file = './out/pro.txt')

! !      write(22,'(a7,   2i12.5,a,/)',    advance = 'yes') 'idx'        , n, top
!       write(22,'(a7,  1000f12.3,a,/)',   advance = 'yes') 'hm_sn'      , hm_sn(:)
!       write(22,'(a7,  1000f12.3,a,/)',   advance = 'yes') 'zm_sn'      , zm_sn(:)
!       write(22,'(a7,  1000f12.3,a,/)',   advance = 'yes') 'dz'         , dzm_sn_now(:)
!       write(22,'(a7,  1000f12.3,a,/)',   advance = 'yes') 't_sn'       , t_sn_now(:)
!       write(22,'(a,   1000f12.3,a,/)',   advance = 'yes') 'theta_i'    , theta_i_now(:)
!       write(22,'(a,   1000f12.3,a,/)',   advance = 'yes') 'theta_w'    , theta_w_now(:)
!       write(22,'(a,   1000f12.3,a,/)',   advance = 'yes') 'theta_a'    , theta_a_now(:)
!       write(22,'(a7,  1000f12.3,a,/)',   advance = 'yes') 'rho_sn'     , rho_sn(:)
!       write(22,'(a7,  1000f12.3,a,/)',   advance = 'yes') 'm_sn'       , m_sn(:)

!       write(22,'(a,   1000f12.3,a,/)',   advance = 'yes') 'hcap_sn'    , hcap_sn(:)
!       write(22,'(a,   1000f12.3,a,/)',   advance = 'yes') 'hcon_sn'    , hcon_sn(:)
!       write(22,'(a,   1000f12.3,a,/)',   advance = 'yes') 'hdif_sn'    , hdif_sn(:)


!       if(n .eq. nsteps) then
!          close(22)
!       endif


!       ! ----------------------
!       ! meteorology
!       ! ----------------------

!       ! open file
!       open(44, file = './out/met.txt')

!       write(44,'(2i6, 7f18.3,a,/)',   advance = 'yes')  n, top, h_snow, t_sn_sfc, swabs_sn(1), lwd_s(n), lwu_s, lh_sn, sh_sn


!       if(n .eq. nsteps) then
!          close(44)
!       endif


! ! ------------------------------------------------------------------------------
! ! - end subroutine write_output
! ! ------------------------------------------------------------------------------

   end subroutine write_output

! =============================================================================
! - end module for ...
! ==============================================================================

end module output


