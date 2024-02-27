module input
! ------------------------------------------------------------------------------
! description: this module reads the meteorological input and prepare it for
!              further use.
! ------------------------------------------------------------------------------

   use mo_kind, only : wp           ! kind-type parameter for real variables
   use mo_nix_config , only : nvec
   use fields

   implicit none

   integer   :: &
      nsteps                        ! count for number of rows in file          (-)

contains

! ------------------------------------------------------------------------------
! + begin subroutine read_input
! ------------------------------------------------------------------------------
   subroutine read_input()
! ------------------------------------------------------------------------------
! description:
!   this routine reads required meteorological input for the snow cover scheme
! ------------------------------------------------------------------------------

      integer   :: &
         i , j          , &   ! loop index
         iostatus

      real (kind=wp), allocatable  :: tmp(:,:)
      character(*), parameter :: file = "./inp/WFJ_forcing.txt"

!    file = "./inp/kenda.txt"
!    file = "./inp/dmo.txt"
!    file = "./inp/snp.txt"

      ! -------------------
      ! count numbers of rows
      ! -------------------

      open(unit=20, status = "old", file=file)
      nsteps = 0
      do
         read(20, *, iostat=iostatus)
         if(iostatus/=0) then ! to avoid end of file error.
            exit
         else
            nsteps = nsteps + 1
         end if
      end do
      close(20)


      ! -------------------
      ! allocate array's
      ! -------------------

      allocate (  tmp     (9,nsteps) )

      allocate ( prr_con  (nvec,nsteps)) ; prr_con = 0.0_wp
      allocate ( prs_con(nvec,nsteps)) ; prs_con = 0.0_wp
      allocate ( prr_gsp(nvec,nsteps)) ; prr_gsp = 0.0_wp
      allocate ( prs_gsp(nvec,nsteps)) ; prs_gsp = 0.0_wp
      allocate ( prg_gsp(nvec,nsteps)) ; prg_gsp = 0.0_wp
      allocate ( u(nvec,nsteps))
      allocate ( v(nvec,nsteps))
      allocate ( t(nvec,nsteps))
      allocate ( qv(nvec,nsteps))
      allocate ( ps(nvec,nsteps))
      allocate ( t_so(nvec,nsteps))
      allocate (iswr(nvec,nsteps))
      allocate (ilwr(nvec,nsteps))
      ! -------------------
      ! read data into array
      ! -------------------

      ! open unit nad read files
      open(unit=20, status="old", file=file)
      do i = 1, nsteps, 1
         read(20, *, iostat=iostatus ) tmp(:,:)
         if(iostatus/=0) then ! to avoid end of file error.
            exit
         endif
      end do

      ! close unit
      close(20)

      ! -------------------
      ! assign data to array
      ! -------------------
      do i = 1, nvec
         t        (i,:)    = tmp(1,:)   ! air temperature
         ps       (i,:)    = tmp(2,:)   ! pressure
         qv       (i,:)    = tmp(3,:)   ! specific humidity
         u        (i,:)    = tmp(4,:)   ! wind speed
         v        (i,:)    = 0.0_wp
         iswr     (i,:)    = tmp(5,:) + tmp(6,:)  ! incoming short-wave radiation ( direct + diffuse)
         ilwr     (i,:)    = tmp(7,:)   ! incoming long-wave radiation
         !tot_prec (i,:)    = tmp(8,:)   ! total precipitation
         t_so     (i,:)    = tmp(9,:)   ! soil temperature

         do j = 1, nsteps
            if ( t(i,j) > 275.15) then
               prr_con(i,j) = tmp(8,j)
            else
               prs_con(i,j) = tmp(8,j)
            end if
         end do

      end do

! ------------------------------------------------------------------------------
! - end subroutine read_input
! ------------------------------------------------------------------------------

   end subroutine read_input


! =============================================================================
! - end module for ...
! ==============================================================================



end module input



































