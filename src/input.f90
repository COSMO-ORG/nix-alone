module input
! ------------------------------------------------------------------------------
! description: this module reads the meteorological input and prepare it for
!              further use.
! ------------------------------------------------------------------------------

   use mo_kind, only : wp           ! kind-type parameter for real variables
   use config , only : nsteps, nvec
   use fields

   implicit none

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
         nrows     , &   ! count for number of rows in file          (-)
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
      nrows = 0
      do
         read(20, *, iostat=iostatus)
         if(iostatus/=0) then ! to avoid end of file error.
            exit
         else
            nrows = nrows + 1
         end if
      end do
      close(20)


      ! -------------------
      ! allocate array's
      ! -------------------

      nsteps = nrows - 1

      allocate (  tmp     (9,nrows) )

      allocate ( prr_con  (nvec,nrows)) ; prr_con = 0.0_wp
      allocate ( prs_con(nvec,nrows)) ; prs_con = 0.0_wp
      allocate ( prr_gsp(nvec,nrows)) ; prr_gsp = 0.0_wp
      allocate ( prs_gsp(nvec,nrows)) ; prs_gsp = 0.0_wp
      allocate ( prg_gsp(nvec,nrows)) ; prg_gsp = 0.0_wp
      allocate ( u(nvec,nrows))
      allocate ( v(nvec,nrows))
      allocate ( t(nvec,nrows))
      allocate ( qv(nvec,nrows))
      allocate ( ps(nvec,nrows))
      allocate ( t_so(nvec,nrows))
      allocate (iswr(nvec,nrows))
      allocate (ilwr(nvec,nrows))
      ! -------------------
      ! read data into array
      ! -------------------

      ! open unit nad read files
      open(unit=20, status="old", file=file)
      do i = 1, nrows, 1
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

         do j = 1, nrows
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



































