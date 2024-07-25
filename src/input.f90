module input
! ------------------------------------------------------------------------------
! description: this module reads the meteorological input and prepare it for
!              further use.
! ------------------------------------------------------------------------------

   use mo_kind, only : wp           ! kind-type parameter for real variables
   use mo_nix_config, only :   &
      ke_snow                , &     ! maximum number of snow layers
      ke_soil                , &     ! number of soil layers
      nvec
   use fields
   use allocation

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

      ! These two variables are needed to prepare the forcing data, but they are not part of the global variables
      real (kind=wp), allocatable  :: &
         iswr_d(:,:),                 &
         tot_prec(:,:)

      character(*), parameter :: file = "./inp/icon_15min_2021.inp"

      ! -------------------
      ! count numbers of rows
      ! -------------------

      open(unit=20, status = "old", file=file)
      read(20,*)  ! Skip one header line

      nsteps = 0

      do
         read(20, *, iostat=iostatus)
         if(iostatus/=0) then ! to avoid end of file error.
            exit
         else
            nsteps = nsteps + 1
         end if
      end do
      rewind(20)


      ! -------------------
      ! allocate array's
      ! -------------------
      allocate ( timestamp(nvec,nsteps) ) ; timestamp = ""
      allocate ( tot_prec(nvec,nsteps) ) ; tot_prec = 0.0_wp
      allocate ( prr_con(nvec,nsteps) ) ; prr_con = 0.0_wp
      allocate ( prs_con(nvec,nsteps) ) ; prs_con = 0.0_wp
      allocate ( prr_gsp(nvec,nsteps) ) ; prr_gsp = 0.0_wp
      allocate ( prs_gsp(nvec,nsteps) ) ; prs_gsp = 0.0_wp
      allocate ( prg_gsp(nvec,nsteps) ) ; prg_gsp = 0.0_wp
      allocate ( u(nvec,nsteps) )
      allocate ( v(nvec,nsteps) )
      allocate ( t(nvec,nsteps) )
      allocate ( qv(nvec,nsteps) )
      allocate ( ps(nvec,nsteps) )
      allocate ( t_so(nvec,nsteps) )
      allocate ( iswr(nvec,nsteps) )
      allocate ( iswr_d(nvec,nsteps) )
      allocate ( ilwr(nvec,nsteps) )

      ! -------------------
      ! read data into array
      ! -------------------

      ! open unit and read files
      read(20,*)  ! Skip one header line
      do i = 1, nsteps, 1
         read (20, *, iostat=iostatus ) timestamp(nvec, i), &
                                        t(nvec, i),         &
                                        ps(nvec, i),        &
                                        qv(nvec, i),        &
                                        u(nvec, i),         &
                                        v(nvec, i),         &
                                        iswr(nvec, i),      &
                                        iswr_d(nvec, i),    &
                                        ilwr(nvec, i),      &
                                        tot_prec(nvec, i),  &
                                        prr_gsp(nvec, i),   &
                                        prs_gsp(nvec, i),   &
                                        prg_gsp(nvec, i),   &
                                        t_so(nvec, i)
         if(iostatus/=0) then ! to avoid end of file error.
            exit
         endif
      end do

      ! close unit
      close(20)

      ! -------------------
      ! further required preprocessing
      ! -------------------
      do i = 1, nvec
         iswr    (i,:)  = iswr   (i,:) + iswr_d (i,:) ! incoming short-wave radiation ( direct + diffuse)
         iswr_d  (i,:)  = 0.0_wp
         prr_con (i,:)  = prr_gsp(i,:)                ! rain rate
         prs_con (i,:)  = prs_gsp(i,:) + prg_gsp(i,:) ! snow rate ( snow + graupel )
         prr_gsp (i,:)  = 0.0_wp
         prs_gsp (i,:)  = 0.0_wp
         prg_gsp (i,:)  = 0.0_wp

         ! Check if precipitation is partitioned for each time step
         do j = 1, nsteps
            if (tot_prec(i,j) > 0.0_wp .AND. prr_con(i,j) == 0.0_wp .AND. prs_con(i,j) == 0.0_wp ) then
               ! If total precipitation is provided, but the phase is undetermined, use a simple air temperature threshold
               if ( t(i,j) > 275.15 ) then
                  prr_con(i,j) = tot_prec(i,j)
               else
                  prs_con(i,j) = tot_prec(i,j)
               end if
            end if
            tot_prec(i,j) = 0.0_wp
         end do

      end do

! ------------------------------------------------------------------------------
! - end subroutine read_input
! ------------------------------------------------------------------------------

   end subroutine read_input


! ------------------------------------------------------------------------------
! + begin subroutine read_state
! ------------------------------------------------------------------------------
   subroutine read_state()
! ------------------------------------------------------------------------------
! description:
!   this routine reads initial state for the snow cover scheme
! ------------------------------------------------------------------------------

      integer   :: &
         i, j, ksn,            &  ! loop index
         nNodes

      integer, parameter :: &
         ivstart = 1, ivend = 1

      character(*), parameter :: state_file = "./inp/icon_15min_2021.state"

      character(len=256) :: line = "", dummy = ""

      IF (nvec .ne. 1) THEN

         write(0,*) "ERROR: computation domain size (nvec) != 1! This is not yet supported for itype_nix_start > 1"

      ENDIF

      open(unit=20, status = "old", file=state_file)

      ! Read first line
      read(20, '(A)') line
      read(line, '(A8,I6)') dummy, ke_snow

      ! Read header row
      read(20,*)

      call allocate_fields()

      ! Read layer data
      DO i = ivstart, ivend ! initial horizontal loop

         DO ksn = 1, ke_snow

            read(20, *) j, dzm_sn(i, ksn), t_sn(i, ksn), theta_a(i, ksn), theta_w(i, ksn), theta_i(i, ksn)

         ENDDO

      ENDDO

      ! Read header rows
      read(20, '(A)') line
      read(line, '(A7,I6)') dummy, nNodes
      read(20,*)

      IF (nNodes .NE. ke_snow+1) THEN

         write(0,*) nNodes, ke_snow, line
         write(0,*) "ERROR: number of nodes != number of snow layers + 1 in snow state file!"
         CALL EXIT(1)

      ENDIF


      ! Read nodal data
      DO i = ivstart, ivend ! initial horizontal loop

         DO ksn = 1, ke_snow+1

            read(20, *) j, t_sn_n(i, ksn)

        ENDDO

      ENDDO

      close(20)

! ------------------------------------------------------------------------------
! - end subroutine read_state
! ------------------------------------------------------------------------------

   end subroutine read_state

! =============================================================================
! - end module for input
! ==============================================================================


end module input

