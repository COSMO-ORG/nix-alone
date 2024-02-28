MODULE output
! ------------------------------------------------------------------------------
!
! description: this module handles the output of the snow cover scheme
!
! ------------------------------------------------------------------------------
   USE fields                     ! contains all required global fields
   USE mo_kind,                    ONLY: wp

! ------------------------------------------------------------------------------
! DECLARATIONS
! ------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

!------------------------------------------------------------------------------
! Anything public?
!------------------------------------------------------------------------------

   PUBLIC :: write_output

CONTAINS

! =============================================================================
! + Begin subroutine: write_output
! ============================================================================

   SUBROUTINE write_output(nvec, ivstart, ivend, n, nout, &
   &                 top, ke_snow, dzm_sn, rho_sn       , &
   &                 theta_i, theta_w, theta_a          , &
   &                 t_sn, t_sn_n, hn_sn, t , alpha_sn)


      ! Subroutine Arguments
      INTEGER, INTENT(IN)   :: &
         nvec          , & ! < array dimensions
         ivstart       , & ! < start index for computations in the parallel program
         ivend         , & ! < end index for computations in the parallel program
         n             , & ! < time step
         nout              ! < output time step (save every nout time step)

      INTEGER, DIMENSION(nvec), INTENT(INOUT) :: &
         top                 ! index of the first (top) layer index       (-)

      INTEGER, INTENT(IN)  ::  &
         ke_snow             ! number of snow layers

      REAL    (KIND = wp), DIMENSION(nvec,1:ke_snow),  INTENT(INOUT) :: &
         dzm_sn        , &   ! layer thickness between main levels        (m)
         rho_sn        , &   ! snow layer density                         (kg/m3)
         theta_i       , &   ! volumetric ice content                     (-)
         theta_w       , &   ! water ice content                          (-)
         theta_a       , &   ! air ice content                            (-)
         t_sn                ! snow temperature (main level)              (K)

      REAL    (KIND = wp), DIMENSION(nvec,1:ke_snow),  INTENT(INOUT) :: &
         t_sn_n

      REAL    (KIND = wp), DIMENSION(nvec),  INTENT(INOUT) :: &
         hn_sn       , &        ! new snow amount                            (m)
         alpha_sn


      REAL    (KIND = wp), DIMENSION(nvec),  INTENT(IN) :: &
         t                ! air temperature first level                (K)

      ! Local variables

      INTEGER ::           &

         i            , &   ! loop index in x-direction
         ksn                ! loop index in z-direction (snow layers)

      REAL    (KIND = wp), DIMENSION(1:ke_snow) :: &
         tmp_sn             ! temporary utility array used for copying

      REAL    (KIND = wp), DIMENSION(1:ke_snow+1) :: &
         tmp_sn_n             ! temporary utility array used for copying
      ! ------------------------------------------------------------------------------
      ! Begin subroutine
      ! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------

!       ! ----------------------
!       ! profile information
!       ! ----------------------

!       ! open file
   DO i=ivstart,ivend
      IF (n .EQ. 1) THEN
         open(unit=22, file='output.pro')
         write(22,'(A)') "[STATION_PARAMETERS]"
         write(22,'(A)') "StationName= NIX"
         write(22,'(A)') "Latitude= -999"
         write(22,'(A)') "Longitude= -999"
         write(22,'(A)') "Altitude= -999"
         write(22,'(A)') "SlopeAngle= 0.00"
         write(22,'(A)') "SlopeAzi= 0.00"
         write(22,'(A)') ""
         write(22,'(A)') "[HEADER]"
         write(22,'(A)') "0500,Date"
         write(22,'(A)') "0501,nElems,height [> 0: top, < 0: bottom of elem.] (cm)"
         write(22,'(A)') "0502,nElems,element density (kg m-3)"
         write(22,'(A)') "0503,nElems,element temperature (degC)"
         write(22,'(A)') "0506,nElems,liquid water content by volume (%)"
         write(22,'(A)') "0515,nElems,ice volume fraction (%)"
         write(22,'(A)') "0516,nElems,air volume fraction (%)"
         write(22,'(A)') ""
         write(22,'(A)') "[DATA]"
                        close(22)
      ENDIF
      IF (MOD(n, nout) == 0) then
         open(unit=22, file='output.pro', status='old', action='write', access='sequential', form='formatted', position='append')
         IF (top(i) .EQ. 0) THEN
            ! Special case with no snow elements
            write(22, '(A,I0)') '0500,', n
            write(22, '(A)') '0501,1,0'
         ELSE
            ! Write time step
            write(22, '(A,I0)', advance='yes') '0500,', n
            ! Layer spacing info
            write(22, '(A,I0)', advance='no') '0501,', top(i)
            DO ksn = 1, top(i), 1
                write(22, '(A,F0.3)', advance='no') ',', 100.*hm_sn(i,ksn)	! Convert to cm
            END DO
            write(22, '(A)') ""
            ! Bulk density
            write(22, '(A,I0)', advance='no') '0502,', top(i)
            DO ksn = 1, top(i), 1
                write(22, '(A,F0.3)', advance='no') ',', rho_sn(i,ksn)
            END DO
            write(22, '(A)') ""
            ! Temperature
            write(22, '(A,I0)', advance='no') '0503,', top(i)
            DO ksn = 1, top(i), 1
                write(22, '(A,F0.3)', advance='no') ',', t_sn(i,ksn) - 273.15	! Convert to degC
            END DO
            write(22, '(A)') ""
            ! Liquid water content
            write(22, '(A,I0)', advance='no') '0506,', top(i)
            DO ksn = 1, top(i), 1
                write(22, '(A,F0.3)', advance='no') ',', 100.*theta_w(i,ksn)	! Convert to %
            END DO
            write(22, '(A)') ""
            ! Ice content
            write(22, '(A,I0)', advance='no') '0515,', top(i)
            DO ksn = 1, top(i), 1
                write(22, '(A,F0.3)', advance='no') ',', 100.*theta_i(i,ksn)	! Convert to %
            END DO
            write(22, '(A)') ""
            ! Air content
            write(22, '(A,I0)', advance='no') '0516,', top(i)
            DO ksn = 1, top(i), 1
                write(22, '(A,F0.3)', advance='no') ',', 100.*theta_a(i,ksn)	! Convert to %
            END DO
            write(22, '(A)') ""
         ENDIf
         CLOSE(22)
      ENDIF
   ENDDO

! ! ------------------------------------------------------------------------------
! ! - end subroutine write_output
! ! ------------------------------------------------------------------------------

   END SUBROUTINE write_output

! =============================================================================
! - end module output
! ==============================================================================

END MODULE output


