MODULE mo_nix_output
! ------------------------------------------------------------------------------
!
! description: this module handles the output of the snow cover scheme
!
! ------------------------------------------------------------------------------
   USE fields                              ! contains all required global fields
   USE mo_kind,       ONLY : wp
   USE mo_nix_config, ONLY : pro_output_file, smet_output_file, &
                             pro_output_freq, smet_output_freq

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

   SUBROUTINE write_output(nvec, ivstart, ivend, n, pro_output_freq, &
   &                 top, ke_snow, dzm_sn, rho_sn                  , &
   &                 theta_i, theta_w, theta_a                     , &
   &                 t_sn, t_sn_n, hn_sn, t , alpha_sn)


      ! Subroutine Arguments
      INTEGER, INTENT(IN)   :: &
         nvec          , & ! < array dimensions
         ivstart       , & ! < start index for computations in the parallel program
         ivend         , & ! < end index for computations in the parallel program
         n             , & ! < time step
         pro_output_freq   ! < output time step (save every nout time step)

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

!       NOTE: in order to be able to visualize the *pro file with niViz.org, it is necessary to include timestamps. Use the following oneliner to achieve this (insert model timestep dt and initial timestamp t)
!       export TZ=UTC; awk -v dt=900 -v t=2020-10-01T00:00 -F, 'BEGIN {data=0; d=mktime(sprintf("%04d %02d %02d %02d %02d %02d 0", substr(t,1,4), substr(t,6,2), substr(t,9,2), substr(t,12,2), substr(t,15,2), substr(t,19,2)))} {if(!data) {print} else {if(/^0500/) {print "0500," strftime("%Y-%m-%dT%H:%M:%S", d+dt*$2)} else {print}}; if(/\[DATA\]/) {data=1}}' output.pro > output2.pro

   IF (pro_output_freq .gt. 0) THEN
      DO i=ivstart,ivend
         IF (n .EQ. 1) THEN
            open(unit=22, file=pro_output_file)
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
            write(22,'(A)') "0521,nElems,thermal conductivity (W K-1 m-1)"
            write(22,'(A)') "0522,nElems,absorbed shortwave radiation (W m-2)"
            write(22,'(A)') ""
            write(22,'(A)') "[DATA]"
                           close(22)
         ENDIF
         IF (MOD(n, pro_output_freq) == 0) then
            open(unit=22, file=pro_output_file, status='old', action='write', access='sequential', &
                 form='formatted', position='append')
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
               ! Thermal conductivity
               write(22, '(A,I0)', advance='no') '0521,', top(i)
               DO ksn = 1, top(i), 1
                  write(22, '(A,F0.3)', advance='no') ',', hcon_sn(i,ksn)
               END DO
               write(22, '(A)') ""
               ! Absorbed shortwave radiation
               write(22, '(A,I0)', advance='no') '0522,', top(i)
               DO ksn = 1, top(i), 1
                  write(22, '(A,F0.3)', advance='no') ',', swflx_sn_abs(i,ksn)
               END DO
               write(22, '(A)') ""
            ENDIf
            CLOSE(22)
         ENDIF
      ENDDO
   ENDIF

!       ! ----------------------
!       ! timeseries information
!       ! ----------------------

!       NOTE: in order to be able to visualize the *pro file with niViz.org, it is necessary to include timestamps. Use the following oneliner to achieve this (insert model timestep dt and initial timestamp t)
!       export TZ=UTC; awk -v dt=900 -v t=2020-10-01T00:00 'BEGIN {data=0; d=mktime(sprintf("%04d %02d %02d %02d %02d %02d 0", substr(t,1,4), substr(t,6,2), substr(t,9,2), substr(t,12,2), substr(t,15,2), substr(t,19,2)))} {if(!data) {if(/^fields/) {gsub(/timestep/, "TIMESTAMP", $0)}; print} else {printf("%s", strftime("%Y-%m-%dT%H:%M:%S", d+dt*$1)); for(i=2; i<=NF; i++) {printf " %s", $i}; printf "\n"}; if(/\[DATA\]/) {data=1}}' output.smet > output2.smet

   IF (smet_output_freq .gt. 0) THEN
      DO i=ivstart,ivend
         IF (n .EQ. 1) THEN
            open(unit=22, file=smet_output_file)
            write(22,'(A)') "SMET 1.1 ASCII"
            write(22,'(A)') "[HEADER]"
            write(22,'(A)') "station_id   = NIX"
            write(22,'(A)') "station_name = NIX"
            write(22,'(A)') "latitude     = -999"
            write(22,'(A)') "longitude    = -999"
            write(22,'(A)') "altitude     = -999"
            write(22,'(A)') "nodata       = -999"
            write(22,'(A)') "fields       = timestep TA QI VW P PSUM PSUM_PH HS &
                                             SWE TSS HN MS_SNOWPACK_RUNOFF SWR_NET &
                                             ISWR RSWR ILWR OLWR ALBEDO SHF LHF EBAL"
            write(22,'(A)') "[DATA]"
            close(22)
         ENDIF
         IF (MOD(n, smet_output_freq) == 0) then
            open(unit=22, file=smet_output_file, status='old', action='write', access='sequential', &
                 form='formatted', position='append')
               write(22, '(I0)', advance="no") n
               ! TA
               write(22, '(AF0.3)', advance="no") " ", t(n)
               ! QI
               write(22, '(AF0.6)', advance="no") " ", qv(i,n)
               ! VW
               write(22, '(AF0.3)', advance="no") " ", sqrt(u(i,n)*u(i,n)+v(i,n)*v(i,n))
               ! P
               write(22, '(AF0.3)', advance="no") " ", ps(i,n)
               ! PSUM
               write(22, '(AF0.3)', advance="no") " ", prr_con(i,n)+prs_con(i,n)+prr_gsp(i,n)+prs_gsp(i,n)+prg_gsp(i,n)
               ! PSUM_PH
               IF (prr_con(i,n)+prs_con(i,n)+prr_gsp(i,n)+prs_gsp(i,n)+prg_gsp(i,n) .gt. 0) THEN
                  write(22, '(AF0.3)', advance="no") " ", (prr_con(i,n)+prr_gsp(i,n)) / &
                                                              (prr_con(i,n)+prs_con(i,n)+prr_gsp(i,n)+prs_gsp(i,n)+prg_gsp(i,n))
               ELSE
                  write(22, '(A)', advance="no") " -999"
               ENDIF
               ! HS
               write(22, '(AF0.3)', advance="no") " ", h_snow(i)
               ! SWE
               write(22, '(AF0.3)', advance="no") " ", swe_sn(i)
               ! Snow surface temperature
               write(22, '(AF0.3)', advance="no") " ", t_sn_sfc(i)
               ! New snow amount
               write(22, '(AF0.3)', advance="no") " ", hn_sn(i)
               ! Runoff
               write(22, '(AF0.3)', advance="no") " ", runoff_sn(i)
               ! Net shortwave
               write(22, '(AF0.3)', advance="no") " ", swflx_sn_net(i)
               ! Incoming shortwave
               write(22, '(AF0.3)', advance="no") " ", iswr(i,n)
               ! Reflected shortwave
               write(22, '(AF0.3)', advance="no") " ", swflx_sn_up(i)
               ! Incoming longwave
               write(22, '(AF0.3)', advance="no") " ", ilwr(i,n)
               ! Outgoing longwave
               write(22, '(AF0.3)', advance="no") " ", lwflx_sn_up(i)
               ! Albedo
               write(22, '(AF0.3)', advance="no") " ", alpha_sn(i)
               ! Sensible heat flux
               write(22, '(AF0.3)', advance="no") " ", shflx_sn(i)
               ! Latent heat flux
               write(22, '(AF0.3)', advance="no") " ", lhflx_sn(i)
               ! Latent heat flux
               write(22, '(AF0.3)', advance="no") " ", for_sn(i)
               write(22, '(A)') ""
            CLOSE(22)
         ENDIF
      ENDDO
   ENDIF
! ! ------------------------------------------------------------------------------
! ! - end subroutine write_output
! ! ------------------------------------------------------------------------------

   END SUBROUTINE write_output

! =============================================================================
! - end module output
! ==============================================================================

END MODULE mo_nix_output


