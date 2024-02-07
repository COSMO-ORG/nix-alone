#
# Settings by the user
#
metfile=icon_15min_2021.inp
target_resolutions="900 20"                                    # List of temporal downscaling resolutions in seconds
meteoio_timeseries_path=~/src/usr/bin/

#
# Additional settings, derived from files
#
begin=$(awk '(NR==2) {print $1 "T" $2}' ${metfile})
end=$(awk 'END {print $1 "T" $2}' ${metfile})
stn_id=$(fgrep CSV1_ID *ini | awk '{gsub(/ /, "", $NF); print $NF}')


#
# Execute MeteoIO's meteoio_timeseries to interpolate temporally for each target_resolution
#
for target_resolution in ${target_resolutions}
do
	echo "IMPORT_BEFORE = ./io_base.ini" > io.ini
	echo "[Interpolations1D]" >> io.ini
	echo "TOT_PREC::accumulate::period    = ${target_resolution}" >> io.ini

	# Execute meteoio_timeseries
	${meteoio_timeseries_path}/meteoio_timeseries -c io.ini -b ${begin} -e ${end} -s $(echo ${target_resolution} | awk '{print $1/60}')

	# Convert output file from meteoio_timseries to inp file for nix
	awk -f convert_smet_to_inp.awk ${stn_id}.smet > WFJ_forcing_${target_resolution}s.txt
	rm ${stn_id}.smet io.ini
done
