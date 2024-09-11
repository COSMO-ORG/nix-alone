resolutions="3600 900 20"                                    # List of temporal resolutions in seconds
numberoflayer="100 10 3"


runNIX() {
    #
    # Run NIX
    #

    if (( ${nl} == 3 )); then
        nlayer=3
        layer_min=0.005
        layer_max=0.6
    fi
    if (( ${nl} == 10 )); then
        nlayer=10
        layer_min=0.002
        layer_max=0.3
    fi
    if (( ${nl} == 14 )); then
        nlayer=14
        layer_min=0.002
        layer_max=0.3
    fi
    if (( ${nl} == 100 )); then
        nlayer=100
        layer_min=0.0002
        layer_max=0.03
    fi
    if (( ${nl} == 200 )); then
        nlayer=200
        layer_min=0.0002
        layer_max=0.03
    fi
    if (( ${nl} == 1000 )); then
        nlayer=1000
        layer_min=0.0002
        layer_max=0.03
    fi

    # Make sure the output frequency is the same in real time independent of chosen model timestep
    output_freq_tgt=3600
    # Note: if output frequency is higher than model time step, the output_freq is forced to 1 (i.e., every timestep)
    output_freq=$(echo ${ts} | awk -v sec_between=${output_freq_tgt} '{printf "%d", (sec_between<$1)?(1):(sec_between/$1)}')

    pushd ../

    # Set the input file
    mv ./src/input.f90 ./src/input.f90.bak
    sed 's/.\/inp\/nix.inp/'${stnfile}'/g' ./src/input.f90.bak > ./src/input.f90

    # Make required changes to the configuration
    mv ./src/mo_nix_config.f90 ./src/mo_nix_config.f90.bak
    sed -e 's/zdt               = 900.0_wp/zdt      = '${ts}'_wp/' \
        -e 's/ke_snow = 10/ke_snow   = '${nlayer}'/' \
        -e 's/min_height_layer  = 0.01_wp/min_height_layer = '${layer_min}'_wp/' \
        -e 's/max_height_layer  = 0.05_wp/max_height_layer = '${layer_max}'_wp/' \
        -e 's/_output_freq = 12/_output_freq = '${output_freq}'/' \
        ./src/mo_nix_config.f90.bak > ./src/mo_nix_config.f90

    # Compile NIX
    make

    mv src/input.f90.bak src/input.f90
    mv src/mo_nix_config.f90.bak src/mo_nix_config.f90

    popd

    /usr/bin/time -a -o 'timings.txt' -f "NIX ${nlayer}L ${ts}s : %e" ../nix > ${stnfile}.${nlayer}layers.out
    mv output.smet ${stnfile}.${nlayer}.smet
    mv output.pro ${stnfile}.${nlayer}.pro

}

runSNOWPACK() {
    #
    # Run SNOWPACK
    #

    metfile=$(awk 'BEGIN {input=0} {if(input) {if(/METEOPATH/) {path=$NF}; if(/STATION1/) {file=$NF}}; if(/\[Input\]/) {input=1}} END {print path "/" file}' io_snowpack_base.ini)
    begin=$(awk '{if(/\[DATA\]/) {getline; print $1; exit}}' ${metfile})
    end=$(tail -n 1 ${metfile} | awk '{print $1}')
    # Take care of specific settings to run at the requested time resolution:
    echo "IMPORT_BEFORE = io_snowpack_base.ini" > io_snowpack_${ts}.ini
    echo "[Input]" >> io_snowpack_${ts}.ini
    echo "STATION1	= WFJ2_${ts}s.smet" >> io_snowpack_${ts}.ini
    echo "[Output]" >> io_snowpack_${ts}.ini
    echo "EXPERIMENT = output_${ts}s" >> io_snowpack_${ts}.ini
    echo "[Snowpack]" >> io_snowpack_${ts}.ini
    echo ${ts} | awk '{printf "CALCULATION_STEP_LENGTH = %f\n", $1/60.}' >> io_snowpack_${ts}.ini
    echo "[Filters]" >> io_snowpack_${ts}.ini
    echo "PSUM::arg1::cst         = ${ts}" >> io_snowpack_${ts}.ini
    echo "[Interpolations1D]" >> io_snowpack_${ts}.ini
    echo "PSUM::arg1::period = ${ts}" >> io_snowpack_${ts}.ini
    # Run SNOWPACK:
    /usr/bin/time -a -o 'timings.txt' -f "SNOWPACK ${ts}s : %e" ../../snowpack/bin/snowpack -c io_snowpack_${ts}.ini -b ${begin} -e ${end}
}


for ts in ${resolutions}
do
    stnfile="WFJ_forcing_${ts}s.txt"
    timespan=($(awk '{ts2=$1; if(/\[DATA\]/) {getline; ts1=$1}} END {print ts1, ts2}' SNOWPACK/WFJ2_${ts}s.smet))
    for nl in ${numberoflayer}
    do
        runNIX
    done
    runSNOWPACK
done
