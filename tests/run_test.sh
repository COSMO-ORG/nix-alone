resolutions="900 20"                                    # List of temporal resolutions in seconds
numberoflayer="200 100 10 14"


runNIX() {
    #
    # Run NIX
    #

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
    pushd ../

    mv ./src/input.f90 ./src/input.f90.bak
    sed 's/.\/inp\/icon_15min_2021.inp/'${stnfile}'/g' ./src/input.f90.bak > ./src/input.f90
    
    mv ./src/mo_nix_config.f90 ./src/mo_nix_config.f90.bak
    sed -e 's/zdt               = 900.0_wp/zdt      = '${ts}'_wp/' -e 's/ke_snow = 10/ke_snow   = '${nlayer}'/' -e 's/min_height_layer  = 0.01_wp/min_height_layer = '${layer_min}'_wp/' -e 's/max_height_layer  = 0.05_wp/max_height_layer = '${layer_max}'_wp/' ./src/mo_nix_config.f90.bak > ./src/mo_nix_config.f90

    make

    mv src/input.f90.bak src/input.f90
    mv src/mo_nix_config.f90.bak src/mo_nix_config.f90

    popd

    ../nix > ${stnfile}.${nlayer}layers.out
    export TZ=UTC; awk -v dt=${ts} -v t=2021-10-01T00:00 -F, 'BEGIN {data=0; d=mktime(sprintf("%04d %02d %02d %02d %02d %02d 0", substr(t,1,4), substr(t,6,2), substr(t,9,2), substr(t,12,2), substr(t,15,2), substr(t,19,2)))} {if(!data) {print} else {if(/^0500/) {print "0500," strftime("%Y-%m-%dT%H:%M:%S", d+dt*$2)} else {print}}; if(/\[DATA\]/) {data=1}}' output.pro > ${stnfile}.${nlayer}.pro
    export TZ=UTC; awk -v dt=${ts} -v t=2021-10-01T00:00 'BEGIN {data=0; d=mktime(sprintf("%04d %02d %02d %02d %02d %02d 0", substr(t,1,4), substr(t,6,2), substr(t,9,2), substr(t,12,2), substr(t,15,2), substr(t,19,2)))} {if(!data) {if(/^fields/) {gsub(/timestep/, "TIMESTAMP", $0)}; print} else {printf("%s", strftime("%Y-%m-%dT%H:%M:%S", d+dt*$1)); for(i=2; i<=NF; i++) {printf " %s", $i}; printf "\n"}; if(/\[DATA\]/) {data=1}}' output.smet > ${stnfile}.${nlayer}.smet
    rm output.pro output.smet
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
    echo "PSUM::accumulate::period = ${ts}" >> io_snowpack_${ts}.ini
    # Run SNOWPACK:
    ../../snowpack/bin/snowpack -c io_snowpack_${ts}.ini -b ${begin} -e ${end}
}


for ts in ${resolutions}
do
    stnfile="WFJ_forcing_${ts}s.txt"
    for nl in ${numberoflayer}
    do
        runNIX
    done
    runSNOWPACK
done
