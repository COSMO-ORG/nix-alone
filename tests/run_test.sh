resolutions="900 20"                                    # List of temporal resolutions in seconds


runNIX() {
    #
    # Run NIX
    #

    nsteps=$(wc -l ${stnfile} | cut -d' ' -f1)

    pushd ../

    mv ./src/input.f90 ./src/input.f90.bak
    sed 's/.\/inp\/WFJ_forcing.txt/'${stnfile}'/g' ./src/input.f90.bak > ./src/input.f90
    
    mv ./src/config.f90 ./src/config.f90.bak
    sed -e 's/zdt      = 900.0_wp/zdt      = '${ts}'_wp/g' -e 's/nsteps    = 1270080/nsteps    = '${nsteps}'/' ./src/config.f90.bak > ./src/config.f90

    make

    mv src/input.f90.bak src/input.f90
    mv src/config.f90.bak src/config.f90

    popd

    ../nix > ${stnfile}.out
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
    runNIX
    runSNOWPACK
done
