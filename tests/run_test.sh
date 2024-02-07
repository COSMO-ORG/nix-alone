resolutions="900 20"                                    # List of temporal resolutions in seconds


runSimulation() {
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


for ts in ${resolutions}
do
    stnfile="WFJ_forcing_${ts}s.txt"
    runSimulation
done
