# nix-alone

# compile nix standalone - tsa

module load cmake
mkdir build && cd build
cmake ..
make

# if one needs to configure
ccmake ..
# now configure by pressing 'c' currently no configure option availble
# one might have to press c twice then generate with 'g' or quit
make

# run nix from run_folder - executable linked to here
cd run_folder
./nix

# currently path to input data input data is hardcoded in /src/input.90
# all other possible setting can be done in /src/config.f90




