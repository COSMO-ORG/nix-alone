cmake_minimum_required(VERSION 3.14)

project(nix LANGUAGES Fortran)

# we default to Release build type
if(NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE "Release")
endif()

# set all flag options
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
 set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall")
 set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback")
 set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -xHOST -unroll-aggressive -qopt-prefetch")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
 set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall")
 set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan")
 set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -march=native")
endif()

file(GLOB_RECURSE sources src/*.f90)

add_executable(nix ${sources})
