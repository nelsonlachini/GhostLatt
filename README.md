## For now compatible with:
gcc-8
clang-10

## For a quick test without use of autotools and make: 
gcc-8 exPlaquette.c src/utilities.c src/algebra.c src/thermalization_hb.c src/measurement.c src/inverters.c  src/gaugefix.c src/instantontools.c src/nrutil.c src/statistics.c src/hoek_custom.c src/fox.c src/center.c -o temp -lm

## To use autotools, follow the usual steps:
- run bootstrap.sh 
- create a build directory 
- run ../configure
- then run make

A minimal set of installation commands would be:
- sh bootstrap.sh
- mkdir build; cd build
- ../configure
- make

## Use deploy.sh and clean.sh to quickly compile from zero
