AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I .buildutils/m4
AM_CXXFLAGS = -fPIC -Wall -Wextra
bin_PROGRAMS = plaquette gluonPropagator
plaquette_SOURCES = examples/plaquette.c src/utilities.c src/algebra.c src/thermalization_hb.c src/measurement.c src/inverters.c
gluonPropagator_SOURCES = examples/gluonPropagator.c src/utilities.c src/algebra.c src/thermalization_hb.c src/measurement.c src/inverters.c src/gaugefix.c src/statistics.c src/nrutil.c
