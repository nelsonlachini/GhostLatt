#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([LGT_C], [0.1], [nelsonlachini@gmail.com])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([src/utilities.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AC_LANG([C])
CC="gcc-8"
CFLAGS="$AM_CFLAGS $CFLAGS -O2"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"
LIBS="$LIBS -lm"

AC_MSG_CHECKING([that a minimal LGT_C program compiles]);
AC_LINK_IFELSE(
[AC_LANG_SOURCE([[
    #include <stdlib.h>
    #include <complex.h>
    #include <stdio.h>

    // global variables (to do: change that to local)
    const int N = 4;

    const int Nt = N;
    const int totalV = N*N*N*Nt;
    const int spatialV = N*N*N;
    const int dimLattice = totalV*4*4;
    const int colorV = totalV*3;
    const int N2 = N*N;


    int main(){

        //LOCAL PARAMETERS
        int N_therm = 1000;		// #thermalization sweeps
        int N_hb 	= 1;		// #heat-bath steps per sweep
        int N_mic 	= N/2;	        // #overrelaxation steps per sweep
        
        double initial_order = 0;			//0 is cold start

        return(0);
    }	
]])],
[AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])]
[AC_MSG_ERROR([Could not compile a minimal LGT_C program])])

AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

