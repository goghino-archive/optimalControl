#!/bin/bash
MKLROOT=/opt/intel/compilers_and_libraries_2017/linux/mkl
PARDISO_DIR=/home/kardos/lib/pardiso

../configure F77=gfortran CC=gcc CXX=g++ \
    --disable-dependency-linking --disable-shared --with-pic --enable-static --without-ma27 \
    ADD_CFLAGS="-fPIC -fno-common -fexceptions"       \
    ADD_CXXFLAGS="-fopenmp -fPIC -fno-common -fexceptions"     \
    ADD_FFLAGS="-fPIC -fexceptions -fbackslash"                       \
     --enable-inexact-solver \
     --with-matlab-home=/opt/MATLAB/R2014b --enable-matlab-static \
     --with-blas-lib="-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -ldl -lpthread -lm"  \
    --with-lapack-lib="${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a" \
    --with-pardiso="-Wl,--start-group ${PARDISO_DIR}/libpardiso.a  ${PARDISO_DIR}/libpils_pardiso.a ${PARDISO_DIR}/libmetis41_pardiso.a  ${PARDISO_DIR}/libmetis41-P_pardiso.a -Wl,--end-group  -ldl -lpthread -lm -lgomp -lgfortran"



#     --with-blas-lib="$-L{MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lgomp" \
#     --with-lapack-lib="-L${MKROOT}/lib/intel64 -lmkl_lapack95_lp64" \
#     --with-pardiso="-L${PARDISO_DIR} -lpardiso500-GNU481-X86-64 -lgfortran"
