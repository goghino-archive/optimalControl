#!/bin/bash
#compilation script to compile ipopt matlab interface
#archimedes
#first call regular make to build object files and then call this script to link everything and create mex file
MKLROOT=/opt/intel/compilers_and_libraries_2017/linux/mkl
SCHUR_LIB=/home/kardos/block_solve

g++ -pthread -Wl,--no-undefined  -shared -O3 -pipe -DNDEBUG -pedantic-errors -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long -fopenmp -fPIC -fno-common -fexceptions -DIPOPT_BUILD -DMATLAB_MEXFILE -static-libgcc -static-libstdc++ -Wl,--version-script,"/opt/MATLAB/R2014b/extern/lib/glnxa64/mexFunction.map" matlabexception.o matlabfunctionhandle.o matlabjournal.o iterate.o ipoptoptions.o options.o sparsematrix.o callbackfunctions.o matlabinfo.o matlabprogram.o ipopt.o\
    -lipopt  -ldl  -lpthread  -lm  -lgomp  -lm  -ldl -L/home/kardos/Ipopt-3.12.4/build_matlab/lib \
    -Wl,--start-group /home/kardos/lib/pardiso/libpardiso.a /home/kardos/lib/pardiso/libpils_pardiso.a /home/kardos/lib/pardiso/libmetis41_pardiso.a /home/kardos/lib/pardiso/libmetis41-P_pardiso.a -Wl,--end-group \
    -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group \
    ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a  \
    -lgfortran \
    -L${SCHUR_LIB}/lib -lschur \
    -L/home/kardos/openmpi/lib -lmpi -lopen-rte  -lopen-pal  -lrt  -libverbs  -lnuma  -lutil  /usr/lib/gcc/x86_64-linux-gnu/5/libquadmath.a \
    -Wl,-rpath-link,/opt/MATLAB/R2014b/bin/glnxa64 -L"/opt/MATLAB/R2014b/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ \
    -o /home/kardos/Ipopt-3.12.4/build_matlab/Ipopt/contrib/MatlabInterface/src/ipopt.mexa64

