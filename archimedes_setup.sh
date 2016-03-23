#!/bin/bash
MKLROOT=/opt/intel/mkl/lib/intel64
MKLIB=-L${MKLROOT}
PARDISOLIB=/home/drosos/Libraries/linuxAMD64

../configure CXXFLAGS="-O2 -fPIC -fopenmp -m64 -fstack-protector-all" F77=gfortran CC=mpicc CXX=mpic++   --with-blas-lib="${MKLIB} -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lgomp" --with-lapack-lib="${MKLIB} -lmkl_lapack95_lp64" --with-pardiso="-L${PARDISOLIB} -lpardiso500-GNU481-X86-64 -lgfortran"
