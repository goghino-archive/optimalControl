#!/bin/bash
#MKLROOT=/opt/intel/mkl/lib/intel64 #old version
#MKLROOT=/opt/intel/compilers_and_libraries_2016.2.181/linux/mkl
MKLROOT=/opt/intel/compilers_and_libraries_2017/linux/mkl
MKLIB=-L${MKLROOT}
PARDISOLIB=/home/kardos/lib/pardiso

#MPI_CC=/home/kardos/openmpi-2.0.0/bin/mpicc
#MPI_CXX=/home/kardos/openmpi-2.0.0/bin/mpic++

MPI_CC=gcc
MPI_CXX=g++

../configure -disable-shared ADD_CFLAGS="-O3 -fPIC -fexceptions -fopenmp" ADD_CXXFLAGS="-O3 -fPIC -fexceptions -fopenmp" F77=gfortran CC=${MPI_CC} CXX=${MPI_CXX} --with-blas-lib="${MKLIB} -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lgomp" --with-lapack-lib="${MKLIB} -lmkl_lapack95_lp64" --with-pardiso="-L${PARDISOLIB} -lpardiso500-GNU481-X86-64 -lgfortran"

#https://projects.coin-or.org/Ipopt/wiki/MatlabInterface
