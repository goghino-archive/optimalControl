/*********************************************************************
 * Demo.cpp
 *
 * This file shows the basics of setting up a mex file to work with
 * Matlab.  This example shows how to use 2D matricies.  This may
 * 
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 * For more information, see my site: www.shawnlankton.com
 * by: Shawn Lankton
 *
 ********************************************************************/
#include <mex.h>   
#include <iostream>

#include <mpi.h>

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex; 

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

inline void mpi_check(int mpi_call)
{
    if ((mpi_call) != 0) { 
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:MPI",
                "MPI error detected\n.");
        return;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //check number of arguments
    if (nrhs != 0)
    {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                 "No inputs allowed.");
         return;
    }

    if (nlhs != 1) {
         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                 "One output required.");
         return;
    }


    double ierr;
    ierr = MPI_Init(NULL, NULL);
    mpi_check(ierr);

    //associate outputs
    plhs[0] = mxCreateDoubleScalar(ierr);

    mexPrintf("MPI_INIT\n");

    return;
}
