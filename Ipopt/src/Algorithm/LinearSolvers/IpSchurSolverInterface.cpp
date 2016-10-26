// Copyright (C) 2005, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: IpSchurSolverInterface.cpp 2594 2015-08-09 14:31:05Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-17
//
//           Olaf Schenk                      Univ of Basel 2005-09-20
//                  - changed options, added PHASE_ flag

/* some useful links:
 * MKL documentation: https://software.intel.com/en-us/intel-mkl/documentation
 * API differences MKL vs Basel PARDISO: http://software.intel.com/en-us/articles/summary-of-api-differences-between-intel-mkl-pardiso-and-university-of-basel-pardiso-400
 */

#include "IpoptConfig.h"
#include "IpSchurSolverInterface.hpp"
#include <mpi.h>
#include <math.h>


#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

#ifdef HAVE_CSTDLIB
# include <cstdlib>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# else
#  error "don't have header file for stdlib"
# endif
#endif

#ifdef HAVE_CSTRING
# include <cstring>
#else
# ifdef HAVE_STRING_H
#  include <string.h>
# else
#  error "don't have header file for string"
# endif
#endif

// determine the correct name of the Pardiso function
#if defined(_MSC_VER) && defined(HAVE_PARDISO)
# define PARDISOINIT_FUNC PARDISOINIT
# define PARDISO_FUNC PARDISO
#else
# define PARDISOINIT_FUNC F77_FUNC(pardisoinit,PARDISOINIT)
# define PARDISO_FUNC F77_FUNC(pardiso,PARDISO)
#endif


/* Prototypes for Pardiso's subroutines */
extern "C"
{
#if defined(HAVE_PARDISO_OLDINTERFACE) || defined(HAVE_PARDISO_MKL)
  void PARDISOINIT_FUNC(void* PT, const ipfint* MTYPE, ipfint* IPARM);
#else
  // The following is a fix to allow linking with Pardiso library under Windows
  void PARDISOINIT_FUNC(void* PT, const ipfint* MTYPE,
                        const ipfint* SOLVER,
                        ipfint* IPARM,
                        double* DPARM,
                        ipfint* ERROR);
#endif
  void PARDISO_FUNC(void** PT, const ipfint* MAXFCT,
                    const ipfint* MNUM, const ipfint* MTYPE,
                    const ipfint* PHASE, const ipfint* N,
                    const double* A, const ipfint* IA,
                    const ipfint* JA, const ipfint* PERM,
                    const ipfint* NRHS, ipfint* IPARM,
                    const ipfint* MSGLVL, double* B, double* X,
                    ipfint* ERROR, double* DPARM);


#ifdef PARDISO_MATCHING_PREPROCESS
  void smat_reordering_pardiso_wsmp_(const ipfint* N, const ipfint* ia, const ipfint* ja, const double* a_, ipfint* a2, ipfint* ja2,  double* a2_,
                                     ipfint* perm2,  double* scale2, ipfint* tmp2_, ipfint preprocess );
#endif

}

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  SchurSolverInterface::SchurSolverInterface()
      :
      a_(NULL),
      negevals_(-1),
      initialized_(false),
      schurSolver(-2, 1, MPI_COMM_WORLD) //equivalent to commented code below
  {
    //const int pardiso_mtype = -2; // symmetric H_i
    //const int schur_factorization = 1;
    //schurSolver = SchurSolve(pardiso_mtype, schur_factorization, MPI_COMM_WORLD);
    DBG_START_METH("SchurSolverInterface::SchurSolverInterface()",dbg_verbosity);
  }

  SchurSolverInterface::~SchurSolverInterface()
  {
    DBG_START_METH("SchurSolverInterface::~SchurSolverInterface()",
                   dbg_verbosity);

    // Tell Pardiso to release all memory

    delete[] a_;
  }

  void SchurSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddIntegerOption(
        "problem_dimension",
        "Dimension of the grid for the single scenario",
        0,
        "The solver requires dimension of the grid N");

    roptions->AddIntegerOption(
        "problem_scenarios",
        "Number of scenarios for stochastic optimization.",
        0,
        "The solver requires the number of the scenarios");
  }

  bool SchurSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetIntegerValue("problem_dimension", N_, prefix);
    options.GetIntegerValue("problem_scenarios", NS_, prefix);

    printf("SchurSolverPGInterface::InitializeImpl N %d, NS %d\n", N_, NS_); 

  
    // Reset all private data
    dim_=0;
    nonzeros_=0;
    initialized_=false;
    delete[] a_;
    a_ = NULL;

    return true;
  }

  ESymSolverStatus SchurSolverInterface::MultiSolve(bool new_matrix,
      const Index* ia,
      const Index* ja,
      Index nrhs,
      double* rhs_vals,
      bool check_NegEVals,
      Index numberOfNegEVals)
  {
    DBG_START_METH("SchurSolverInterface::MultiSolve",dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    DBG_ASSERT(initialized_);

    CSRdouble* KKT = new CSRdouble((int)dim_, (int)dim_, (int)nonzeros_, ia, ja, a_);
    KKT->matrixType = SYMMETRIC;

#ifdef DEBUG
    char buffer[200];
    sprintf (buffer, "~/Ipopt-3.12.4/build_mpi/Ipopt/examples/ScalableProblems/PardisoMat_%d_%d.csr", N_, NS_);
    KKT->writeToFile(buffer);
#endif
    
    static int initialized = 0;
    if(!initialized)
    {
        //printf("-------------Initializing Schur Solver at master-----------\n");
        schurSolver.initSystem_OptimalControl(KKT, N_, NS_);
        initialized = 1;
        new_matrix = 0; // new matrix is set by calling initSystem, no need to call update
    }

    //broadcast flag that the solution is not yet found
    int terminate = 0;
    MPI_Bcast(&terminate, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //broadcast flag new_matrix to child processes
    MPI_Bcast(&new_matrix, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    if (new_matrix) 
    {
        //printf("--------------Updating system at master-----------\n");
        schurSolver.updateSystem(a_);
    }
    else
    {
        //printf("----------------REUSING MATRIX---------------------\n");
    }

    // Only master contains RHS with actual data at this point
    // it is communicated to children inside the solve
    double* X = new double [dim_ * nrhs];
    schurSolver.solveSystem(X, rhs_vals, nrhs);
    //schurSolver.errorReport(nrhs, *KKT, rhs_vals, X);
    schurSolver.timingReport();  

    // overwrite rhs by the solution
    for (int i=0; i<dim_*nrhs ;i++)
    {
      rhs_vals[i] = X[i];
    }

    // clean up
    delete [] X;
    KKT->clear();
    delete KKT;

    return SYMSOLVER_SUCCESS;  
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus SchurSolverInterface::InitializeStructure
  (Index dim, Index nonzeros,
   const Index* ia,
   const Index* ja)
  {
    DBG_START_METH("SchurSolverInterface::InitializeStructure",dbg_verbosity);
    dim_ = dim;
    nonzeros_ = nonzeros;

    // Make space for storing the matrix elements
    delete[] a_;
    a_ = NULL;
    a_ = new double[nonzeros_];

    initialized_ = true;

    return SYMSOLVER_SUCCESS;
  }

  double* SchurSolverInterface::GetValuesArrayPtr()
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(a_);
    return a_;
  }


  Index SchurSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("SchurSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(negevals_>=0);
    return negevals_;
  }

  bool SchurSolverInterface::IncreaseQuality()
  {
    // At the moment, I don't see how we could tell Pardiso to do better
    // (maybe switch from IPARM[20]=1 to IPARM[20]=2?)
    return false;
  }

} // namespace Ipopt
