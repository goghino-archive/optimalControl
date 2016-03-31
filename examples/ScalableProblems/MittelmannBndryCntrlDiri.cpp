// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: MittelmannBndryCntrlDiri.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Andreas Waechter             IBM    2005-10-18
//                  based on MyNLP.cpp

#include "MittelmannBndryCntrlDiri.hpp"

#ifdef HAVE_CASSERT
# include <cassert>
#else
# ifdef HAVE_ASSERT_H
#  include <assert.h>
# else
#  error "don't have header file for assert"
# endif
#endif

using namespace Ipopt;

/* Constructor. */
MittelmannBndryCntrlDiriBase::MittelmannBndryCntrlDiriBase()
    :
    y_d_(NULL)
{}

MittelmannBndryCntrlDiriBase::~MittelmannBndryCntrlDiriBase()
{
  delete [] y_d_;
}

void
MittelmannBndryCntrlDiriBase::SetBaseParameters(Index NS, Index N, Number alpha, Number lb_y,
    Number ub_y, Number lb_u, Number ub_u,
    Number d_const)
{
  NS_ = NS; // number of scenarios for stochastic control 
  N_ = N;   // number of meshpoints (excluding boundary)
  h_ = 1./(N+1); // step size
  hh_ = h_*h_;
  lb_y_ = lb_y;  // bound on state variable
  ub_y_ = ub_y;
  lb_u_ = lb_u;  // bound on control variable
  ub_u_ = ub_u;
  d_const_ = d_const; // rhs of equality (laplace) constraint
  alpha_ = alpha;

  // Initialize the target state  variables
  delete [] y_d_;
  y_d_ = new Number[(N_)*(N_)*NS_]; //TODO include also control U?
  for (Index k=0; k<NS_; k++) {
    for (Index i=0; i< N_; i++) {
      for (Index j=0; j< N_; j++) {
        y_d_[y_index(i,j,k)] = y_d_cont(x1_grid(i+1),x2_grid(j+1)); //TODO: add noise, sum error also on boundary to dirichlet cond.
      }
    }
  }
}

bool MittelmannBndryCntrlDiriBase::get_nlp_info(
  Index& n, Index& m, Index& nnz_jac_g,
  Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // We for each of the N_+2 times N_+2 mesh points we have the value
  // of the functions y, including the control parameters on the boundary
  // # of mesh points including boundary 
  n = N_*N_*NS_ + (4*N_ + 4);

  // For each of the N_ times N_ interior mesh points we have the
  // discretized PDE.
  // # of equality constraints
  m = N_*N_*NS_;

  // y(i,j), y(i-1,j), y(i+1,j), y(i,j-1), y(i,j+1) for each
  // of the N_*N_ discretized PDEs
  nnz_jac_g = 5*N_*N_*NS_;

  // diagonal entry for each y(i,j) in the interior
  nnz_h_lag = N_*N_*NS_;
  if (alpha_>0.) {
    // and one entry for u(i,j) in the bundary if alpha is not zero
    nnz_h_lag += 4*N_ + 4;
  }

  // We use the C indexing style for row/col entries (corresponding to
  // the C notation, starting at 0)
  index_style = C_STYLE;

  return true;
}

bool
MittelmannBndryCntrlDiriBase::get_bounds_info(Index n, Number* x_l, Number* x_u,
    Index m, Number* g_l, Number* g_u)
{
  //[Y1, Y2, ..., Yns, U] of size [NN, NN, ..., NN, 4N+4]

  // Set overall bounds on the state variables y
  for (Index k=0; k<NS_; k++) {
    for (Index i=0; i<N_; i++) {
      for (Index j=0; j<N_; j++) {
        Index iy = y_index(i,j,k);
        x_l[iy] = lb_y_;
        x_u[iy] = ub_y_;
      }
    }
  }

  Index offset = NS_ * N_*N_;
  // Set the overall bounds on the control variables u
    for (Index i=0; i<4*N_; i++) {
        x_l[offset+i] = lb_u_;
        x_u[offset+i] = ub_u_;
    }
  
  // The values of y on the corners doesn't appear anywhere, so we fix
  // them to zero
  offset += 4*N_;
  x_l[offset]   = x_u[offset]   = 0.;
  x_l[offset+1] = x_u[offset+1] = 0.;
  x_l[offset+2] = x_u[offset+2] = 0.;
  x_l[offset+3] = x_u[offset+3] = 0.;

  // all discretized PDE constraints have right hand side equal to
  // minus the constant value of the function d  
  // m is already N*N*NS_
  for (Index i=0; i<m; i++)
  {
    g_l[i] = -hh_*d_const_;
    g_u[i] = -hh_*d_const_;
  }

  return true;
}

bool
MittelmannBndryCntrlDiriBase::get_starting_point(Index n, bool init_x, Number* x,
    bool init_z, Number* z_L, Number* z_U,
    Index m, bool init_lambda,
    Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // set all y's to the perfect match with y_d
  for (Index k=0; k<NS_; k++) {
      for (Index i=0; i< N_; i++) { 
          for (Index j=0; j<N_; j++) {
              x[y_index(i,j,k)] = y_d_[y_index(i,j,k)]; //TODO: size of the x = new [size???], where is it initialized?
          }
      }
  }

  // set initial control on the boundary from the interior of the constraints
  Number umid = (ub_u_ + lb_u_)/2.;
  Index offset = NS_ * N_*N_;
  for (Index i=0; i<4*N_ + 4; i++) {
      x[offset + i] = umid;
  }

  return true;
}

bool
MittelmannBndryCntrlDiriBase::get_scaling_parameters(Number& obj_scaling,
    bool& use_x_scaling, Index n, Number* x_scaling,
    bool& use_g_scaling, Index m, Number* g_scaling)
{
  obj_scaling = 1./hh_;
  use_x_scaling = false;
  use_g_scaling = false;
  return true;
}

bool
MittelmannBndryCntrlDiriBase::eval_f(Index n, const Number* x,
                                     bool new_x, Number& obj_value)
{
  // return the value of the objective function
  obj_value = 0.;

  for (Index k=0; k<NS_; k++) {
      Number obj_value_i = 0.;
      
      // First the integration of y-td over the interior
      for (Index i=0; i<N_; i++) {
          for (Index j=0; j< N_; j++) {
              Index iy = y_index(i,j,k);
              Number tmp = x[iy] - y_d_[iy];
              obj_value_i += tmp*tmp;
          }
      }
      obj_value += obj_value_i * hh_/2.;
    }

    // Now the integration of u over the boundary
    Index offset = NS_* N_ * N_;
    if (alpha_>0.) {
        Number usum = 0.;
        for (Index i=0; i<4*N_+4; i++) {
            usum += x[offset+i]*x[offset+i]; //TODO - single control for all scenarios
        }
        obj_value += alpha_*h_/2.*usum;
    }
  
  //obj_value = obj_value / NS_; // if scaled the emphasis is put on the lin. constraints

  return true;
}

bool
MittelmannBndryCntrlDiriBase::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    // return the gradient of the objective function grad_{x} f(x)

    for (Index k=0; k<NS_; k++) {
        // now let's take care of the nonzero values coming from the
        // integrand over the interior
        for (Index i=0; i<N_; i++) {
            for (Index j=0; j< N_; j++) {
                Index iy = y_index(i,j,k);
                grad_f[iy] = hh_*(x[iy] - y_d_[iy]);
            }
        }
    }

    Index offset = NS_*N_*N_;
    // The values for variables on the boundary
    if (alpha_ > 0.) {
        for (Index i=0; i<4*N_+4; i++) {
            grad_f[offset+i] = alpha_*h_*x[offset+i];
        }
    }
    else {
        for (Index i=0; i< 4*N_+4; i++) {
            grad_f[offset+i] = 0.;
        }
    }

  return true;
}

bool MittelmannBndryCntrlDiriBase::eval_g(Index n, const Number* x, bool new_x,
    Index m, Number* g)
{
  // return the value of the constraints: g(x)

  // compute the discretized PDE for each interior grid point
  Index ig = 0;

  for (Index k=0; k<NS_; k++) {
      for (Index i=0; i<N_; i++) {
          for (Index j=0; j<N_; j++) {
              Number val;

              // Start with the discretized Laplacian operator
              val = 4.* x[y_index(i,j,k)]
                  - x[y_index(i-1,j,k)] - x[y_index(i+1,j,k)] 
                  - x[y_index(i,j-1,k)] - x[y_index(i,j+1,k)]; //!!! if i,j<0 || i,j>=N_

              g[ig] = val;
              ig++;
          }
      }
  }

  DBG_ASSERT(ig==m);

  return true;
}

bool MittelmannBndryCntrlDiriBase::eval_jac_g(Index n, const Number* x, bool new_x,
    Index m, Index nele_jac, Index* iRow, Index *jCol,
    Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian of the constraints, laplace for scenarios

      Index ijac = 0;
      Index ig = 0;

      for (Index k=0; k<NS_; k++) {
          for (Index i=0; i< N_; i++) {
              for (Index j=0; j< N_; j++) {

                  // y(i,j)
                  iRow[ijac] = ig;
                  jCol[ijac] = y_index(i,j,k);
                  ijac++;

                  // y(i-1,j)
                  iRow[ijac] = ig;
                  jCol[ijac] = y_index(i-1,j,k);
                  ijac++;

                  // y(i+1,j)
                  iRow[ijac] = ig;
                  jCol[ijac] = y_index(i+1,j,k);
                  ijac++;

                  // y(i,j-1)
                  iRow[ijac] = ig;
                  jCol[ijac] = y_index(i,j-1,k);
                  ijac++;

                  // y(i,j+1)
                  iRow[ijac] = ig;
                  jCol[ijac] = y_index(i,j+1,k);
                  ijac++;

                  ig++;
              }
          }
      }

    DBG_ASSERT(ijac==nele_jac);
  }
  else {
    // return the values of the jacobian of the constraints
    Index ijac = 0;

    for (Index k=0; k<NS_; k++) {
        for (Index i=0; i< N_; i++) {
            for (Index j=0; j< N_; j++) {
                // y(i,j)
                values[ijac] = 4.;
                ijac++;

                // y(i-1,j)
                values[ijac] = -1.;
                ijac++;

                // y(i+1,j)
                values[ijac] = -1.;
                ijac++;

                // y(1,j-1)
                values[ijac] = -1.;
                ijac++;

                // y(1,j+1)
                values[ijac] = -1.;
                ijac++;
            }
        }
    }

    DBG_ASSERT(ijac==nele_jac);
  }

  return true;
}

bool
MittelmannBndryCntrlDiriBase::eval_h(Index n, const Number* x, bool new_x,
                                     Number obj_factor, Index m,
                                     const Number* lambda,
                                     bool new_lambda, Index nele_hess, Index* iRow,
                                     Index* jCol, Number* values)
{
    if (values == NULL) {
        // return the structure. This is a symmetric matrix, fill the lower left
        // triangle only.

        Index ihes=0;

        for (Index k=0; k<NS_; k++) {
            // First the diagonal entries for y(i,j)
            for (Index i=0; i< N_; i++) {
                for (Index j=0; j< N_; j++) {
                    iRow[ihes] = y_index(i,j,k);
                    jCol[ihes] = y_index(i,j,k);
                    ihes++;
                }
            }
        }

        Index offset = NS_ * N_ * N_;
        if (alpha_>0.) {
            // Now the diagonal entries for u at the boundary
            for (Index i=0; i<4*N_+4; i++) {
                iRow[ihes] = offset + i;
                jCol[ihes] = offset + i;
                ihes++;
            }
        }

        DBG_ASSERT(ihes==nele_hess);
    }
    else {
        // return the values

        Index ihes=0;

        for (Index k=0; k<NS_; k++) {
            // First the diagonal entries for y(i,j)
            for (Index i=0; i< N_; i++) {
                for (Index j=0; j< N_; j++) {
                    // Contribution from the objective function
                    values[ihes] = obj_factor*hh_;

                    ihes++;
                }
            }
        }

        // Now the diagonal entries for u(i,j)
        if (alpha_>0.) {
            // Now the diagonal entries for u at the boundary
            for (Index i=0; i<4*N_+4; i++) {
                values[ihes] = obj_factor*h_*alpha_;
                ihes++;
            }
        }

        DBG_ASSERT(ihes==nele_hess);
    }

  return true;
}

void
MittelmannBndryCntrlDiriBase::finalize_solution(SolverReturn status,
    Index n, const Number* x, const Number* z_L, const Number* z_U,
    Index m, const Number* g, const Number* lambda, Number obj_value,
    const IpoptData* ip_data,
    IpoptCalculatedQuantities* ip_cq)
{
    FILE* fp1 = fopen("solution1.txt", "w+"); 
    FILE* fp2 = fopen("solution2.txt", "w+");
    FILE* fp3 = fopen("solution3.txt", "w+");
    FILE* fp4 = fopen("solution4.txt", "w+");


    Index offset = NS_*N_*N_ ;
    for (Index k=0; k<NS_; k++) {
        for (Index i=0; i<N_; i++) {
            fprintf(fp1, "%15.8e %15.8e\n", i*h_ + 1.*k, x[offset + i]);      // North
            fprintf(fp2, "%15.8e %15.8e\n", i*h_ + 1.*k, x[offset+N_ + i]);   // South
            fprintf(fp3, "%15.8e %15.8e\n", i*h_ + 1.*k, x[offset+2*N_ + i]); // East
            fprintf(fp4, "%15.8e %15.8e\n", i*h_ + 1.*k, x[offset+3*N_ + i]); // West
            //fprintf(fp, "y[%6d,%6d] = %15.8e\n", i, j, x[y_index(i,j,k)]);
        }
    }

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
}