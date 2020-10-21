/**
 * @brief      MEX wrapper for TEASER's PMC functionality
 * @author     Parker Lusk <plusk@mit.edu>
 * @date       21 Oct 2020
 */

/**
 * Copyright 2020, Massachusetts Institute of Technology,
 * Cambridge, MA 02139
 * All Rights Reserved
 * Authors: Jingnan Shi, et al. (see THANKS for the full author list)
 * See LICENSE for the license information
 */

#include <algorithm>
#include <iostream>
#include <chrono>

#include "mex.h"
#include <Eigen/Core>

#include "teaser/graph.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // nlhs   number of expected outputs
  // plhs   array to be populated by outputs (data passed back to matlab)
  // nrhs   number of inputs
  // prhs   array poplulated by inputs (data passed from matlab)

  if (nrhs != 2) {
    mexErrMsgIdAndTxt("pmc:nargin", "Two arguments (type, data) required.");
  }

  teaser::MaxCliqueSolver::Params params;
  teaser::MaxCliqueSolver clique_solver(params);
  std::vector<int> max_clique;

  const auto t1 = std::chrono::high_resolution_clock::now();
  const double type = *mxGetPr(prhs[0]);
  if (type == 0) {
    // data is a mtx filename
    std::string tmp(mxArrayToString(prhs[1]));
    std::string filename = tmp;

    max_clique = clique_solver.findMaxClique(filename);

  } else if (type == 1) {

  }

  std::sort(max_clique.begin(), max_clique.end());
  const auto t2 = std::chrono::high_resolution_clock::now();
  const double t_e2e = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) / 1e6;

  const teaser::MaxCliqueSolver::Info pmcinfo = clique_solver.getSolutionInfo();


  //
  // PMC solution info
  //

  if (nlhs > 1) {
    mexErrMsgIdAndTxt("pmc:nargout", "Only one output supported (pmcinfo)");
  }

  // create struct
  static constexpr int FIELDS = 8;
  const char * fieldnames[] = {"t_heu", "omega_heu", "t_exact", "omega_exact", "exact_ran", "num_vertices", "num_edges", "density"};
  plhs[0] = mxCreateStructMatrix(1,1,FIELDS,fieldnames);

  // struct data
  mxArray * t_heu = mxCreateDoubleScalar(pmcinfo.t_heu);
  mxArray * omega_heu = mxCreateDoubleScalar(pmcinfo.omega_heu);
  mxArray * t_exact = mxCreateDoubleScalar(pmcinfo.t_exact);
  mxArray * omega_exact = mxCreateDoubleScalar(pmcinfo.omega_exact);
  mxArray * pmc_exact = mxCreateDoubleScalar(((pmcinfo.exact_ran)?1.0:0.0));
  mxArray * num_vertices = mxCreateDoubleScalar(pmcinfo.num_vertices);
  mxArray * num_edges = mxCreateDoubleScalar(pmcinfo.num_edges);
  mxArray * density = mxCreateDoubleScalar(pmcinfo.density);

  // set struct data
  mxSetFieldByNumber(plhs[0],0,0,t_heu);
  mxSetFieldByNumber(plhs[0],0,1,omega_heu);
  mxSetFieldByNumber(plhs[0],0,2,t_exact);
  mxSetFieldByNumber(plhs[0],0,3,omega_exact);
  mxSetFieldByNumber(plhs[0],0,4,pmc_exact);
  mxSetFieldByNumber(plhs[0],0,5,num_vertices);
  mxSetFieldByNumber(plhs[0],0,6,num_edges);
  mxSetFieldByNumber(plhs[0],0,7,density);
}