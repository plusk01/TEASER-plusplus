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
#include <chrono>
#include <iostream>

#include "mex.h"
#include <Eigen/Core>

#include <teaser/graph.h>

#include "mexutils.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // nlhs   number of expected outputs
  // plhs   array to be populated by outputs (data passed back to matlab)
  // nrhs   number of inputs
  // prhs   array poplulated by inputs (data passed from matlab)

  if (nrhs != 1 && nrhs != 2) {
    mexErrMsgIdAndTxt("pmc:nargin", "One or two arguments (data, params) required.");
  }

  teaser::MaxCliqueSolver::Params params;
  if (nrhs >= 2) {
    if (mxIsStruct(prhs[1])) {
      const mxArray * field = nullptr;
      if (field = mxGetField(prhs[1], 0, "num_threads"))
        params.num_threads = static_cast<int>(*mxGetPr(field));
      if (field = mxGetField(prhs[1], 0, "solver_mode"))
        params.solver_mode = static_cast<teaser::MaxCliqueSolver::CLIQUE_SOLVER_MODE>(*mxGetPr(field));
    } else {
      mexErrMsgIdAndTxt("pmc:params", "Second argument must be a struct.");
    }
  }

  std::chrono::high_resolution_clock::time_point t1, t2;
  teaser::MaxCliqueSolver clique_solver(params);
  std::vector<int> maxclique;

  if (mxIsChar(prhs[0])) {
    // data is a mtx filename
    std::string tmp(mxArrayToString(prhs[0]));
    std::string filename = tmp;

    t1 = std::chrono::high_resolution_clock::now();
    maxclique = clique_solver.findMaxClique(filename);
    t2 = std::chrono::high_resolution_clock::now();

  } else if (mxIsSparse(prhs[0])) {
    mexErrMsgIdAndTxt("pmc:params", "First arg as sparse not yet implemented.");
  } else if (mxIsNumeric(prhs[0])) {

    Eigen::MatrixXd M;
    mexMatrixToEigen(prhs[0], &M);

    const int n = M.cols();

    teaser::Graph G;
    G.populateVertices(n);
    for (size_t i=0; i<n; ++i) {
      for (size_t j=i+1; j<n; ++j) {
        if (M(i,j) > 0) {
          G.addEdge(i, j);
        }
      }
    }

    t1 = std::chrono::high_resolution_clock::now();
    maxclique = clique_solver.findMaxClique(G);
    t2 = std::chrono::high_resolution_clock::now();

  } else {
    mexErrMsgIdAndTxt("pmc:params", "First arg may only be mtxfilename string or sparse / full numeric matrix.");
  }

  const double time_e2e = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) / 1e6;

  std::sort(maxclique.begin(), maxclique.end());
  const teaser::MaxCliqueSolver::Info pmcinfo = clique_solver.getSolutionInfo();

  // transform maxclique indices into doubles for matlab (with 1-based indexing)
  std::vector<double> mcidx; mcidx.reserve(maxclique.size());
  std::transform(maxclique.begin(), maxclique.end(), std::back_inserter(mcidx),
          [](int x){ return x+1; }); // for matlab 1-based indexing

  //
  // PMC solution info
  //

  if (nlhs > 1) {
    mexErrMsgIdAndTxt("pmc:nargout", "Only one output supported (pmcinfo)");
  }

  // create struct
  static constexpr int FIELDS = 10;
  const char * fieldnames[] = {"inliers", "t_e2e", "t_heu", "omega_heu", "t_exact", "omega_exact", "exact_ran", "num_vertices", "num_edges", "density"};
  plhs[0] = mxCreateStructMatrix(1,1,FIELDS,fieldnames);

  // struct data
  mxArray * inliers = mxCreateDoubleMatrix(1, mcidx.size(), mxREAL);
  memcpy(mxGetData(inliers), &mcidx[0], mcidx.size()*sizeof(double));
  mxArray * t_e2e = mxCreateDoubleScalar(time_e2e);
  mxArray * t_heu = mxCreateDoubleScalar(pmcinfo.t_heu);
  mxArray * omega_heu = mxCreateDoubleScalar(pmcinfo.omega_heu);
  mxArray * t_exact = mxCreateDoubleScalar(pmcinfo.t_exact);
  mxArray * omega_exact = mxCreateDoubleScalar(pmcinfo.omega_exact);
  mxArray * pmc_exact = mxCreateDoubleScalar(((pmcinfo.exact_ran)?1.0:0.0));
  mxArray * num_vertices = mxCreateDoubleScalar(pmcinfo.num_vertices);
  mxArray * num_edges = mxCreateDoubleScalar(pmcinfo.num_edges);
  mxArray * density = mxCreateDoubleScalar(pmcinfo.density);

  // set struct data
  mxSetFieldByNumber(plhs[0],0,0,inliers);
  mxSetFieldByNumber(plhs[0],0,1,t_e2e);
  mxSetFieldByNumber(plhs[0],0,2,t_heu);
  mxSetFieldByNumber(plhs[0],0,3,omega_heu);
  mxSetFieldByNumber(plhs[0],0,4,t_exact);
  mxSetFieldByNumber(plhs[0],0,5,omega_exact);
  mxSetFieldByNumber(plhs[0],0,6,pmc_exact);
  mxSetFieldByNumber(plhs[0],0,7,num_vertices);
  mxSetFieldByNumber(plhs[0],0,8,num_edges);
  mxSetFieldByNumber(plhs[0],0,9,density);
}