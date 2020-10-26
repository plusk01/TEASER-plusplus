/**
 * Copyright 2020, Massachusetts Institute of Technology,
 * Cambridge, MA 02139
 * All Rights Reserved
 * Authors: Jingnan Shi, et al. (see THANKS for the full author list)
 * See LICENSE for the license information
 */

#include <iostream>
#include <chrono>

#include "teaser/graph.h"
#include "pmc/pmc.h"

vector<int> teaser::MaxCliqueSolver::findMaxClique(teaser::Graph graph) {

  // Create a PMC graph from the TEASER graph
  vector<int> edges;
  vector<long long> vertices;
  vertices.push_back(0);

  const std::vector<int> all_vertices = graph.getVertices();
  for (const auto& i : all_vertices) {
    const std::vector<int>& c_edges = graph.getEdges(i);
    edges.insert(edges.end(), c_edges.begin(), c_edges.end());
    vertices.push_back(edges.size());
  }

  // Use PMC to calculate
  pmc::pmc_graph G(vertices, edges);

  return callPMC(G);
}

vector<int> teaser::MaxCliqueSolver::findMaxClique(const std::string& mtxfilename) {

  // Use PMC to read matrix market file
  const auto old_buffer = std::cout.rdbuf(nullptr);
  pmc::pmc_graph G(mtxfilename);
  std::cout.rdbuf(old_buffer);

  if (G.num_vertices() <= 0 || G.num_edges() == 0) return {};

  return callPMC(G);
}

std::vector<int> teaser::MaxCliqueSolver::callPMC(pmc::pmc_graph& G) {
  // Handle deprecated field
  if (!params_.solve_exactly) {
    params_.solver_mode = CLIQUE_SOLVER_MODE::PMC_HEU;
  }

  // Prepare PMC input
  // TODO: Incorporate this to the constructor
  pmc::input in;
  in.algorithm = 0;
  in.threads = params_.num_threads;
  in.experiment = 0;
  in.lb = 0;
  in.ub = 0;
  in.param_ub = 0;
  in.adj_limit = 20000;
  in.time_limit = params_.time_limit;
  in.remove_time = 4;
  in.graph_stats = false;
  in.verbose = false;
  in.help = false;
  in.MCE = false;
  in.decreasing_order = false;
  in.heu_strat = "kcore";
  in.vertex_search_order = "deg";

  TEASER_DEBUG_INFO_MSG("Using " << in.threads << " threads");
  omp_set_num_threads(in.threads); // keep in sync with PMC

  // vector to represent max clique
  vector<int> C;

  // upper-bound of max clique
  auto t1 = std::chrono::high_resolution_clock::now();
  G.compute_cores();
  auto max_core = G.get_max_core();
  auto t2 = std::chrono::high_resolution_clock::now();
  const double t_kcore = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) / 1e6;

  TEASER_DEBUG_INFO_MSG("Max core number: " << max_core);
  TEASER_DEBUG_INFO_MSG("Num vertices: " << G.num_vertices());

  // capture some info about the input graph
  solninfo_.num_vertices = G.num_vertices();
  solninfo_.num_edges = G.num_edges();
  solninfo_.density = G.density();

  // check for k-core heuristic threshold
  // check whether threshold equals 1 to short circuit the comparison
  if (params_.solver_mode == CLIQUE_SOLVER_MODE::KCORE_HEU &&
      params_.kcore_heuristic_threshold != 1 &&
      max_core > static_cast<int>(params_.kcore_heuristic_threshold *
                                  static_cast<double>(G.num_vertices()))) {
    TEASER_DEBUG_INFO_MSG("Using K-core heuristic finder.");
    // remove all nodes with core number less than max core number
    // k_cores is a vector saving the core number of each vertex
    auto k_cores = G.get_kcores();
    for (int i = 1; i < k_cores->size(); ++i) {
      // Note: k_core has size equals to num vertices + 1
      if ((*k_cores)[i] >= max_core) {
        C.push_back(i-1);
      }
    }
    return C;
  }

  if (in.ub == 0) {
    in.ub = max_core + 1;
  }

  // lower-bound of max clique
  if (in.lb == 0 && in.heu_strat != "0") { // skip if given as input
    const auto old_buffer = std::cout.rdbuf(nullptr);
    t1 = std::chrono::high_resolution_clock::now();
    pmc::pmc_heu maxclique(G, in);
    in.lb = maxclique.search(G, C);
    t2 = std::chrono::high_resolution_clock::now();
    solninfo_.t_heu = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) / 1e6;
    solninfo_.t_heu += t_kcore;
    std::cout.rdbuf(old_buffer);

    solninfo_.omega_heu = C.size();
  }

  assert(in.lb != 0);
  if (in.lb == 0) {
    // This means that max clique has a size of one
    TEASER_DEBUG_ERROR_MSG("Max clique lower bound equals to zero. Abort.");
    return C;
  }

  if (in.lb == in.ub) {
    return C;
  }

  // Optional exact max clique finding
  if (params_.solver_mode == CLIQUE_SOLVER_MODE::PMC_EXACT) {
    // The following methods are used:
    // 1. k-core pruning
    // 2. neigh-core pruning/ordering
    // 3. dynamic coloring bounds/sort
    // see the original PMC paper and implementation for details:
    // R. A. Rossi, D. F. Gleich, and A. H. Gebremedhin, “Parallel Maximum Clique Algorithms with
    // Applications to Network Analysis,” SIAM J. Sci. Comput., vol. 37, no. 5, pp. C589–C616, Jan.
    // 2015.
    const auto old_buffer = std::cout.rdbuf(nullptr);
    t1 = std::chrono::high_resolution_clock::now();
    if (G.num_vertices() < in.adj_limit) {
      G.create_adj();
      pmc::pmcx_maxclique finder(G, in);
      finder.search_dense(G, C);
    } else {
      pmc::pmcx_maxclique finder(G, in);
      finder.search(G, C);
    }
    t2 = std::chrono::high_resolution_clock::now();
    solninfo_.t_exact = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) / 1e6;
    solninfo_.t_exact += t_kcore;
    std::cout.rdbuf(old_buffer);

    solninfo_.exact_ran = true;
    solninfo_.omega_exact = C.size();
  }

  return C;
}