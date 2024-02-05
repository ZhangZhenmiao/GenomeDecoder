#pragma once
#include <string>
#include <unordered_map>
#include "dot_graph.hpp"

// merge non-branching paths
void merge_non_branching_paths(std::unordered_map<std::string, Node>& graph, int k);

// remove simple bulges, allow arbitrary multiplicity of legs, can replace simple_bulge_removal
void multi_bulge_removal(std::unordered_map<std::string, Node>& graph, int k, unsigned& removed_bulge);

// remove general whirls, may result in inaccurate edge lengths, can replace multi_whirl_removal
void general_whirl_removal(std::unordered_map<std::string, Node>& graph, int k, unsigned& removed_whirls);

// remove bulges with a single edge and a path of at most x edges
void resolving_bulge_with_single_edge_and_x_edge_path(std::unordered_map<std::string, Node>& graph, int k, unsigned& removed_paths, int x = 4, double identity = 0.9);

// remove bulges with two paths of multiple edges (at most x)
void resolving_bulge_with_two_multi_edge_paths(std::unordered_map<std::string, Node>& graph, int k, unsigned& removed_paths, int x = 4, double identity = 0.9);

// remove bulges with two paths of multiple edges (at most x)
void resolving_bulge_with_two_multi_edge_paths_rerouting(std::unordered_map<std::string, Node>& graph, int k, unsigned& removed_paths, int x = 4, double identity = 0.9);