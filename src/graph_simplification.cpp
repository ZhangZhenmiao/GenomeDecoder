#include <string>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <queue>
#include <set>
#include "dot_graph.hpp"
#include "edlib.h"
#include "graph_simplification.hpp"

double Graph::percent_identity(std::string cigar) {
    int matches = 0;
    int mismatches = 0;
    int num = 0;
    for (char c : cigar) {
        if (std::isdigit(c))
            num = num * 10 + (c - '0');
        else {
            if (c == 'M')
                matches += num;
            else if (c == 'X')
                mismatches += num;
            else if (c == 'I')
                mismatches += num;
            else if (c == 'D')
                mismatches += num;
            else
                std::cout << "Unknown character in cigar" << std::endl;
            num = 0;
        }
    }
    return 1.0 * matches / (matches + mismatches);
}

std::string Graph::reverse_complementary_node(std::string node) {
    if (node.at(0) == '-')
        return node.substr(1);
    else
        return '-' + node;
}

bool Graph::check_non_branching(std::string node) {
    // there is no self-loop
    if (this->graph[node].incoming_edges.size() == 1 && this->graph[node].outgoing_edges.size() == 1) {
        bool flag = true;
        Edge edge_in, edge_out;
        unsigned multiplicity1, multiplicity2;
        std::string sequence_in, sequence_out;

        // get color and multiplicity of incoming edge
        for (auto&& i : this->graph[node].incoming_edges) {
            // pass bulges
            if (i.second.size() != 1)
                flag = false;
            edge_in = i.second.at(0);
        }
        // get color and multiplicity of outgoing edge
        for (auto&& j : this->graph[node].outgoing_edges) {
            // pass bulges
            if (j.second.size() != 1)
                flag = false;
            edge_out = j.second.at(0);
        }
        // pass cases of different multiplicities
        if (!edge_in.check_equal_multi(edge_out))
            flag = false;
        // the forward non-branching path does not pass check
        return flag;
    }

    // there is self-loop
    else if (this->graph[node].incoming_edges.size() == 2 && this->graph[node].outgoing_edges.size() == 2 && this->graph[node].incoming_edges.find(node) != this->graph[node].incoming_edges.end()) {
        return false;
        bool flag = true;
        Edge edge_in, edge_out, edge_loop;
        // get color and multiplicity of incoming edge and loop
        for (auto&& i : this->graph[node].incoming_edges) {
            // pass bulges
            if (i.second.size() != 1)
                flag = false;
            if (i.first != node)
                edge_in = i.second.at(0);
            else
                edge_loop = i.second.at(0);
        }
        // get color and multiplicity of outgoing edge
        for (auto&& j : this->graph[node].outgoing_edges) {
            // pass bulges
            if (j.second.size() != 1)
                flag = false;
            if (j.first != node)
                edge_out = j.second.at(0);
        }
        // pass cases of different multiplicities
        if (!edge_in.check_equal_multi(edge_out))
            flag = false;
        if (!edge_in.check_folds_multi(edge_loop))
            flag = false;
        // the forward non-branching path does not pass check
        return flag;
    }
    // the node is not part of a non-branching path
    else
        return false;
}

bool Graph::add_node_to_path(Path& path, std::string node, int bulge_leg) {
    if (path.nodes.empty()) {
        // the node should exist
        if (graph.find(node) == graph.end())
            return false;
        path.nodes.push_back(node);
        return true;
    }
    std::string prev_node = path.nodes.at(path.nodes.size() - 1);

    path.nodes.push_back(node);
    path.bulge_legs.push_back(bulge_leg);

    // the edge should exist
    if (this->graph[prev_node].outgoing_edges.find(node) == this->graph[prev_node].outgoing_edges.end()) {
        return false;
    }

    // a path should not have two colors
    if (!path.update_min_multi(this->graph[prev_node].outgoing_edges[node].at(bulge_leg))) {
        return false;
    }

    // a path should be consecutive in the genome path; this will update the multiplicity of the path (maybe use the min multi as alternative)
    if (!this->find_path_from_genome_paths(path)) {
        return false;
    }

    // update path sequence
    path.length += this->graph[prev_node].outgoing_edges[node].at(bulge_leg).length;
    if (path.sequence.empty())
        path.sequence = this->graph[prev_node].outgoing_edges[node].at(bulge_leg).sequence;
    else
        path.sequence += this->graph[prev_node].outgoing_edges[node].at(bulge_leg).sequence.substr(this->k);

    // for palindromic bulge
    if (prev_node == reverse_complementary_node(node)) {
        path.index2palindromic_seq[path.nodes.size() - 2] = this->graph[prev_node].outgoing_edges[node].at(bulge_leg).sequence;
    }

    return true;
}

void Graph::merge_edges(std::string node) {
    std::string node1, node2;
    Edge edge1, edge2;
    // there is no self-loop
    if (this->graph[node].incoming_edges.size() == 1) {
        // std::cout << "Non-branching (no self-loop): " << node << std::endl;
        // get info from edge 1
        for (auto&& i : graph[node].incoming_edges) {
            assert(i.second.size() == 1);
            node1 = i.first;
            edge1 = i.second.at(0);
        }
        // get info from edge 2
        for (auto&& j : graph[node].outgoing_edges) {
            assert(j.second.size() == 1);
            node2 = j.first;
            edge2 = j.second.at(0);
        }
        Edge new_edge(edge1.start_base, edge1.length + edge2.length, edge1.sequence + edge2.sequence.substr(this->k), edge1.multi11, edge1.multi12, edge1.multi21, edge1.multi22, edge1.multi31, edge1.multi32);
        assert(new_edge.length == new_edge.sequence.size() - this->k);

        this->graph[node1].outgoing_edges.erase(node);
        this->graph[node1].outgoing_edges[node2].push_back(new_edge);
        this->graph[node2].incoming_edges.erase(node);
        this->graph[node2].incoming_edges[node1].push_back(new_edge);
        this->graph.erase(node);
        this->erase_node_from_genome_paths(node);
    }
    // there is self-loop
    else {
        Path new_path;
        for (auto&& i : graph[node].incoming_edges) {
            if (i.first != node)
                node1 = i.first;
        }
        for (auto&& j : graph[node].outgoing_edges) {
            if (j.first != node)
                node2 = j.first;
        }

        assert(add_node_to_path(new_path, node1));
        assert(add_node_to_path(new_path, node));

        bool flag = false;
        for (int i = 0; i < this->graph[node].outgoing_edges[node].at(0).multiplicity / this->graph[node].outgoing_edges[node2].at(0).multiplicity; ++i) {
            if (!add_node_to_path(new_path, node)) {
                flag = true;
                break;
            }
        }
        if (flag || !add_node_to_path(new_path, node2))
            return;

        if (this->graph[node].outgoing_edges[node2].at(0).multiplicity != new_path.multiplicity)
            return;

        Edge new_edge(new_path.sequence.at(k), new_path.length, new_path.sequence, new_path.multi11, new_path.multi12, new_path.multi21, new_path.multi22, new_path.multi31, new_path.multi32);

        this->graph[node1].outgoing_edges.erase(node);
        this->graph[node1].outgoing_edges[node2].push_back(new_edge);
        this->graph[node2].incoming_edges.erase(node);
        this->graph[node2].incoming_edges[node1].push_back(new_edge);
        this->graph.erase(node);
        this->erase_node_from_genome_paths(node);
        std::cout << "Non-branching (self-loop): " << node1 << "->" << node << "->" << node2 << ", multi: " << new_path.multiplicity << std::endl;
    }
}

void Graph::merge_non_branching_paths() {
    std::set<std::string> nodes_to_remove;
    for (auto&& node : this->graph) {
        if (nodes_to_remove.find(node.first) != nodes_to_remove.end() || nodes_to_remove.find(reverse_complementary_node(node.first)) != nodes_to_remove.end())
            continue;
        if (check_non_branching(node.first) && check_non_branching(reverse_complementary_node(node.first)))
            nodes_to_remove.insert(node.first);
    }

    for (auto&& node : nodes_to_remove) {
        merge_edges(node);
        merge_edges(reverse_complementary_node(node));
    }
}

bool Graph::check_reverse_bulge(std::string node1, std::string node2, bool check_overlap) {
    std::vector<Edge>& forward_edges = this->graph[node1].outgoing_edges[node2];
    assert(forward_edges.size() >= 2);

    if (this->graph.find(reverse_complementary_node(node2)) == this->graph.end() || this->graph.find(reverse_complementary_node(node1)) == this->graph.end()) {
        std::cout << "Simple bulge (reverse not exist): " << node1 << ", " << node2 << std::endl;
        return false;
    }
    std::vector<Edge>& reverse_edges = this->graph[reverse_complementary_node(node2)].outgoing_edges[reverse_complementary_node(node1)];

    // check whether reverse nodes form a bulge
    if (reverse_edges.size() < 2) {
        std::cout << "Simple bulge (reverse not a bulge): " << node1 << ", " << node2 << std::endl;
        return false;
    }
    // check whether forward and reverse nodes share nodes
    if (check_overlap) {
        if (node1 == reverse_complementary_node(node2) || node2 == reverse_complementary_node(node1))
            return false;
    }

    return true;
}

std::string Graph::collapse_bulge(std::string node1, std::string node2, unsigned& removed_bulges, double similarity) {
    auto& edges = this->graph[node1].outgoing_edges[node2];
    auto& edges_reverse = graph[node2].incoming_edges[node1];

    while (true) {
        int leg1 = -1, leg2 = -1;
        double sim = 0;
        for (int i = 0; i < edges.size() - 1; ++i) {
            for (int j = i + 1; j < edges.size(); ++j) {
                double sim1 = PI_edlib(edges.at(i).sequence, edges.at(j).sequence);
                double sim2 = PI_edlib(reverse_complementary(edges.at(i).sequence), reverse_complementary(edges.at(j).sequence));
                sim = std::min(sim1, sim2);
                if (sim >= similarity) {
                    leg1 = i;
                    leg2 = j;
                    break;
                }
            }
            if (leg1 != -1 && leg2 != -1)
                break;
        }
        if (leg1 == -1 || leg2 == -1) {
            for (int i = 0; i < edges.size() - 1; ++i) {
                for (int j = i + 1; j < edges.size(); ++j) {
                    double sim1 = PI_edlib(edges.at(i).sequence, edges.at(j).sequence);
                    double sim2 = PI_edlib(reverse_complementary(edges.at(i).sequence), reverse_complementary(edges.at(j).sequence));
                    sim = std::min(sim1, sim2);
                    int color1 = edges.at(i).color, color2 = edges.at(j).color;
                    std::cout << "(Similarity check failed) Simple bulge: " << node1 << "->" << node2 << ", sim (length) " << sim << ", color1 " << color1 << ", color2 " << color2 << std::endl;
                }
            }
            break;
        }

        // keep the longer edge
        int color1 = edges.at(leg1).color, color2 = edges.at(leg2).color;
        if (edges.at(leg2).length > edges.at(leg1).length) {
            edges.at(leg2).add_multi_from_edge_or_path(edges.at(leg1));
            edges.erase(edges.begin() + leg1);

            edges_reverse.at(leg2).add_multi_from_edge_or_path(edges_reverse.at(leg1));
            edges_reverse.erase(edges_reverse.begin() + leg1);
        }
        else {
            edges.at(leg1).add_multi_from_edge_or_path(edges.at(leg2));
            edges.erase(edges.begin() + leg2);

            edges_reverse.at(leg1).add_multi_from_edge_or_path(edges_reverse.at(leg2));
            edges_reverse.erase(edges_reverse.begin() + leg2);
        }
        std::cout << "Simple bulge: " << node1 << "->" << node2 << ", sim (length) " << sim << ", color1 " << color1 << ", color2 " << color2 << std::endl;
        removed_bulges += 1;
    }
    return edges.at(0).sequence;
}

// remove simple bulges, allow arbitrary multiplicity of legs
void Graph::multi_bulge_removal(unsigned& removed_bulges, bool check_overlap, double similarity) {
    removed_bulges = 0;
    for (auto&& node : this->graph) {
        for (auto&& i : node.second.outgoing_edges) {
            if (i.second.size() >= 2) {
                std::string node1 = node.first;
                std::string node2 = i.first;

                // if cannot pass reverse complementary check, this operation cannot be performed simultaneously
                if (this->check_reverse_bulge(node1, node2, check_overlap) == false)
                    continue;
                // collapse forward bulge
                std::string seq1 = this->collapse_bulge(node1, node2, removed_bulges, similarity);

                // collapse reverse complementary bulge
                if (reverse_complementary_node(node2) != node1 || reverse_complementary_node(node1) != node2) {
                    std::string seq2 = this->collapse_bulge(reverse_complementary_node(node2), reverse_complementary_node(node1), removed_bulges, similarity);
                    if (seq1 != reverse_complementary(seq2))
                        std::cout << "The two bulges did not result in reverse complementary sequence" << std::endl;
                }
            }
        }
    }
    this->merge_non_branching_paths();
}

void Graph::remove_whirl(Path& unambiguous_path, std::vector<std::string>& nodes_to_remove) {
    Path path_replace;
    this->add_node_to_path(path_replace, unambiguous_path.nodes[0]);
    this->add_node_to_path(path_replace, unambiguous_path.nodes[0]);
    if (!this->replace_genome_path_with_path(unambiguous_path, path_replace, true)) {
        std::cout << "General whirl (min multi " << unambiguous_path.multiplicity << ", color " << unambiguous_path.color << ", length " << unambiguous_path.length << "): " << unambiguous_path.nodes.at(0);
        for (int i = 1; i < unambiguous_path.nodes.size(); ++i) {
            std::cout << "->" << unambiguous_path.nodes.at(i);
        }
        std::cout << std::endl;
        return;
    }

    std::cout << "General whirl (min multi " << unambiguous_path.multiplicity << ", color " << unambiguous_path.color << ", length " << unambiguous_path.length << "): " << unambiguous_path.nodes.at(0);
    // change multiplicity of cyclic unambiguous path
    for (int i = 1; i < unambiguous_path.nodes.size(); ++i) {
        std::cout << "->" << unambiguous_path.nodes.at(i);
        int bulge_leg = unambiguous_path.bulge_legs[i - 1];

        graph[unambiguous_path.nodes[i - 1]].outgoing_edges[unambiguous_path.nodes[i]].at(bulge_leg).remove_multi_from_path(unambiguous_path);
        graph[unambiguous_path.nodes[i]].incoming_edges[unambiguous_path.nodes[i - 1]].at(bulge_leg).remove_multi_from_path(unambiguous_path);

        // if multiplicity is reduced to 0, remove edge
        if (graph[unambiguous_path.nodes[i - 1]].outgoing_edges[unambiguous_path.nodes[i]].at(bulge_leg).multiplicity == 0) {
            graph[unambiguous_path.nodes[i - 1]].outgoing_edges[unambiguous_path.nodes[i]].erase(graph[unambiguous_path.nodes[i - 1]].outgoing_edges[unambiguous_path.nodes[i]].begin() + bulge_leg);
            graph[unambiguous_path.nodes[i]].incoming_edges[unambiguous_path.nodes[i - 1]].erase(graph[unambiguous_path.nodes[i]].incoming_edges[unambiguous_path.nodes[i - 1]].begin() + bulge_leg);
            if (graph[unambiguous_path.nodes[i - 1]].outgoing_edges[unambiguous_path.nodes[i]].empty()) {
                graph[unambiguous_path.nodes[i - 1]].outgoing_edges.erase(unambiguous_path.nodes[i]);
                graph[unambiguous_path.nodes[i]].incoming_edges.erase(unambiguous_path.nodes[i - 1]);
                // add nodes without edges to remove list; we should never remove node.first because there is a self-loop
                if (graph[unambiguous_path.nodes[i - 1]].outgoing_edges.size() == 0 && graph[unambiguous_path.nodes[i - 1]].incoming_edges.size() == 0)
                    nodes_to_remove.push_back(unambiguous_path.nodes[i - 1]);
            }
        }
    }
    std::cout << std::endl;

    // add a self-loop to node.first
    Edge edge = Edge(unambiguous_path.sequence.at(k), unambiguous_path.length, unambiguous_path.sequence, unambiguous_path.multi11, unambiguous_path.multi12, unambiguous_path.multi21, unambiguous_path.multi22, unambiguous_path.multi31, unambiguous_path.multi32);
    graph[unambiguous_path.nodes[0]].outgoing_edges[unambiguous_path.nodes[0]].push_back(edge);
    graph[unambiguous_path.nodes[0]].incoming_edges[unambiguous_path.nodes[0]].push_back(edge);
    // update genome paths
    this->replace_genome_path_with_path(unambiguous_path, path_replace);
}

// remove general whirls
void Graph::general_whirl_removal(unsigned& removed_whirls, bool simple_whirl) {
    unsigned removed_bulges = 0;
    this->multi_bulge_removal(removed_bulges, true, sim_simple);
    std::vector<std::string> nodes_to_remove;

    removed_whirls = 0;
    for (auto&& node : graph) {
        // std::cout << node.first << std::endl;
        std::vector<std::string> outgoing_edges;

        // skip self loops on the start node
        if (graph[node.first].outgoing_edges.find(node.first) != graph[node.first].outgoing_edges.end())
            continue;

        for (auto&& successor : node.second.outgoing_edges) {
            // skip bulges for now
            if (successor.second.size() == 1)
                outgoing_edges.push_back(successor.first);
        }

        for (auto&& successor : outgoing_edges) {
            Path unambiguous_path;
            this->add_node_to_path(unambiguous_path, node.first);
            std::string successor_node = successor;

            std::string prev_node = node.first;
            // successor_node has unambiguous outgoing edge without or with a self-loop
            while ((graph[successor_node].outgoing_edges.size() == 1 && graph[successor_node].outgoing_edges.find(successor_node) == graph[successor_node].outgoing_edges.end()) || (graph[successor_node].outgoing_edges.size() == 2 && graph[successor_node].outgoing_edges.find(successor_node) != graph[successor_node].outgoing_edges.end())) {
                std::string node_after_successor, self_loop;
                // process self-loop
                for (auto&& n : graph[successor_node].outgoing_edges) {
                    if (n.first != successor_node)
                        node_after_successor = n.first;
                    else
                        self_loop = n.first;
                }
                // the unambiguous path cannot pass check
                if (!this->add_node_to_path(unambiguous_path, successor_node))
                    break;

                Path path_check = unambiguous_path;
                if (!this->add_node_to_path(path_check, node_after_successor) && !self_loop.empty())
                    node_after_successor = self_loop;

                // skip bulges for now
                if (this->graph[successor_node].outgoing_edges[node_after_successor].size() > 1)
                    break;

                // find a cyclic unambiguous path
                if (node_after_successor == node.first) {
                    if (!this->add_node_to_path(unambiguous_path, node_after_successor))
                        break;

                    if (simple_whirl && unambiguous_path.nodes.size() > 3)
                        break;

                    this->remove_whirl(unambiguous_path, nodes_to_remove);

                    // process the complementary non-whirl regardless of whether it forms a whirl
                    Path unambiguous_path_reverse;
                    assert(get_reverse_path(unambiguous_path, unambiguous_path_reverse));
                    this->remove_whirl(unambiguous_path_reverse, nodes_to_remove);

                    removed_whirls += 2;
                    break;
                }

                // update previous node and successor node
                prev_node = successor_node;
                successor_node = node_after_successor;
            }
        }
    }

    for (auto&& n : nodes_to_remove) {
        this->graph.erase(n);
        this->erase_node_from_genome_paths(n);
    }
    this->merge_non_branching_paths();
}

double Graph::PI_edlib(std::string sequence1, std::string sequence2) {
    EdlibAlignResult result = edlibAlign(sequence1.c_str(), sequence1.size(), sequence2.c_str(), sequence2.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
        std::string cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
        edlibFreeAlignResult(result);
        for (auto&& c : cigar) {
            if (c == '=')
                c = 'M';
        }
        double identity = percent_identity(cigar);
        return identity;
    }
    else {
        std::cout << "edlib failed" << std::endl;
        return -1;
    }
}

std::string Graph::collapse_complex_bulge_two_multi_edge_paths(Path& p1, Path& p2, double max_identity, std::vector<std::string>& nodes_to_remove, bool verbose, std::string output) {
    Path path_replace;
    this->add_node_to_path(path_replace, p1.nodes.at(0));
    this->add_node_to_path(path_replace, p1.nodes.at(p1.nodes.size() - 1));
    if (!this->replace_genome_path_with_path(p1, path_replace, true) || !this->replace_genome_path_with_path(p2, path_replace, true))
        return "";

    std::cout << "Complex bulge: path1 (length " << p1.length << ", min multi " << p1.multiplicity << ", color " << p1.color << "): " << p1.nodes[0];
    for (int i = 1; i < p1.nodes.size();++i)
        std::cout << "->" << p1.nodes[i];
    std::cout << "; path2 (length " << p2.length << ", min multi " << p2.multiplicity << ", color " << p2.color << "): " << p2.nodes[0];
    for (int i = 1; i < p2.nodes.size();++i)
        std::cout << "->" << p2.nodes[i];
    std::cout << "; " << "identity " << max_identity << std::endl;

    // decrease multiplicity for p1
    for (int i = 1; i < p1.nodes.size();++i) {
        int bulge_leg = p1.bulge_legs.at(i - 1);
        graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].at(bulge_leg).remove_multi_from_path(p1);
        graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].at(bulge_leg).remove_multi_from_path(p1);
        if (graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].at(bulge_leg).multiplicity == 0) {
            graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].erase(graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].begin() + bulge_leg);
            graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].erase(graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].begin() + bulge_leg);
        }
        if (graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].empty()) {
            graph[p1.nodes.at(i - 1)].outgoing_edges.erase(p1.nodes.at(i));
            graph[p1.nodes.at(i)].incoming_edges.erase(p1.nodes.at(i - 1));
            if (p1.nodes.at(i - 1) != p1.nodes.at(0) && graph[p1.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p1.nodes.at(i - 1)].outgoing_edges.size() == 0)
                nodes_to_remove.emplace_back(p1.nodes.at(i - 1));
        }
    }

    // decrease multiplicity for p2
    for (int i = 1; i < p2.nodes.size();++i) {
        int bulge_leg = p2.bulge_legs.at(i - 1);
        graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].at(bulge_leg).remove_multi_from_path(p2);
        graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].at(bulge_leg).remove_multi_from_path(p2);
        if (graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].at(bulge_leg).multiplicity == 0) {
            graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].erase(graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].begin() + bulge_leg);
            graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].erase(graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].begin() + bulge_leg);
        }
        if (graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].empty()) {
            graph[p2.nodes.at(i - 1)].outgoing_edges.erase(p2.nodes.at(i));
            graph[p2.nodes.at(i)].incoming_edges.erase(p2.nodes.at(i - 1));
            if (p2.nodes.at(i - 1) != p2.nodes.at(0) && graph[p2.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p2.nodes.at(i - 1)].outgoing_edges.size() == 0)
                nodes_to_remove.emplace_back(p2.nodes.at(i - 1));
        }
    }

    // add an edge from node.first to node.last with multiplicity p1 + p2
    std::string seq;
    if (p1.length > p2.length) {
        Edge edge(p1.sequence.at(k), p1.length, p1.sequence, p1.multi11, p1.multi12, p1.multi21, p1.multi22, p1.multi31, p1.multi32);
        edge.add_multi_from_edge_or_path(p2);
        graph[p1.nodes.at(0)].outgoing_edges[p1.nodes.at(p1.nodes.size() - 1)].emplace_back(edge);
        graph[p1.nodes.at(p1.nodes.size() - 1)].incoming_edges[p1.nodes.at(0)].emplace_back(edge);
        seq = edge.sequence;
    }
    else {
        Edge edge(p2.sequence.at(k), p2.length, p2.sequence, p2.multi11, p2.multi12, p2.multi21, p2.multi22, p2.multi31, p2.multi32);
        edge.add_multi_from_edge_or_path(p1);
        graph[p1.nodes.at(0)].outgoing_edges[p1.nodes.at(p1.nodes.size() - 1)].emplace_back(edge);
        graph[p1.nodes.at(p1.nodes.size() - 1)].incoming_edges[p1.nodes.at(0)].emplace_back(edge);
        seq = edge.sequence;
    }

    // if it forms a simple bulge, remove it immediately
    if (graph[p1.nodes.at(0)].outgoing_edges[p1.nodes.at(p1.nodes.size() - 1)].size() == 2 && p1.nodes.at(0) != reverse_complementary_node(p1.nodes.at(p1.nodes.size() - 1))) {
        unsigned removed_bulges = 0;
        seq = this->collapse_bulge(p1.nodes.at(0), p1.nodes.at(p1.nodes.size() - 1), removed_bulges, sim_simple);
    }

    this->replace_genome_path_with_path(p1, path_replace);
    this->replace_genome_path_with_path(p2, path_replace);

    if (verbose)
        write_graph(output + "/graph.multi_edge." + p1.nodes.at(0) + "_" + p1.nodes.at(p1.nodes.size() - 1));
    return seq;
}

// collapse bulges with two paths of multiple edges (at most x)
void Graph::resolving_bulge_with_two_multi_edge_paths(unsigned& removed_paths, int x, double identity, bool verbose, std::string output) {
    unsigned removed_bulges = 0;
    this->multi_bulge_removal(removed_bulges, true, sim_simple);

    removed_paths = 0;
    std::vector<std::string> nodes_to_remove;
    std::vector<Bulge> bulges;
    for (auto&& node : graph) {
        if (node.second.outgoing_edges.size() <= 1)
            continue;
        //store all paths here
        std::unordered_map<std::string, std::vector<Path>> all_paths_ending_at_key;
        std::queue<Path> paths_bfs;
        Path path;
        this->add_node_to_path(path, node.first);
        paths_bfs.push(path);

        while (paths_bfs.front().nodes.size() <= x) {
            if (paths_bfs.empty()) break;
            Path prev_path = paths_bfs.front();
            paths_bfs.pop();

            std::string prev_node = prev_path.nodes.at(prev_path.nodes.size() - 1);
            for (auto&& successor : this->graph[prev_node].outgoing_edges) {
                // ignore back edges, this will ignore self - loops as well
                // bool flag = false;
                // for (auto&& i : prev_path.nodes)
                //     if (successor.first == i)
                //         flag = true;
                // if (flag) continue;

                for (int i = 0; i < successor.second.size(); ++i) {
                    Path current_path = prev_path;
                    if (!this->add_node_to_path(current_path, successor.first, i))
                        continue;

                    paths_bfs.push(current_path);
                    all_paths_ending_at_key[successor.first].emplace_back(current_path);
                }
            }
        }

        // process the paths
        for (auto&& ending_node : all_paths_ending_at_key) {
            // check whether there are multiple paths
            if (ending_node.second.size() < 2)
                continue;

            // store all the remaining paths
            std::vector<Path> paths;
            for (auto&& path : ending_node.second) {
                // check the existence of path
                bool flag = true;
                Path path_new;
                for (int i = 0; i < path.nodes.size(); ++i) {
                    if (i == 0) {
                        if (!this->add_node_to_path(path_new, path.nodes[i]))
                            flag = false;
                    }
                    else if (!this->add_node_to_path(path_new, path.nodes[i], path.bulge_legs[i - 1]))
                        flag = false;
                }
                if (flag)
                    paths.emplace_back(path_new);
            }
            // check whether there are multiple remaining paths
            if (paths.size() < 2)
                continue;

            // select two paths with the highest identity
            double max_identity = 0;
            int max_lcs_len = 0;
            Path p1, p2;
            for (int i = 0;i < paths.size();++i) {
                std::set<std::string> s1;
                for (int k = 1; k < paths[i].nodes.size() - 1; ++k)
                    s1.insert(paths[i].nodes[k]);
                for (int j = i + 1; j < paths.size(); ++j) {
                    std::string sequence1 = paths[i].sequence, sequence2 = paths[j].sequence;
                    // the two paths should have no shared nodes
                    bool flag = false;
                    for (auto&& n : paths[j].nodes) {
                        if (s1.find(n) != s1.end())
                            flag = true;
                    }
                    if (flag) continue;

                    std::set<std::string> s2;
                    for (int k = 1; k < paths[j].nodes.size() - 1; ++k)
                        s2.insert(paths[j].nodes[k]);
                    for (auto&& n : paths[i].nodes) {
                        if (s2.find(n) != s2.end())
                            flag = true;
                    }
                    if (flag) continue;

                    // if the two sequences are reverse complementary, pass collapsing
                    if (paths[i].sequence == reverse_complementary(paths[j].sequence))
                        continue;

                    std::cout << "Run edlib for path1 (length " << paths[i].length << "): " << paths[i].nodes[0];
                    for (int k = 1; k < paths[i].nodes.size();++k)
                        std::cout << "->" << paths[i].nodes[k];
                    std::cout << "; path2 (length " << paths[j].length << "): " << paths[j].nodes[0];
                    for (int k = 1; k < paths[j].nodes.size();++k)
                        std::cout << "->" << paths[j].nodes[k];

                    double alignment_identity = 0;
                    if (sequence1.size() > sequence2.size() && 1.0 * (sequence1.size() - sequence2.size()) / sequence1.size() > 1 - identity) {
                        std::cout << std::endl;
                        continue;
                    }
                    if (sequence2.size() > sequence1.size() && 1.0 * (sequence2.size() - sequence1.size()) / sequence2.size() > 1 - identity) {
                        std::cout << std::endl;
                        continue;
                    }

                    alignment_identity = PI_edlib(sequence1, sequence2);

                    std::cout << "; identity " << alignment_identity << std::endl;
                    if (alignment_identity > max_identity) {
                        p1 = paths[i];
                        p2 = paths[j];
                        max_identity = alignment_identity;
                    }
                }
            }
            // the two paths cannot pass similarity check
            if (p1.nodes.empty() || p2.nodes.empty())
                continue;
            if (max_identity < identity)
                continue;

            // skip palindromic multi-edge bulges
            if (p1.sequence == reverse_complementary(p2.sequence))
                continue;

            if (p1.nodes.size() == 2 && p2.nodes.size() == 2)
                continue;

            Bulge b(p1, p2, max_identity);
            bulges.emplace_back(b);
        }
    }

    std::sort(bulges.begin(), bulges.end(),
        [](const Bulge& lhs, const Bulge& rhs) {
            if (lhs.num_edges < rhs.num_edges)
                return true;
            if (lhs.num_edges > rhs.num_edges)
                return false;
            if (lhs.max_identity > rhs.max_identity)
                return true;
            if (lhs.max_identity < rhs.max_identity)
                return false;
            return lhs.leg1.nodes < rhs.leg1.nodes;
        }
    );

    for (auto&& bulge : bulges) {
        Path p1, p2;
        bool flag = true;
        // check forward paths
        for (int i = 0; i < bulge.leg1.nodes.size(); ++i) {
            if (i == 0) {
                if (!this->add_node_to_path(p1, bulge.leg1.nodes[i])) {
                    flag = false;
                    break;
                }
            }
            else if (!this->add_node_to_path(p1, bulge.leg1.nodes[i], bulge.leg1.bulge_legs[i - 1])) {
                flag = false;
                break;
            }
        }
        for (int i = 0; i < bulge.leg2.nodes.size(); ++i) {
            if (i == 0) {
                if (!this->add_node_to_path(p2, bulge.leg2.nodes[i])) {
                    flag = false;
                    break;
                }
            }
            else if (!this->add_node_to_path(p2, bulge.leg2.nodes[i], bulge.leg2.bulge_legs[i - 1])) {
                flag = false;
                break;
            }
        }

        if (!flag)
            continue;

        std::string seq1 = this->collapse_complex_bulge_two_multi_edge_paths(p1, p2, bulge.max_identity, nodes_to_remove, verbose, output);

        Path p1_reverse, p2_reverse;
        assert(get_reverse_path(p1, p1_reverse));
        assert(get_reverse_path(p2, p2_reverse));
        std::string seq2 = this->collapse_complex_bulge_two_multi_edge_paths(p1_reverse, p2_reverse, bulge.max_identity, nodes_to_remove, verbose, output);

        if (seq1 != reverse_complementary(seq2))
            std::cout << "The two bulges did not result in reverse complementary sequence" << std::endl;

        if (seq1 != "" && seq2 != "")
            removed_paths += 2;

        for (auto&& node : nodes_to_remove) {
            this->graph.erase(node);
            this->erase_node_from_genome_paths(node);
        }
        nodes_to_remove.clear();
        merge_non_branching_paths();
    }
}