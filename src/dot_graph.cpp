#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "dot_graph.hpp"
#include <cassert>
#include <algorithm>
#include <iomanip>

Edge::Edge() {}

Edge::Edge(char& base, const unsigned& length, const std::string& sequence, int multi11, int multi12, int multi21, int multi22, int multi31, int multi32) {
    this->start_base = base;
    this->length = length;
    this->sequence = sequence;
    this->multi11 = multi11;
    this->multi12 = multi12;
    this->multi21 = multi21;
    this->multi22 = multi22;
    this->multi31 = multi31;
    this->multi32 = multi32;
    this->update_color_by_multi();
    this->original_edge = false;
}

void Edge::update_color_by_multi() {
    bool flag1 = false, flag2 = false, flag3 = false;
    this->multiplicity = this->multi11 + this->multi12 + this->multi21 + this->multi22 + this->multi31 + this->multi32;
    if (this->multi11 || this->multi12)
        flag1 = true;
    if (this->multi21 || this->multi22)
        flag2 = true;
    if (this->multi31 || this->multi32)
        flag3 = true;
    if (flag1 + flag2 + flag3 >= 2)
        this->color = 0;
    else if (flag1)
        this->color = 1;
    else if (flag2)
        this->color = 2;
    else if (flag3)
        this->color = 3;
    else {
        if (this->multiplicity != 0)
            std::cout << "Error: edge cannot be assigned color (update edge)" << std::endl;
        this->color = -1;
    }
}

void Edge::add_multi_from_edge_or_path(Edge& edge_or_path) {
    this->multi11 += edge_or_path.multi11;
    this->multi12 += edge_or_path.multi12;
    this->multi21 += edge_or_path.multi21;
    this->multi22 += edge_or_path.multi22;
    this->multi31 += edge_or_path.multi31;
    this->multi32 += edge_or_path.multi32;
    this->update_color_by_multi();
    this->original_edge = false;
}

void Edge::add_multi_from_edge_or_path(Path& edge_or_path) {
    this->multi11 += edge_or_path.multi11;
    this->multi12 += edge_or_path.multi12;
    this->multi21 += edge_or_path.multi21;
    this->multi22 += edge_or_path.multi22;
    this->multi31 += edge_or_path.multi31;
    this->multi32 += edge_or_path.multi32;
    this->update_color_by_multi();
    this->original_edge = false;
}

void Edge::remove_multi_from_path(Path& edge_or_path) {
    this->multi11 -= edge_or_path.multi11;
    this->multi12 -= edge_or_path.multi12;
    this->multi21 -= edge_or_path.multi21;
    this->multi22 -= edge_or_path.multi22;
    this->multi31 -= edge_or_path.multi31;
    this->multi32 -= edge_or_path.multi32;
    assert(this->multi11 >= 0 && this->multi12 >= 0 && this->multi21 >= 0 && this->multi22 >= 0 && this->multi31 >= 0 && this->multi32 >= 0);
    this->update_color_by_multi();
    this->original_edge = false;
}

bool Edge::check_equal_multi(Edge& edge_or_path) {
    return this->multi11 == edge_or_path.multi11 && this->multi12 == edge_or_path.multi12 && this->multi21 == edge_or_path.multi21 && this->multi22 == edge_or_path.multi22 && this->multi31 == edge_or_path.multi31 && this->multi32 == edge_or_path.multi32;
}

bool Edge::check_equal_multi(Path& edge_or_path) {
    return this->multi11 == edge_or_path.multi11 && this->multi12 == edge_or_path.multi12 && this->multi21 == edge_or_path.multi21 && this->multi22 == edge_or_path.multi22 && this->multi31 == edge_or_path.multi31 && this->multi32 == edge_or_path.multi32;
}

bool Edge::check_contain_multi(Path& edge_or_path) {
    return this->multi11 >= edge_or_path.multi11 && this->multi12 >= edge_or_path.multi12 && this->multi21 >= edge_or_path.multi21 && this->multi22 >= edge_or_path.multi22 && this->multi31 >= edge_or_path.multi31 && this->multi32 >= edge_or_path.multi32;
}

bool Edge::check_folds_multi(Edge& edge_or_path) {
    if (this->multi11 != 0 && edge_or_path.multi11 % this->multi11 != 0)
        return false;
    if (this->multi12 != 0 && edge_or_path.multi12 % this->multi12 != 0)
        return false;
    if (this->multi21 != 0 && edge_or_path.multi21 % this->multi21 != 0)
        return false;
    if (this->multi22 != 0 && edge_or_path.multi22 % this->multi22 != 0)
        return false;
    if (this->multi31 != 0 && edge_or_path.multi31 % this->multi31 != 0)
        return false;
    if (this->multi32 != 0 && edge_or_path.multi32 % this->multi32 != 0)
        return false;

    int fold = -1;
    if (this->multi11 != 0) {
        int folds1 = edge_or_path.multi11 / this->multi11;
        fold = folds1;
    }
    else if (edge_or_path.multi11 != 0)
        return false;

    if (this->multi12 != 0) {
        int folds2 = edge_or_path.multi12 / this->multi12;
        if (fold != -1 && folds2 != fold)
            return false;
        else if (fold == -1)
            fold = folds2;
    }
    else if (edge_or_path.multi12 != 0)
        return false;

    if (this->multi21 != 0) {
        int folds3 = edge_or_path.multi21 / this->multi21;
        if (fold != -1 && folds3 != fold)
            return false;
        else if (fold == -1)
            fold = folds3;
    }
    else if (edge_or_path.multi21 != 0)
        return false;

    if (this->multi22 != 0) {
        int folds4 = edge_or_path.multi22 / this->multi22;
        if (fold != -1 && folds4 != fold)
            return false;
        else if (fold == -1)
            fold = folds4;
    }
    else if (edge_or_path.multi22 != 0)
        return false;

    if (this->multi31 != 0) {
        int folds5 = edge_or_path.multi31 / this->multi31;
        if (fold != -1 && folds5 != fold)
            return false;
        else if (fold == -1)
            fold = folds5;
    }
    else if (edge_or_path.multi31 != 0)
        return false;

    if (this->multi32 != 0) {
        int folds6 = edge_or_path.multi32 / this->multi32;
        if (fold != -1 && folds6 != fold)
            return false;
        else if (fold == -1)
            fold = folds6;
    }
    else if (edge_or_path.multi32 != 0)
        return false;

    return true;
}

Node::Node() {}

Path::Path() {}

void Path::update_color_by_multi() {
    bool flag1 = false, flag2 = false, flag3 = false;
    if (this->multi11 || this->multi12)
        flag1 = true;
    if (this->multi21 || this->multi22)
        flag2 = true;
    if (this->multi31 || this->multi32)
        flag3 = true;
    if (flag1 + flag2 + flag3 >= 2)
        this->color = 0;
    else if (flag1)
        this->color = 1;
    else if (flag2)
        this->color = 2;
    else if (flag3)
        this->color = 3;
    else {
        if (this->multiplicity != 0)
            std::cout << "Error: edge cannot be assigned color (update path)" << std::endl;
        this->color = -1;
    }
    this->multiplicity = this->multi11 + this->multi12 + this->multi21 + this->multi22 + this->multi31 + this->multi32;
}

bool Path::update_min_multi(Edge edge_or_path) {
    if (this->multiplicity == -1) {
        this->multi11 = edge_or_path.multi11;
        this->multi12 = edge_or_path.multi12;
        this->multi21 = edge_or_path.multi21;
        this->multi22 = edge_or_path.multi22;
        this->multi31 = edge_or_path.multi31;
        this->multi32 = edge_or_path.multi32;
    }
    else {
        this->multi11 = std::min(this->multi11, edge_or_path.multi11);
        this->multi12 = std::min(this->multi12, edge_or_path.multi12);
        this->multi21 = std::min(this->multi21, edge_or_path.multi21);
        this->multi22 = std::min(this->multi22, edge_or_path.multi22);
        this->multi31 = std::min(this->multi31, edge_or_path.multi31);
        this->multi32 = std::min(this->multi32, edge_or_path.multi32);
    }

    this->multiplicity = this->multi11 + this->multi12 + this->multi21 + this->multi22 + this->multi31 + this->multi32;

    if (this->multiplicity == 0)
        return false;
    this->update_color_by_multi();
    return true;
}

Bulge::Bulge() {}

Bulge::Bulge(Path p1, Path p2, double identity) {
    this->leg1 = p1;
    this->leg2 = p2;
    this->max_identity = identity;
    this->num_edges = p1.nodes.size() + p2.nodes.size();
}

Genome::Genome() {}

void Genome::add_node_to_end(std::string node) {
    this->genome_path.push_back(node);
}

bool Genome::find_path(Path& path, std::vector<int>& pos) {
    pos.clear();
    if (path.nodes.empty())
        return false;
    std::string node = path.nodes.at(0);
    auto iter = std::find(this->genome_path.begin(), this->genome_path.end(), node);
    while (iter != this->genome_path.end()) {
        bool flag = true;
        for (int i = 1; i < path.nodes.size(); ++i) {
            if (iter + i == this->genome_path.end()) {
                flag = false;
                break;
            }
            if (path.nodes.at(i) != *(iter + i)) {
                flag = false;
                break;
            }
        }
        if (flag) {
            pos.push_back(iter - this->genome_path.begin());
            iter = std::find(iter + path.nodes.size() - 1, this->genome_path.end(), node);
        }
        else
            iter = std::find(iter + 1, this->genome_path.end(), node);
    }
    if (pos.empty())
        return false;
    else
        return true;
}

void Genome::construct_edge_orders() {
    this->node2order.clear();
    for (int i = 0; i + 1 < this->genome_path.size(); ++i) {
        this->node2order[this->genome_path.at(i)][this->genome_path.at(i + 1)].push_back(i);
    }
}

void Genome::erase_node(std::string node) {
    auto iter = std::find(this->genome_path.begin(), this->genome_path.end(), node);
    while (iter != this->genome_path.end()) {
        this->genome_path.erase(iter);
        iter = std::find(this->genome_path.begin(), this->genome_path.end(), node);
    }
}

Graph::Graph() {}

int Graph::get_num_nodes() {
    return this->graph.size();
}

void Graph::erase_node_from_genome_paths(std::string node) {
    this->genome_path_haplome1_strand1.erase_node(node);
    this->genome_path_haplome1_strand2.erase_node(node);
    this->genome_path_haplome2_strand1.erase_node(node);
    this->genome_path_haplome2_strand2.erase_node(node);
    this->genome_path_haplome3_strand1.erase_node(node);
    this->genome_path_haplome3_strand2.erase_node(node);
}

void Graph::get_locations(std::string& ref, std::string& seq, std::vector<int>& locations) {
    size_t pos = 0;
    while (ref.find(seq, pos) != std::string::npos) {
        pos = ref.find(seq, pos);
        locations.push_back(pos);
        pos += 1;
    }
    // // test if the seq is start
    // if (ref.find(seq.substr(k)) == 0) {
    //     locations.push_back(0);
    // }
    // // test if the seq is end
    // if (ref.find_last_of(seq.substr(0, seq.size() - k)) == ref.size() - seq.size() + k) {
    //     locations.push_back(ref.size() - seq.size() + k);
    // }
}

bool Graph::replace_genome_path_with_path(Path& path, Path& path_replace, bool only_check) {
    this->find_path_from_genome_paths(path);
    int maximum_replace = std::min(path.multiplicity, int(path.index_haplome1_strand1.size() + path.index_haplome1_strand2.size() + path.index_haplome2_strand1.size() + path.index_haplome2_strand2.size() + path.index_haplome3_strand1.size() + path.index_haplome3_strand2.size()));
    assert(maximum_replace == path.multiplicity);
    if (maximum_replace < int(path.index_haplome1_strand1.size() + path.index_haplome1_strand2.size() + path.index_haplome2_strand1.size() + path.index_haplome2_strand2.size() + path.index_haplome3_strand1.size() + path.index_haplome3_strand2.size())) {
        std::cout << "Path multiplicity less than genome locations: path multiplicity " << path.multiplicity << " genome locatons " << path.index_haplome1_strand1.size() << " " << path.index_haplome1_strand2.size() << " " << path.index_haplome2_strand1.size() << " " << path.index_haplome2_strand2.size() << " " << path.index_haplome3_strand1.size() << " " << path.index_haplome3_strand2.size() << std::endl;
        return false;
    }
    if (only_check)
        return true;
    if (path.index_haplome1_strand1.size()) {
        int off_set = 0;
        for (auto&& pos : path.index_haplome1_strand1) {
            if (maximum_replace == 0)
                break;
            maximum_replace -= 1;
            auto iter_start = this->genome_path_haplome1_strand1.genome_path.begin() + pos + off_set;
            auto iter_end = iter_start + path.nodes.size();
            this->genome_path_haplome1_strand1.genome_path.erase(iter_start, iter_end);
            this->genome_path_haplome1_strand1.genome_path.insert(iter_start, path_replace.nodes.begin(), path_replace.nodes.end());
            off_set = off_set + path_replace.nodes.size() - path.nodes.size();
        }
    }
    if (path.index_haplome1_strand2.size()) {
        int off_set = 0;
        for (auto&& pos : path.index_haplome1_strand2) {
            if (maximum_replace == 0)
                break;
            maximum_replace -= 1;
            auto iter_start = this->genome_path_haplome1_strand2.genome_path.begin() + pos + off_set;
            auto iter_end = iter_start + path.nodes.size();
            this->genome_path_haplome1_strand2.genome_path.erase(iter_start, iter_end);
            this->genome_path_haplome1_strand2.genome_path.insert(iter_start, path_replace.nodes.begin(), path_replace.nodes.end());
            off_set = off_set + path_replace.nodes.size() - path.nodes.size();
        }
    }
    if (path.index_haplome2_strand1.size()) {
        int off_set = 0;
        for (auto&& pos : path.index_haplome2_strand1) {
            if (maximum_replace == 0)
                break;
            maximum_replace -= 1;
            auto iter_start = this->genome_path_haplome2_strand1.genome_path.begin() + pos + off_set;
            auto iter_end = iter_start + path.nodes.size();
            this->genome_path_haplome2_strand1.genome_path.erase(iter_start, iter_end);
            this->genome_path_haplome2_strand1.genome_path.insert(iter_start, path_replace.nodes.begin(), path_replace.nodes.end());
            off_set = off_set + path_replace.nodes.size() - path.nodes.size();
        }
    }
    if (path.index_haplome2_strand2.size()) {
        int off_set = 0;
        for (auto&& pos : path.index_haplome2_strand2) {
            if (maximum_replace == 0)
                break;
            maximum_replace -= 1;
            auto iter_start = this->genome_path_haplome2_strand2.genome_path.begin() + pos + off_set;
            auto iter_end = iter_start + path.nodes.size();
            this->genome_path_haplome2_strand2.genome_path.erase(iter_start, iter_end);
            this->genome_path_haplome2_strand2.genome_path.insert(iter_start, path_replace.nodes.begin(), path_replace.nodes.end());
            off_set = off_set + path_replace.nodes.size() - path.nodes.size();
        }
    }
    if (path.index_haplome3_strand1.size()) {
        int off_set = 0;
        for (auto&& pos : path.index_haplome3_strand1) {
            if (maximum_replace == 0)
                break;
            maximum_replace -= 1;
            auto iter_start = this->genome_path_haplome3_strand1.genome_path.begin() + pos + off_set;
            auto iter_end = iter_start + path.nodes.size();
            this->genome_path_haplome3_strand1.genome_path.erase(iter_start, iter_end);
            this->genome_path_haplome3_strand1.genome_path.insert(iter_start, path_replace.nodes.begin(), path_replace.nodes.end());
            off_set = off_set + path_replace.nodes.size() - path.nodes.size();
        }
    }
    if (path.index_haplome3_strand2.size()) {
        int off_set = 0;
        for (auto&& pos : path.index_haplome3_strand2) {
            if (maximum_replace == 0)
                break;
            maximum_replace -= 1;
            auto iter_start = this->genome_path_haplome3_strand2.genome_path.begin() + pos + off_set;
            auto iter_end = iter_start + path.nodes.size();
            this->genome_path_haplome3_strand2.genome_path.erase(iter_start, iter_end);
            this->genome_path_haplome3_strand2.genome_path.insert(iter_start, path_replace.nodes.begin(), path_replace.nodes.end());
            off_set = off_set + path_replace.nodes.size() - path.nodes.size();
        }
    }
    return true;
}

bool Graph::get_reverse_path(Path& path, Path& path_reverse) {
    bool flag = true;
    this->add_node_to_path(path_reverse, reverse_complementary_node(path.nodes[path.nodes.size() - 1]));

    for (int i = int(path.nodes.size()) - 2; i >= 0; --i) {
        // to process palindromic simple bulge
        if (graph[path_reverse.nodes[path_reverse.nodes.size() - 1]].outgoing_edges[reverse_complementary_node(path.nodes[i])].size() > 1 && path_reverse.nodes[path_reverse.nodes.size() - 1] == path.nodes[i]) {
            std::vector<Edge>& edges = graph[path_reverse.nodes[path_reverse.nodes.size() - 1]].outgoing_edges[reverse_complementary_node(path.nodes[i])];

            std::string seq_forward = path.index2palindromic_seq.at(i);
            std::string seq_reverse = reverse_complementary(seq_forward);

            int bulge_leg = 0;
            bool find_reverse = false;
            for (int e = 0; e < edges.size(); ++e) {
                if (edges.at(e).sequence == seq_reverse) {
                    bulge_leg = e;
                    find_reverse = true;
                    break;
                }
            }
            if (!find_reverse) {
                flag = false;
                std::cout << "Error: reverse edge does not exist" << std::endl;
                break;
            }

            if (!this->add_node_to_path(path_reverse, reverse_complementary_node(path.nodes[i]), bulge_leg)) {
                flag = false;
                break;
            }

            path_reverse.index2palindromic_seq[path_reverse.nodes.size() - 2] = seq_reverse;
        }
        else if (!this->add_node_to_path(path_reverse, reverse_complementary_node(path.nodes[i]), path.bulge_legs[i])) {
            flag = false;
            break;
        }
    }
    if (!flag)
        return false;

    assert(path.multiplicity == path_reverse.multiplicity);
    assert(path.color == path_reverse.color);
    return true;
}

void Graph::find_haplome(std::string& haplome1_seq, std::string& haplome1_seq_r,
    std::string& haplome2_seq, std::string& haplome2_seq_r, std::string& haplome3_seq, std::string& haplome3_seq_r,
    std::string edge_sequence, std::string& edge, std::string& edge_r, int& haplome_index,
    std::vector<int>& positions_haplome1_strand1, std::vector<int>& positions_haplome1_strand2,
    std::vector<int>& positions_haplome2_strand1, std::vector<int>& positions_haplome2_strand2,
    std::vector<int>& positions_haplome3_strand1, std::vector<int>& positions_haplome3_strand2
) {
    bool haplome_index1 = haplome1_seq.find(edge_sequence) != std::string::npos || haplome1_seq_r.find(edge_sequence) != std::string::npos;
    bool haplome_index2 = haplome2_seq.find(edge_sequence) != std::string::npos || haplome2_seq_r.find(edge_sequence) != std::string::npos;
    bool haplome_index3 = haplome3_seq.find(edge_sequence) != std::string::npos || haplome3_seq_r.find(edge_sequence) != std::string::npos;

    if (haplome_index1 + haplome_index2 + haplome_index3 >= 2)
        haplome_index = 0;
    else if (haplome_index1)
        haplome_index = 1;
    else if (haplome_index2)
        haplome_index = 2;
    else if (haplome_index3)
        haplome_index = 3;
    else {
        std::cout << "Edge " << edge << " and " << edge_r << " cannot be found in the input haplomes." << std::endl;
        haplome_index = -1;
    }

    this->get_locations(haplome1_seq, edge_sequence, positions_haplome1_strand1);
    this->get_locations(haplome1_seq_r, edge_sequence, positions_haplome1_strand2);
    this->get_locations(haplome2_seq, edge_sequence, positions_haplome2_strand1);
    this->get_locations(haplome2_seq_r, edge_sequence, positions_haplome2_strand2);
    this->get_locations(haplome3_seq, edge_sequence, positions_haplome3_strand1);
    this->get_locations(haplome3_seq_r, edge_sequence, positions_haplome3_strand2);

    if (haplome_index == 1)
        assert(positions_haplome1_strand1.size() + positions_haplome1_strand2.size() > 0 && positions_haplome2_strand1.size() + positions_haplome2_strand2.size() == 0 && positions_haplome3_strand1.size() + positions_haplome3_strand2.size() == 0);
    else if (haplome_index == 2)
        assert(positions_haplome1_strand1.size() + positions_haplome1_strand2.size() == 0 && positions_haplome2_strand1.size() + positions_haplome2_strand2.size() > 0 && positions_haplome3_strand1.size() + positions_haplome3_strand2.size() == 0);
    else if (haplome_index == 3)
        assert(positions_haplome1_strand1.size() + positions_haplome1_strand2.size() == 0 && positions_haplome2_strand1.size() + positions_haplome2_strand2.size() == 0 && positions_haplome3_strand1.size() + positions_haplome3_strand2.size() > 0);
    // bool haplome_index1 = positions_haplome1_strand1.size() + positions_haplome1_strand2.size() > 0;
    // bool haplome_index2 = positions_haplome2_strand1.size() + positions_haplome2_strand2.size() > 0;
    // bool haplome_index3 = positions_haplome3_strand1.size() + positions_haplome3_strand2.size() > 0;

    // if (haplome_index1 + haplome_index2 + haplome_index3 >= 2)
    //     haplome_index = 0;
    // else if (haplome_index1)
    //     haplome_index = 1;
    // else if (haplome_index2)
    //     haplome_index = 2;
    // else if (haplome_index3)
    //     haplome_index = 3;
    // else {
    //     std::cout << "Edge " << edge << " and " << edge_r << " cannot be found in the input haplomes." << std::endl;
    //     haplome_index = -1;
    // }
}

bool Graph::find_path_from_genome_paths(Path& path) {
    this->genome_path_haplome1_strand1.find_path(path, path.index_haplome1_strand1);
    this->genome_path_haplome1_strand2.find_path(path, path.index_haplome1_strand2);
    this->genome_path_haplome2_strand1.find_path(path, path.index_haplome2_strand1);
    this->genome_path_haplome2_strand2.find_path(path, path.index_haplome2_strand2);
    this->genome_path_haplome3_strand1.find_path(path, path.index_haplome3_strand1);
    this->genome_path_haplome3_strand2.find_path(path, path.index_haplome3_strand2);

    path.multi11 = std::min(int(path.index_haplome1_strand1.size()), path.multi11);
    path.multi12 = std::min(int(path.index_haplome1_strand2.size()), path.multi12);
    path.multi21 = std::min(int(path.index_haplome2_strand1.size()), path.multi21);
    path.multi22 = std::min(int(path.index_haplome2_strand2.size()), path.multi22);
    path.multi31 = std::min(int(path.index_haplome3_strand1.size()), path.multi31);
    path.multi32 = std::min(int(path.index_haplome3_strand2.size()), path.multi32);

    path.multiplicity = path.multi11 + path.multi12 + path.multi21 + path.multi22 + path.multi31 + path.multi32;

    if (path.multiplicity > 0)
        path.update_color_by_multi();

    return path.multiplicity;
}

void Graph::clear_buffer(auto& buffer, const int omp_size) {
    buffer.clear();
    buffer.resize(omp_size);
}

void Graph::sort_locations(std::unordered_map<std::string, std::vector<int>>& edge2positions, Genome& genome) {
    std::vector<unsigned> postions;
    std::unordered_map<unsigned, std::string> pos2edge;
    for (auto&& edge : edge2positions) {
        std::sort(edge.second.begin(), edge.second.end());
        for (auto&& pos : edge.second) {
            postions.push_back(pos);
            pos2edge[pos] = edge.first;
        }
    }

    // sort positions
    std::unordered_map<unsigned, int> pos2id;
    std::sort(postions.begin(), postions.end());
    for (int i = 0; i < postions.size(); ++i) {
        pos2id[postions[i]] = i;
        std::string edge_index = pos2edge[postions[i]];
        size_t pos1 = edge_index.find(" -> ");
        std::string node1 = edge_index.substr(0, pos1);
        size_t pos2 = edge_index.find(" L ", pos1);
        std::string node2 = edge_index.substr(pos1 + 4, pos2 - pos1 - 4);
        if (!genome.genome_path.empty()) {
            assert(genome.genome_path.at(genome.genome_path.size() - 1) == node1);
            genome.add_node_to_end(node2);
        }
        else {
            genome.add_node_to_end(node1);
            genome.add_node_to_end(node2);
        }
    }
}

// Read graph from DOT file
void Graph::read_from_dot(const std::string& graph_dot, const std::string& graph_fasta, const std::string& haplome1, const std::string& haplome2, const std::string& haplome3, int k) {
    // load fasta sequence
    std::string line, node1, node2, length, multiplicity, start_base;
    std::unordered_map<std::string, std::string> edge2sequence;
    std::unordered_map<std::string, int> edge2haplome;
    std::unordered_map<std::string, std::vector<int>> edge2positions_haplome1_strand1;
    std::unordered_map<std::string, std::vector<int>> edge2positions_haplome1_strand2;
    std::unordered_map<std::string, std::vector<int>> edge2positions_haplome2_strand1;
    std::unordered_map<std::string, std::vector<int>> edge2positions_haplome2_strand2;
    std::unordered_map<std::string, std::vector<int>> edge2positions_haplome3_strand1;
    std::unordered_map<std::string, std::vector<int>> edge2positions_haplome3_strand2;

    this->k = k;

    // read input haplomes
    std::string haplome1_seq, haplome2_seq, haplome3_seq;
    if (!haplome1.empty()) {
        std::ifstream haplome1_f(haplome1);
        while (getline(haplome1_f, line)) {
            if (line.at(0) != '>') haplome1_seq += line;
        };
        haplome1_f.close();
        for (auto&& c : haplome1_seq) c = toupper(c);
    }

    if (!haplome2.empty()) {
        std::ifstream haplome2_f(haplome2);

        while (getline(haplome2_f, line)) {
            if (line.at(0) != '>') haplome2_seq += line;
        };
        haplome2_f.close();
        for (auto&& c : haplome2_seq) c = toupper(c);
    }

    if (!haplome3.empty()) {
        std::ifstream haplome3_f(haplome3);

        while (getline(haplome3_f, line)) {
            if (line.at(0) != '>') haplome3_seq += line;
        };
        haplome3_f.close();
        for (auto&& c : haplome3_seq) c = toupper(c);
    }

    std::string haplome1_seq_r = reverse_complementary(haplome1_seq), haplome2_seq_r = reverse_complementary(haplome2_seq), haplome3_seq_r = reverse_complementary(haplome3_seq);
    std::cout << "Finished reading haplome sequences." << std::endl;

    int node_len = 0;
    std::ifstream dot_file_check_node_len(graph_dot);
    while (getline(dot_file_check_node_len, line)) {
        // line is neithor a node nor an edge
        if (line.find('[') == std::string::npos)
            continue;
        // line is an edge, in this case, all the node should already be loaded
        if (line.find("->") != std::string::npos)
            continue;
        // line is an node
        else {
            std::string node_name = line.substr(0, line.find('[') - 1);
            if (node_name.at(0) == '-')
                node_len = std::max(node_len, int(node_name.size() - 1));
            else
                node_len = std::max(node_len, int(node_name.size()));
        }
    }
    dot_file_check_node_len.close();

    // read input fasta
    unsigned omp_index = 0;
    const int omp_size = 1000;
    // input buffers
    std::vector<std::string> edges_openmp(1000);
    std::vector<std::string> edges_r_openmp(1000);
    std::vector<std::string> edge_sequences_openmp(1000);
    // output buffers
    std::vector<int> edge_indices_openmp(1000);
    std::vector<int> edge_r_indices_openmp(1000);
    std::vector<std::vector<int>> edge_positions_haplome1_strand1_openmp(1000);
    std::vector<std::vector<int>> edge_positions_haplome1_strand2_openmp(1000);
    std::vector<std::vector<int>> edge_positions_haplome2_strand1_openmp(1000);
    std::vector<std::vector<int>> edge_positions_haplome2_strand2_openmp(1000);
    std::vector<std::vector<int>> edge_positions_haplome3_strand1_openmp(1000);
    std::vector<std::vector<int>> edge_positions_haplome3_strand2_openmp(1000);
    std::vector<std::vector<int>> edge_r_positions_haplome1_strand1_openmp(1000);
    std::vector<std::vector<int>> edge_r_positions_haplome1_strand2_openmp(1000);
    std::vector<std::vector<int>> edge_r_positions_haplome2_strand1_openmp(1000);
    std::vector<std::vector<int>> edge_r_positions_haplome2_strand2_openmp(1000);
    std::vector<std::vector<int>> edge_r_positions_haplome3_strand1_openmp(1000);
    std::vector<std::vector<int>> edge_r_positions_haplome3_strand2_openmp(1000);

    unsigned cnt_line = 0, cnt_0 = 0, cnt_1 = 0, cnt_2 = 0;
    std::ifstream fasta_file(graph_fasta);
    while (getline(fasta_file, line)) {
        // line is contig name
        if (line.at(0) == '>') {
            size_t pos1 = line.find_first_of('_');
            size_t pos2 = line.find_first_of('_', pos1 + 1);
            if (pos2 - pos1 - 1 == 2) {
                if (line.at(pos2 - 1) == '0')
                    node1 = "-0";
                else node1 = "0";
            }
            else {
                node1 = line.substr(pos2 - node_len - 1, node_len);
                node1 = node1.substr(node1.find_first_not_of('0'));
                if (line.at(pos2 - 1) == '0')
                    node1 = '-' + node1;
            }

            pos1 = pos2;
            pos2 = line.find_first_of('_', pos2 + 1);
            if (pos2 - pos1 - 1 == 2) {
                if (line.at(pos2 - 1) == '0')
                    node2 = "-0";
                else node2 = "0";
            }
            else {
                node2 = line.substr(pos2 - node_len - 1, node_len);
                node2 = node2.substr(node2.find_first_not_of('0'));
                if (line.at(pos2 - 1) == '0')
                    node2 = '-' + node2;
            }

            pos1 = pos2;
            pos2 = line.find_first_of('_', pos2 + 1);
            length = line.substr(pos1 + 1, pos2 - pos1 - 1);
            multiplicity = std::to_string(std::atoi(line.substr(pos2 + 1).c_str()));
        }
        // line is a contig
        else {
            std::string edge_sequence = line;
            std::string node1_r, node2_r;
            if (node1.at(0) == '-') node1_r = node1.substr(1);
            else node1_r = '-' + node1;
            if (node2.at(0) == '-') node2_r = node2.substr(1);
            else node2_r = '-' + node2;
            std::string edge_sequence_r = reverse_complementary(edge_sequence);

            // extract edge labels
            std::string edge = node1 + " -> " + node2 + " L " + edge_sequence.at(this->k) + std::to_string(edge_sequence.size() - this->k) + '(' + multiplicity + ')';
            std::string edge_r = node2_r + " -> " + node1_r + " L " + edge_sequence_r.at(this->k) + std::to_string(edge_sequence_r.size() - this->k) + '(' + multiplicity + ')';
            assert(edge2sequence.find(edge) == edge2sequence.end());
            assert(edge2sequence.find(edge_r) == edge2sequence.end());
            edge2sequence[edge] = edge_sequence;
            edge2sequence[edge_r] = edge_sequence_r;
            edges_openmp[omp_index] = edge;
            edges_r_openmp[omp_index] = edge_r;
            edge_sequences_openmp[omp_index] = edge_sequence;

            // start parallel
            if (++omp_index == omp_size) {
#pragma omp parallel for
                for (int i = 0; i < omp_size; i++) {
                    this->find_haplome(haplome1_seq, haplome1_seq_r, haplome2_seq, haplome2_seq_r, haplome3_seq, haplome3_seq_r, edge_sequences_openmp[i], edges_openmp[i], edges_r_openmp[i],
                        edge_indices_openmp[i], edge_positions_haplome1_strand1_openmp[i], edge_positions_haplome1_strand2_openmp[i], edge_positions_haplome2_strand1_openmp[i], edge_positions_haplome2_strand2_openmp[i], edge_positions_haplome3_strand1_openmp[i], edge_positions_haplome3_strand2_openmp[i]);
                    this->find_haplome(haplome1_seq, haplome1_seq_r, haplome2_seq, haplome2_seq_r, haplome3_seq, haplome3_seq_r, reverse_complementary(edge_sequences_openmp[i]), edges_openmp[i], edges_r_openmp[i],
                        edge_r_indices_openmp[i], edge_r_positions_haplome1_strand1_openmp[i], edge_r_positions_haplome1_strand2_openmp[i], edge_r_positions_haplome2_strand1_openmp[i], edge_r_positions_haplome2_strand2_openmp[i], edge_r_positions_haplome3_strand1_openmp[i], edge_r_positions_haplome3_strand2_openmp[i]);
                }

                for (int i = 0; i < omp_size; ++i) {
                    edge2haplome[edges_openmp[i]] = edge_indices_openmp[i];
                    edge2haplome[edges_r_openmp[i]] = edge_r_indices_openmp[i];
                    assert(edge2haplome[edges_openmp[i]] == edge2haplome[edges_r_openmp[i]]);
                    edge2positions_haplome1_strand1[edges_openmp[i]] = edge_positions_haplome1_strand1_openmp[i];
                    edge2positions_haplome1_strand2[edges_openmp[i]] = edge_positions_haplome1_strand2_openmp[i];
                    edge2positions_haplome2_strand1[edges_openmp[i]] = edge_positions_haplome2_strand1_openmp[i];
                    edge2positions_haplome2_strand2[edges_openmp[i]] = edge_positions_haplome2_strand2_openmp[i];
                    edge2positions_haplome3_strand1[edges_openmp[i]] = edge_positions_haplome3_strand1_openmp[i];
                    edge2positions_haplome3_strand2[edges_openmp[i]] = edge_positions_haplome3_strand2_openmp[i];
                    edge2positions_haplome1_strand1[edges_r_openmp[i]] = edge_r_positions_haplome1_strand1_openmp[i];
                    edge2positions_haplome1_strand2[edges_r_openmp[i]] = edge_r_positions_haplome1_strand2_openmp[i];
                    edge2positions_haplome2_strand1[edges_r_openmp[i]] = edge_r_positions_haplome2_strand1_openmp[i];
                    edge2positions_haplome2_strand2[edges_r_openmp[i]] = edge_r_positions_haplome2_strand2_openmp[i];
                    edge2positions_haplome3_strand1[edges_r_openmp[i]] = edge_r_positions_haplome3_strand1_openmp[i];
                    edge2positions_haplome3_strand2[edges_r_openmp[i]] = edge_r_positions_haplome3_strand2_openmp[i];

                    if (edge_indices_openmp[i] == 0)
                        cnt_0 += 1;
                    else if (edge_indices_openmp[i] == 1)
                        cnt_1 += 1;
                    else if (edge_indices_openmp[i] == 2)
                        cnt_2 += 1;
                }
                omp_index = 0;
                this->clear_buffer(edges_openmp, omp_size);
                this->clear_buffer(edges_r_openmp, omp_size);
                this->clear_buffer(edge_sequences_openmp, omp_size);
                this->clear_buffer(edge_indices_openmp, omp_size);
                this->clear_buffer(edge_r_indices_openmp, omp_size);
                this->clear_buffer(edge_positions_haplome1_strand1_openmp, omp_size);
                this->clear_buffer(edge_r_positions_haplome1_strand1_openmp, omp_size);
                this->clear_buffer(edge_positions_haplome1_strand2_openmp, omp_size);
                this->clear_buffer(edge_r_positions_haplome1_strand2_openmp, omp_size);
                this->clear_buffer(edge_positions_haplome2_strand1_openmp, omp_size);
                this->clear_buffer(edge_r_positions_haplome2_strand1_openmp, omp_size);
                this->clear_buffer(edge_positions_haplome2_strand2_openmp, omp_size);
                this->clear_buffer(edge_r_positions_haplome2_strand2_openmp, omp_size);
                this->clear_buffer(edge_positions_haplome3_strand1_openmp, omp_size);
                this->clear_buffer(edge_r_positions_haplome3_strand1_openmp, omp_size);
                this->clear_buffer(edge_positions_haplome3_strand2_openmp, omp_size);
                this->clear_buffer(edge_r_positions_haplome3_strand2_openmp, omp_size);
            }

            if (++cnt_line % 1000 == 0)
                std::cout << "Processed " << cnt_line << " edge sequences." << std::endl;
        }
    }
    // dealing with remaining edges in the buffer
#pragma omp parallel for
    for (int i = 0; i < omp_index; i++) {
        this->find_haplome(haplome1_seq, haplome1_seq_r, haplome2_seq, haplome2_seq_r, haplome3_seq, haplome3_seq_r, edge_sequences_openmp[i], edges_openmp[i], edges_r_openmp[i],
            edge_indices_openmp[i], edge_positions_haplome1_strand1_openmp[i], edge_positions_haplome1_strand2_openmp[i], edge_positions_haplome2_strand1_openmp[i], edge_positions_haplome2_strand2_openmp[i], edge_positions_haplome3_strand1_openmp[i], edge_positions_haplome3_strand2_openmp[i]);
        this->find_haplome(haplome1_seq, haplome1_seq_r, haplome2_seq, haplome2_seq_r, haplome3_seq, haplome3_seq_r, reverse_complementary(edge_sequences_openmp[i]), edges_openmp[i], edges_r_openmp[i],
            edge_r_indices_openmp[i], edge_r_positions_haplome1_strand1_openmp[i], edge_r_positions_haplome1_strand2_openmp[i], edge_r_positions_haplome2_strand1_openmp[i], edge_r_positions_haplome2_strand2_openmp[i], edge_r_positions_haplome3_strand1_openmp[i], edge_r_positions_haplome3_strand2_openmp[i]);
    }

    for (int i = 0; i < omp_index; ++i) {
        edge2haplome[edges_openmp[i]] = edge_indices_openmp[i];
        edge2haplome[edges_r_openmp[i]] = edge_r_indices_openmp[i];
        assert(edge2haplome[edges_openmp[i]] == edge2haplome[edges_r_openmp[i]]);
        edge2positions_haplome1_strand1[edges_openmp[i]] = edge_positions_haplome1_strand1_openmp[i];
        edge2positions_haplome1_strand2[edges_openmp[i]] = edge_positions_haplome1_strand2_openmp[i];
        edge2positions_haplome2_strand1[edges_openmp[i]] = edge_positions_haplome2_strand1_openmp[i];
        edge2positions_haplome2_strand2[edges_openmp[i]] = edge_positions_haplome2_strand2_openmp[i];
        edge2positions_haplome3_strand1[edges_openmp[i]] = edge_positions_haplome3_strand1_openmp[i];
        edge2positions_haplome3_strand2[edges_openmp[i]] = edge_positions_haplome3_strand2_openmp[i];
        edge2positions_haplome1_strand1[edges_r_openmp[i]] = edge_r_positions_haplome1_strand1_openmp[i];
        edge2positions_haplome1_strand2[edges_r_openmp[i]] = edge_r_positions_haplome1_strand2_openmp[i];
        edge2positions_haplome2_strand1[edges_r_openmp[i]] = edge_r_positions_haplome2_strand1_openmp[i];
        edge2positions_haplome2_strand2[edges_r_openmp[i]] = edge_r_positions_haplome2_strand2_openmp[i];
        edge2positions_haplome3_strand1[edges_r_openmp[i]] = edge_r_positions_haplome3_strand1_openmp[i];
        edge2positions_haplome3_strand2[edges_r_openmp[i]] = edge_r_positions_haplome3_strand2_openmp[i];

        if (edge_indices_openmp[i] == 0)
            cnt_0 += 1;
        else if (edge_indices_openmp[i] == 1)
            cnt_1 += 1;
        else if (edge_indices_openmp[i] == 2)
            cnt_2 += 1;
    }

    this->sort_locations(edge2positions_haplome1_strand1, this->genome_path_haplome1_strand1);
    this->sort_locations(edge2positions_haplome1_strand2, this->genome_path_haplome1_strand2);
    this->sort_locations(edge2positions_haplome2_strand1, this->genome_path_haplome2_strand1);
    this->sort_locations(edge2positions_haplome2_strand2, this->genome_path_haplome2_strand2);
    this->sort_locations(edge2positions_haplome3_strand1, this->genome_path_haplome3_strand1);
    this->sort_locations(edge2positions_haplome3_strand2, this->genome_path_haplome3_strand2);
    std::cout << "Read sequences from " << graph_fasta << " finished (loaded " << edge2sequence.size() << " edges, including reverse complements)." << std::endl;
    std::cout << "Edges belonging to haplome1: " << cnt_1 << ", haplome2: " << cnt_2 << ", both haplomes: " << cnt_0 << std::endl;

    // load graph
    std::ifstream dot_file(graph_dot);
    std::vector<std::string> edge_lines;
    while (getline(dot_file, line)) {
        // line is neithor a node nor an edge
        if (line.find('[') == std::string::npos)
            continue;
        // line is an edge
        if (line.find("->") != std::string::npos) {
            edge_lines.emplace_back(line);
        }
        // line is an node
        else {
            std::string node_name = line.substr(0, line.find('[') - 1);
            Node node;
            this->graph[node_name] = node;
        }
    }
    dot_file.close();

    struct Edge_RC {
        std::string line1;
        std::string line2;
    };
    std::vector<Edge_RC> edge_rcs;
    for (int i = 0; i < edge_lines.size(); ++i) {
        if (i % 2 == 0) {
            Edge_RC emp_edge;
            emp_edge.line1 = edge_lines[i];
            edge_rcs.emplace_back(emp_edge);
        }
        else {
            edge_rcs[i / 2].line2 = edge_lines[i];
        }
    }
    std::sort(edge_rcs.begin(), edge_rcs.end(), [](const Edge_RC& a, const Edge_RC& b) {
        return a.line1 < b.line1;
        });
    // edge_lines.clear();
    // for (auto&& e : edge_rcs) {
    //     edge_lines.emplace_back(e.line1);
    //     edge_lines.emplace_back(e.line2);
    // }
    // all the node should already be loaded
    for (auto&& line : edge_lines) {
        // extract node name
        size_t pos1 = line.find("->");
        std::string start_name = line.substr(1, pos1 - 3);
        size_t pos2 = line.find('[');
        std::string end_name = line.substr(pos1 + 4, pos2 - pos1 - 6);
        assert(this->graph.find(start_name) != this->graph.end() && this->graph.find(end_name) != this->graph.end());

        // parse edge label: starting base, length, multiplicity and sequence
        pos1 = line.find(", ");
        std::string edge_label = line.substr(pos2 + 8, pos1 - pos2 - 9);
        pos1 = edge_label.find('(');
        assert(pos1 != std::string::npos);
        char start_base = edge_label[0];
        unsigned length = std::atoi(edge_label.substr(2, pos1 - 2).c_str());
        unsigned multiplicity = std::atoi(edge_label.substr(pos1 + 1, edge_label.size() - pos1 - 2).c_str());

        std::string edge_index = start_name + " -> " + end_name + " L " + start_base + std::to_string(length) + '(' + std::to_string(multiplicity) + ')';
        assert(edge2sequence.find(edge_index) != edge2sequence.end());
        if (multiplicity != edge2positions_haplome1_strand1[edge_index].size() + edge2positions_haplome1_strand2[edge_index].size() + edge2positions_haplome2_strand1[edge_index].size() + edge2positions_haplome2_strand2[edge_index].size() + edge2positions_haplome3_strand1[edge_index].size() + edge2positions_haplome3_strand2[edge_index].size()) {
            std::cout << "Correct multiplicity of " << edge_index << " from " << multiplicity << " to " << edge2positions_haplome1_strand1[edge_index].size() << '+' << edge2positions_haplome1_strand2[edge_index].size() << '+' << edge2positions_haplome2_strand1[edge_index].size() << '+' << edge2positions_haplome2_strand2[edge_index].size() << '+' << edge2positions_haplome3_strand1[edge_index].size() << '+' << edge2positions_haplome3_strand2[edge_index].size() << std::endl;
            multiplicity = edge2positions_haplome1_strand1[edge_index].size() + edge2positions_haplome1_strand2[edge_index].size() + edge2positions_haplome2_strand1[edge_index].size() + edge2positions_haplome2_strand2[edge_index].size() + edge2positions_haplome3_strand1[edge_index].size() + edge2positions_haplome3_strand2[edge_index].size();
        }
        // assert(multiplicity == edge2positions_haplome1_strand1[edge_index].size() + edge2positions_haplome1_strand2[edge_index].size() + edge2positions_haplome2_strand1[edge_index].size() + edge2positions_haplome2_strand2[edge_index].size() + edge2positions_haplome3_strand1[edge_index].size() + edge2positions_haplome3_strand2[edge_index].size());
        Edge edge = Edge(start_base, length, edge2sequence[edge_index], int(edge2positions_haplome1_strand1[edge_index].size()), int(edge2positions_haplome1_strand2[edge_index].size()), int(edge2positions_haplome2_strand1[edge_index].size()), int(edge2positions_haplome2_strand2[edge_index].size()), int(edge2positions_haplome3_strand1[edge_index].size()), int(edge2positions_haplome3_strand2[edge_index].size()));
        edge.multiplicity11 = edge2positions_haplome1_strand1[edge_index];
        edge.multiplicity12 = edge2positions_haplome1_strand2[edge_index];
        edge.multiplicity21 = edge2positions_haplome2_strand1[edge_index];
        edge.multiplicity22 = edge2positions_haplome2_strand2[edge_index];
        edge.multiplicity31 = edge2positions_haplome3_strand1[edge_index];
        edge.multiplicity32 = edge2positions_haplome3_strand2[edge_index];
        edge.original_edge = true;

        this->graph[start_name].outgoing_edges[end_name].push_back(edge);
        this->graph[end_name].incoming_edges[start_name].push_back(edge);
    }
}

void Graph::mask_edges_with_multi_x(std::unordered_map<std::string, Node>& graph, int x, int& start_id) {
    std::unordered_set<std::string> nodes_scanned;
    for (auto&& node1 : graph) {
        std::string n1 = node1.first;
        if (nodes_scanned.find(n1) != nodes_scanned.end())
            continue;
        nodes_scanned.insert(n1);
        // nodes_scanned.insert(reverse_complementary_node(n1));
        for (auto&& node2 : node1.second.outgoing_edges) {
            std::string n2 = node2.first;
            for (int i = 0; i < node2.second.size(); ++i) {
                auto& edge = node2.second.at(i);
                bool flag = false;
                if (edge.block_id != 0)
                    continue;
                if (x == 1) {
                    if (edge.multi11 + edge.multi12 <= x && edge.multi21 + edge.multi22 <= x && edge.multi31 + edge.multi32 <= x)
                        flag = true;
                }
                else if (edge.max_multi_in_non_branching <= x)
                    flag = true;
                if (flag) {
                    std::cout << "Single edge: " << n1 << " " << n2 << " multi " << edge.multiplicity << " length " << edge.length << " id " << start_id << std::endl;
                    edge.block_id = start_id;
                    graph[n2].incoming_edges[n1].at(i).block_id = start_id;

                    int reverse_index = i;
                    if (n1 == reverse_complementary_node(n2)) {
                        bool flag = false;
                        for (int i_r = 0; i_r < graph[reverse_complementary_node(n2)].outgoing_edges[reverse_complementary_node(n1)].size();++i_r) {
                            // if (i_r != i && graph[reverse_complementary_node(n2)].outgoing_edges[reverse_complementary_node(n1)].at(i_r).sequence == reverse_complementary(edge.sequence)) {
                            if (graph[reverse_complementary_node(n2)].outgoing_edges[reverse_complementary_node(n1)].at(i_r).sequence == reverse_complementary(edge.sequence)) {
                                reverse_index = i_r;
                                flag = true;
                            }
                        }
                        assert(flag);
                    }
                    auto& edge_r = graph[reverse_complementary_node(n2)].outgoing_edges[reverse_complementary_node(n1)].at(reverse_index);
                    if (!(n1 == reverse_complementary_node(n2) && reverse_index == i))
                        assert(edge_r.block_id == 0);

                    if (x == 1)
                        assert(edge_r.multi11 + edge_r.multi12 <= x && edge_r.multi21 + edge_r.multi22 <= x && edge_r.multi31 + edge_r.multi32 <= x);
                    else
                        assert(edge_r.max_multi_in_non_branching <= x);
                    std::cout << "Single edge: " << reverse_complementary_node(n2) << " " << reverse_complementary_node(n1) << " multi " << edge.multiplicity << " length " << edge.length << " id " << -start_id << std::endl;

                    graph[reverse_complementary_node(n2)].outgoing_edges[reverse_complementary_node(n1)].at(reverse_index).block_id = -start_id;
                    graph[reverse_complementary_node(n1)].incoming_edges[reverse_complementary_node(n2)].at(reverse_index).block_id = -start_id;

                    this->synteny_blocks[start_id].length = edge.length;
                    this->synteny_blocks[start_id].multiplicity = edge.multiplicity;
                    this->synteny_blocks[start_id].num_edges = 1;

                    this->synteny_blocks[-start_id].length = edge.length;
                    this->synteny_blocks[-start_id].multiplicity = edge.multiplicity;
                    this->synteny_blocks[-start_id].num_edges = 1;

                    start_id++;
                }
            }
        }
    }
}

void Graph::change_multi_non_branching_path(std::unordered_map<std::string, Node>& graph, std::string node, std::vector<std::string>& non_branching_path, int& max_multi) {
    std::vector<std::string> nodes_before;
    std::vector<std::string> nodes_after;

    std::string node_current = node;
    bool is_circular = false;
    while (true) {
        std::string node_before;
        int cnt_existing = 0;
        for (auto&& n : graph[node_current].incoming_edges) {
            for (auto e : n.second)
                if (e.block_id == 0) {
                    node_before = n.first;
                    cnt_existing += 1;
                }
        }
        if (cnt_existing != 1)
            break;

        cnt_existing = 0;
        for (auto&& n : graph[node_before].outgoing_edges) {
            for (auto e : n.second)
                if (e.block_id == 0) {
                    cnt_existing += 1;
                }
        }

        nodes_before.emplace_back(node_before);
        node_current = node_before;
        if (node_current == node) {
            is_circular = true;
            break;
        }
        if (cnt_existing != 1)
            break;
    }

    if (!is_circular) {
        node_current = node;
        while (true) {
            std::string node_after;
            int cnt_existing = 0;
            for (auto&& n : graph[node_current].outgoing_edges) {
                for (auto e : n.second)
                    if (e.block_id == 0) {
                        node_after = n.first;
                        cnt_existing += 1;
                    }
            }
            if (cnt_existing != 1)
                break;

            cnt_existing = 0;
            for (auto&& n : graph[node_after].incoming_edges) {
                for (auto e : n.second)
                    if (e.block_id == 0) {
                        cnt_existing += 1;
                    }
            }

            nodes_after.emplace_back(node_after);
            node_current = node_after;
            if (cnt_existing != 1)
                break;
        }
    }

    non_branching_path.clear();
    for (int i = nodes_before.size() - 1; i >= 0; --i)
        non_branching_path.emplace_back(nodes_before.at(i));
    non_branching_path.emplace_back(node);
    for (int i = 0; i < nodes_after.size(); ++i)
        non_branching_path.emplace_back(nodes_after.at(i));

    max_multi = 0;
    for (int i = 0; i < non_branching_path.size() - 1; ++i) {
        for (auto&& e : graph[non_branching_path.at(i)].outgoing_edges[non_branching_path.at(i + 1)]) {
            if (e.block_id == 0) {
                if (e.max_multi_in_non_branching > max_multi)
                    max_multi = e.max_multi_in_non_branching;
            }
        }
        for (auto&& e : graph[non_branching_path.at(i + 1)].incoming_edges[non_branching_path.at(i)]) {
            if (e.block_id == 0) {
                if (e.max_multi_in_non_branching > max_multi)
                    max_multi = e.max_multi_in_non_branching;
            }
        }
    }

    for (int i = 0; i < non_branching_path.size() - 1; ++i) {
        for (auto&& e : graph[non_branching_path.at(i)].outgoing_edges[non_branching_path.at(i + 1)]) {
            if (e.block_id == 0) {
                e.max_multi_in_non_branching = max_multi;
            }
        }
        for (auto&& e : graph[non_branching_path.at(i + 1)].incoming_edges[non_branching_path.at(i)]) {
            if (e.block_id == 0) {
                e.max_multi_in_non_branching = max_multi;
            }
        }
    }
}

bool Graph::mask_non_branching_path(std::unordered_map<std::string, Node>& graph, int x, int& start_id) {
    int max_multi = 0;
    std::vector<std::string> non_branching_path, non_branching_path_r;
    std::vector<std::vector<std::string>> all_non_branching_paths;
    std::unordered_set<std::string> all_nodes_scanned;
    bool all_circular_linear = true;
    for (auto&& node : graph) {
        if (all_nodes_scanned.find(node.first) != all_nodes_scanned.end())
            continue;
        change_multi_non_branching_path(graph, node.first, non_branching_path, max_multi);
        change_multi_non_branching_path(graph, reverse_complementary_node(node.first), non_branching_path_r, max_multi);

        assert(non_branching_path.size() == non_branching_path_r.size());
        for (int i = 0;i < non_branching_path.size();++i)
            assert(non_branching_path.at(i) == reverse_complementary_node(non_branching_path_r.at(non_branching_path_r.size() - 1 - i)));

        all_nodes_scanned.insert(non_branching_path.begin(), non_branching_path.end());
        all_nodes_scanned.insert(non_branching_path_r.begin(), non_branching_path_r.end());

        std::string start_node = non_branching_path.at(0);
        std::string end_node = non_branching_path.at(non_branching_path.size() - 1);
        // circular
        if (non_branching_path.size() > 1 && start_node == end_node) {
            std::cout << "Circular: ";
            for (auto&& n : non_branching_path)
                std::cout << n << " ";
            std::cout << "multi " << max_multi << std::endl;
            all_non_branching_paths.push_back(non_branching_path);
        }
        else {
            bool flag = true;
            for (auto&& n : graph[start_node].incoming_edges) {
                for (auto&& e : n.second) {
                    if (e.block_id == 0) {
                        flag = false;
                        break;
                    }
                }
            }
            for (auto&& n : graph[end_node].outgoing_edges) {
                for (auto&& e : n.second) {
                    if (e.block_id == 0) {
                        flag = false;
                        break;
                    }
                }
            }
            // linear
            if (flag && non_branching_path.size() > 1) {
                std::cout << "Linear: ";
                for (auto&& n : non_branching_path)
                    std::cout << n << " ";
                std::cout << "multi " << max_multi << std::endl;
                all_non_branching_paths.push_back(non_branching_path);
            }
            else if (!flag) {
                all_circular_linear = false;
                if (max_multi <= x) {
                    std::cout << "Non-branching: ";
                    for (auto&& n : non_branching_path)
                        std::cout << n << " ";
                    std::cout << "multi " << max_multi << std::endl;
                    all_non_branching_paths.push_back(non_branching_path);
                }
            }
        }
    }

    for (auto&& path : all_non_branching_paths) {
        int block_length = 0;
        int multiplicity = 0;
        for (int i = 0; i < path.size() - 1; ++i) {
            bool flag = false;
            for (auto&& e : graph[path.at(i)].outgoing_edges[path.at(i + 1)]) {
                if (e.block_id == 0) {
                    if (!flag) flag = true;
                    else {
                        flag = false;
                        break;
                    }
                    block_length += e.length;
                    multiplicity = e.max_multi_in_non_branching;
                }
            }
            assert(flag);
        }

        this->synteny_blocks[start_id].length = block_length;
        this->synteny_blocks[start_id].multiplicity = multiplicity;
        this->synteny_blocks[start_id].num_edges = path.size() - 1;

        this->synteny_blocks[-start_id].length = block_length;
        this->synteny_blocks[-start_id].multiplicity = multiplicity;
        this->synteny_blocks[-start_id].num_edges = path.size() - 1;

        std::cout << "Assign id: ";
        std::cout << path.at(0);
        for (int i = 0; i < path.size() - 1; ++i) {
            for (auto&& e : graph[path.at(i)].outgoing_edges[path.at(i + 1)]) {
                std::cout << " " << path.at(i + 1);
                if (e.block_id == 0) {
                    e.block_id = start_id;
                }
            }

            for (auto&& e : graph[path.at(i + 1)].incoming_edges[path.at(i)]) {
                if (e.block_id == 0) {
                    e.block_id = start_id;
                }
            }
        }
        std::cout << " multi " << this->synteny_blocks[start_id].multiplicity << " id " << start_id << std::endl;

        std::cout << "Assign id: ";
        std::vector<std::string> path_r;
        for (int i = path.size() - 1; i >= 0;--i)
            path_r.push_back(reverse_complementary_node(path.at(i)));
        std::cout << path_r.at(0);
        for (int i = 0; i < path_r.size() - 1; ++i) {
            for (auto&& e : graph[path_r.at(i)].outgoing_edges[path_r.at(i + 1)]) {
                std::cout << " " << path_r.at(i + 1);
                if (e.block_id == 0) {
                    e.block_id = -start_id;
                }
            }

            for (auto&& e : graph[path_r.at(i + 1)].incoming_edges[path_r.at(i)]) {
                if (e.block_id == 0) {
                    e.block_id = -start_id;
                }
            }
        }
        std::cout << " multi " << this->synteny_blocks[start_id].multiplicity << " id " << -start_id << std::endl;

        start_id++;
    }
    return all_circular_linear;
}

void Graph::merge_blocks(std::vector<Block_in_Genome_Path>& blocks) {
    int cnt_edge_in_block = 0;
    int start = -1, end = -1;
    std::string node_s, node_e;
    bool include_left, include_right;
    int prev_id = 0;
    std::vector<Block_in_Genome_Path> blocks_1_final;
    for (int i = 0; i < blocks.size(); ++i) {
        if (blocks[i].id == prev_id) {
            cnt_edge_in_block += 1;
            end = blocks[i].end;
            node_e = blocks[i].node_e;
            include_right = blocks[i].include_right;
        }
        else {
            if (cnt_edge_in_block >= 1) {
                Block_in_Genome_Path tmp_path(blocks[i - 1].id, blocks[i - 1].multiplicity, blocks[i - 1].length, node_s, node_e, start, end, include_left, include_right, cnt_edge_in_block);
                tmp_path.is_short_version = true;
                blocks_1_final.emplace_back(tmp_path);
                cnt_edge_in_block = 0;
            }
            cnt_edge_in_block += 1;
            start = blocks[i].start;
            node_s = blocks[i].node_s;
            include_left = blocks[i].include_left;
            end = blocks[i].end;
            node_e = blocks[i].node_e;
            include_right = blocks[i].include_right;
            prev_id = blocks[i].id;
        }

        if (cnt_edge_in_block == blocks[i].num_edges) {
            end = blocks[i].end;
            node_e = blocks[i].node_e;
            include_right = blocks[i].include_right;

            Block_in_Genome_Path tmp_path(blocks[i].id, blocks[i].multiplicity, blocks[i].length, node_s, node_e, start, end, include_left, include_right, cnt_edge_in_block);
            blocks_1_final.emplace_back(tmp_path);
            cnt_edge_in_block = 0;
            prev_id = 0;
        }
    }
    if (cnt_edge_in_block >= 1) {
        Block_in_Genome_Path tmp_path(blocks[blocks.size() - 1].id, blocks[blocks.size() - 1].multiplicity, blocks[blocks.size() - 1].length, node_s, node_e, start, end, include_left, include_right, cnt_edge_in_block);
        tmp_path.is_short_version = true;
        blocks_1_final.emplace_back(tmp_path);
        cnt_edge_in_block = 0;
    }

    blocks.clear();
    for (auto&& b : blocks_1_final)
        blocks.emplace_back(b);
}

void Graph::synteny_block_generation(const std::string& prefix) {
    int start_id = 1;
    auto graph_synteny = graph;
    // initialize multiplicity for synteny blocks
    for (auto&& n : graph_synteny) {
        for (auto&& n_in : n.second.incoming_edges) {
            for (auto&& e : n_in.second)
                e.max_multi_in_non_branching = e.multiplicity;
        }
        for (auto&& n_out : n.second.outgoing_edges) {
            for (auto&& e : n_out.second)
                e.max_multi_in_non_branching = e.multiplicity;
        }
    }
    // remove the lowest multiplicity edges
    std::cout << "Synteny blocks multi (1+1)" << std::endl;
    mask_edges_with_multi_x(graph_synteny, 1, start_id);
    write_graph(graph_synteny, prefix + "_duplication");

    int i = 2;
    while (true) {
        std::cout << "Synteny blocks multi (" << i << ")" << std::endl;
        if (mask_non_branching_path(graph_synteny, i, start_id)) {
            break;
        }
        mask_edges_with_multi_x(graph_synteny, i, start_id);
        // write_graph(graph_synteny, prefix + "_duplication" + std::to_string(i));
        ++i;
    }
    std::cout << "Iteration stops at multi " << i << '\n' << std::endl;

    auto graph_synteny_gp = graph_synteny;

    // Haplome 1
    int len_genome1 = 0;
    int previous_multi = -1;
    std::vector<Block_in_Genome_Path> blocks_1;
    for (int i = 0; i + 1 < this->genome_path_haplome1_strand1.genome_path.size(); ++i) {
        std::string node1 = this->genome_path_haplome1_strand1.genome_path[i], node2 = this->genome_path_haplome1_strand1.genome_path[i + 1];

        // find edge index
        int edge_index = -1;
        int cnt = 0;
        for (int e = 0; e < graph_synteny_gp[node1].outgoing_edges[node2].size(); ++e) {
            if (graph_synteny_gp[node1].outgoing_edges[node2].at(e).multi11 > 0) {
                edge_index = e;
                break;
            }
        }
        int min_loc = -1;
        for (int e = 0; e < graph_synteny_gp[node1].outgoing_edges[node2].size(); ++e) {
            if (graph_synteny_gp[node1].outgoing_edges[node2].at(e).original_edge && graph_synteny_gp[node1].outgoing_edges[node2].at(e).multiplicity11.size() > 0) {
                cnt += 1;
                if (min_loc == -1 || graph_synteny_gp[node1].outgoing_edges[node2].at(e).multiplicity11.at(0) < min_loc) {
                    min_loc = graph_synteny_gp[node1].outgoing_edges[node2].at(e).multiplicity11.at(0);
                    edge_index = e;
                }
            }
        }
        assert(edge_index != -1);

        Edge current_edge = graph_synteny[node1].outgoing_edges[node2].at(edge_index);
        if (current_edge.multiplicity >= previous_multi) {
            Block_in_Genome_Path block(current_edge.block_id, this->synteny_blocks[current_edge.block_id].multiplicity, this->synteny_blocks[current_edge.block_id].length, node1, node2, len_genome1, len_genome1 + current_edge.length, true, false, this->synteny_blocks[current_edge.block_id].num_edges);
            blocks_1.emplace_back(block);
        }
        else {
            Block_in_Genome_Path block(current_edge.block_id, this->synteny_blocks[current_edge.block_id].multiplicity, this->synteny_blocks[current_edge.block_id].length, node1, node2, len_genome1 + k, len_genome1 + current_edge.length, false, false, this->synteny_blocks[current_edge.block_id].num_edges);
            blocks_1[blocks_1.size() - 1].include_right = true;
            blocks_1[blocks_1.size() - 1].end += k;
            blocks_1.emplace_back(block);
        }
        previous_multi = current_edge.multiplicity;
        len_genome1 += current_edge.length;
        graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity -= 1;
        graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multi11 -= 1;
        if (min_loc != -1) {
            graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity11.erase(graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity11.begin());
        }
    }
    if (blocks_1.size()) {
        blocks_1[blocks_1.size() - 1].include_right = true;
        blocks_1[blocks_1.size() - 1].end += k;
        merge_blocks(blocks_1);
    }

    // Haplome 2
    int len_genome2 = 0;
    previous_multi = -1;
    std::vector<Block_in_Genome_Path> blocks_2;
    for (int i = 0; i + 1 < this->genome_path_haplome2_strand1.genome_path.size(); ++i) {
        std::string node1 = this->genome_path_haplome2_strand1.genome_path[i], node2 = this->genome_path_haplome2_strand1.genome_path[i + 1];

        // find edge index
        int edge_index = -1;
        int cnt = 0;
        for (int e = 0; e < graph_synteny_gp[node1].outgoing_edges[node2].size(); ++e) {
            if (graph_synteny_gp[node1].outgoing_edges[node2].at(e).multi21 > 0) {
                edge_index = e;
                break;
            }
        }
        int min_loc = -1;
        for (int e = 0; e < graph_synteny_gp[node1].outgoing_edges[node2].size(); ++e) {
            if (graph_synteny_gp[node1].outgoing_edges[node2].at(e).original_edge && graph_synteny_gp[node1].outgoing_edges[node2].at(e).multiplicity21.size() > 0) {
                cnt += 1;
                if (min_loc == -1 || graph_synteny_gp[node1].outgoing_edges[node2].at(e).multiplicity21.at(0) < min_loc) {
                    min_loc = graph_synteny_gp[node1].outgoing_edges[node2].at(e).multiplicity21.at(0);
                    edge_index = e;
                }
            }
        }
        assert(edge_index != -1);
        Edge current_edge = graph_synteny[node1].outgoing_edges[node2].at(edge_index);
        if (current_edge.multiplicity >= previous_multi) {
            Block_in_Genome_Path block(current_edge.block_id, this->synteny_blocks[current_edge.block_id].multiplicity, this->synteny_blocks[current_edge.block_id].length, node1, node2, len_genome2, len_genome2 + current_edge.length, true, false, this->synteny_blocks[current_edge.block_id].num_edges);
            blocks_2.emplace_back(block);
        }
        else {
            Block_in_Genome_Path block(current_edge.block_id, this->synteny_blocks[current_edge.block_id].multiplicity, this->synteny_blocks[current_edge.block_id].length, node1, node2, len_genome2 + k, len_genome2 + current_edge.length, false, false, this->synteny_blocks[current_edge.block_id].num_edges);
            blocks_2[blocks_2.size() - 1].include_right = true;
            blocks_2[blocks_2.size() - 1].end += k;
            blocks_2.emplace_back(block);
        }
        previous_multi = current_edge.multiplicity;
        len_genome2 += current_edge.length;
        graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity -= 1;
        graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multi21 -= 1;
        if (min_loc != -1) {
            graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity21.erase(graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity21.begin());
        }
    }
    if (blocks_2.size()) {
        blocks_2[blocks_2.size() - 1].include_right = true;
        blocks_2[blocks_2.size() - 1].end += k;
        merge_blocks(blocks_2);
    }

    // Haplome 3
    int len_genome3 = 0;
    previous_multi = -1;
    std::vector<Block_in_Genome_Path> blocks_3;
    for (int i = 0; i + 1 < this->genome_path_haplome3_strand1.genome_path.size(); ++i) {
        std::string node1 = this->genome_path_haplome3_strand1.genome_path[i], node2 = this->genome_path_haplome3_strand1.genome_path[i + 1];

        // find edge index
        int edge_index = -1;
        int cnt = 0;
        for (int e = 0; e < graph_synteny_gp[node1].outgoing_edges[node2].size(); ++e) {
            if (graph_synteny_gp[node1].outgoing_edges[node2].at(e).multi31 > 0) {
                edge_index = e;
                break;
            }
        }
        int min_loc = -1;
        for (int e = 0; e < graph_synteny_gp[node1].outgoing_edges[node2].size(); ++e) {
            if (graph_synteny_gp[node1].outgoing_edges[node2].at(e).original_edge && graph_synteny_gp[node1].outgoing_edges[node2].at(e).multiplicity31.size() > 0) {
                cnt += 1;
                if (min_loc == -1 || graph_synteny_gp[node1].outgoing_edges[node2].at(e).multiplicity31.at(0) < min_loc) {
                    min_loc = graph_synteny_gp[node1].outgoing_edges[node2].at(e).multiplicity31.at(0);
                    edge_index = e;
                }
            }
        }
        assert(edge_index != -1);
        Edge current_edge = graph_synteny[node1].outgoing_edges[node2].at(edge_index);
        if (current_edge.multiplicity >= previous_multi) {
            Block_in_Genome_Path block(current_edge.block_id, this->synteny_blocks[current_edge.block_id].multiplicity, this->synteny_blocks[current_edge.block_id].length, node1, node2, len_genome3, len_genome3 + current_edge.length, true, false, this->synteny_blocks[current_edge.block_id].num_edges);
            blocks_3.emplace_back(block);
        }
        else {
            Block_in_Genome_Path block(current_edge.block_id, this->synteny_blocks[current_edge.block_id].multiplicity, this->synteny_blocks[current_edge.block_id].length, node1, node2, len_genome3 + k, len_genome3 + current_edge.length, false, false, this->synteny_blocks[current_edge.block_id].num_edges);
            blocks_3[blocks_3.size() - 1].include_right = true;
            blocks_3[blocks_3.size() - 1].end += k;
            blocks_3.emplace_back(block);
        }
        previous_multi = current_edge.multiplicity;
        len_genome3 += current_edge.length;
        graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity -= 1;
        graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multi31 -= 1;
        if (min_loc != -1) {
            graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity31.erase(graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity31.begin());
        }
    }
    if (blocks_3.size()) {
        blocks_3[blocks_3.size() - 1].include_right = true;
        blocks_3[blocks_3.size() - 1].end += k;
        merge_blocks(blocks_3);
    }

    std::unordered_map<int, int> id_map;
    int id_new = 1;
    for (auto&& b : blocks_1) {
        if (id_map.find(b.id) != id_map.end())
            b.id = id_map.at(b.id);
        else {
            id_map[b.id] = id_new;
            id_map[-b.id] = -id_new;
            b.id = id_new;
            id_new += 1;
        }
    }
    for (auto&& b : blocks_2) {
        if (id_map.find(b.id) != id_map.end())
            b.id = id_map.at(b.id);
        else {
            id_map[b.id] = id_new;
            id_map[-b.id] = -id_new;
            b.id = id_new;
            id_new += 1;
        }
    }
    for (auto&& b : blocks_3) {
        if (id_map.find(b.id) != id_map.end())
            b.id = id_map.at(b.id);
        else {
            id_map[b.id] = id_new;
            id_map[-b.id] = -id_new;
            b.id = id_new;
            id_new += 1;
        }
    }

    if (!blocks_1.empty()) {
        // Output synteny blocks
        std::string haplome1_blocks = prefix + ".haplome1.blocks";
        std::ofstream file_haplome1_blocks(haplome1_blocks);
        for (int i = 0; i < blocks_1.size(); ++i) {
            file_haplome1_blocks << blocks_1[i].id << '(' << blocks_1[i].multiplicity << ')' << std::setprecision(1) << std::fixed << 1.0 * (blocks_1[i].end - blocks_1[i].start) / 1000 << '\t' << blocks_1[i].node_s << '\t' << blocks_1[i].node_e << '\t' << blocks_1[i].start << '\t' << blocks_1[i].end << '\t' << blocks_1[i].include_left << '\t' << blocks_1[i].include_right << '\t' << blocks_1[i].is_short_version << std::endl;
        }
        file_haplome1_blocks.close();

    }
    if (!blocks_2.empty()) {
        std::string haplome2_blocks = prefix + ".haplome2.blocks";
        std::ofstream file_haplome2_blocks(haplome2_blocks);
        for (int i = 0; i < blocks_2.size(); ++i) {
            file_haplome2_blocks << blocks_2[i].id << '(' << blocks_2[i].multiplicity << ')' << std::setprecision(1) << std::fixed << 1.0 * (blocks_2[i].end - blocks_2[i].start) / 1000 << '\t' << blocks_2[i].node_s << '\t' << blocks_2[i].node_e << '\t' << blocks_2[i].start << '\t' << blocks_2[i].end << '\t' << blocks_2[i].include_left << '\t' << blocks_2[i].include_right << '\t' << blocks_2[i].is_short_version << std::endl;
        }
        file_haplome2_blocks.close();
    }
    if (!blocks_3.empty()) {
        std::string haplome3_blocks = prefix + ".haplome3.blocks";
        std::ofstream file_haplome3_blocks(haplome3_blocks);
        for (int i = 0; i < blocks_3.size(); ++i) {
            file_haplome3_blocks << blocks_3[i].id << '(' << blocks_3[i].multiplicity << ')' << std::setprecision(1) << std::fixed << 1.0 * (blocks_3[i].end - blocks_3[i].start) / 1000 << '\t' << blocks_3[i].node_s << '\t' << blocks_3[i].node_e << '\t' << blocks_3[i].start << '\t' << blocks_3[i].end << '\t' << blocks_3[i].include_left << '\t' << blocks_3[i].include_right << '\t' << blocks_3[i].is_short_version << std::endl;
        }
        file_haplome3_blocks.close();
    }
}

void Graph::collapse_non_branching_path(std::unordered_map<std::string, Node>& graph, std::string node, std::unordered_set<std::string>& nodes_to_remove) {
    if (graph[node].incoming_edges.size() == 1 && graph[node].outgoing_edges.size() == 1) {
        std::string node_before, node_after;
        Edge edge1, edge2;
        for (auto&& n_b : graph[node].incoming_edges) {
            node_before = n_b.first;
            edge1 = n_b.second.at(0);
        }
        for (auto&& n_a : graph[node].outgoing_edges) {
            node_after = n_a.first;
            edge2 = n_a.second.at(0);
        }
        if (node_before == node || node_after == node)
            return;
        if (graph[node].incoming_edges[node_before].size() == 1 && graph[node].outgoing_edges[node_after].size() == 1) {
            if (edge1.multiplicity > edge2.multiplicity) {
                Edge new_edge(edge1.start_base, edge1.length + edge2.length, edge1.sequence + edge2.sequence.substr(this->k), edge1.multi11, edge1.multi12, edge1.multi21, edge1.multi22, edge1.multi31, edge1.multi32);
                new_edge.multiplicity = edge1.multiplicity;
                graph[node_before].outgoing_edges[node_after].push_back(new_edge);
                graph[node_after].incoming_edges[node_before].push_back(new_edge);
            }
            else {
                Edge new_edge(edge1.start_base, edge1.length + edge2.length, edge1.sequence + edge2.sequence.substr(this->k), edge2.multi11, edge2.multi12, edge2.multi21, edge2.multi22, edge2.multi31, edge2.multi32);
                new_edge.multiplicity = edge2.multiplicity;
                graph[node_before].outgoing_edges[node_after].push_back(new_edge);
                graph[node_after].incoming_edges[node_before].push_back(new_edge);
            }

            graph[node].incoming_edges.erase(node_before);
            graph[node_before].outgoing_edges.erase(node);
            graph[node].outgoing_edges.erase(node_after);
            graph[node_after].incoming_edges.erase(node);
            nodes_to_remove.insert(node);
        }
    }
}

void Graph::write_graph_abnormal(std::unordered_map<std::string, Node>& graph, const std::string& prefix) {
    std::string graph_dot = prefix + ".dot";
    std::ofstream file_dot(graph_dot);
    file_dot << "digraph {\nnodesep = 0.5;\n";
    for (auto&& node : graph) {
        file_dot << node.first << " [style=filled fillcolor=\"white\"]\n";
    }
    for (auto&& node : graph) {
        std::string start_node = node.first;
        for (auto&& edges : node.second.outgoing_edges) {
            std::string end_node = edges.first;
            for (auto&& edge : edges.second) {
                if (this->genome_path_haplome3_strand1.genome_path.size())
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << '+' << edge.multi21 + edge.multi22 << '+' << edge.multi31 + edge.multi32 << ")";
                else if (this->genome_path_haplome2_strand1.genome_path.size())
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << '+' << edge.multi21 + edge.multi22 << ")";
                else
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << ")";
                assert(edge.multiplicity == edge.multi11 + edge.multi12 + edge.multi21 + edge.multi22 + edge.multi31 + edge.multi32);
                if (edge.multiplicity == 1) {
                    if (edge.multi11) {
                        if (this->genome_path_haplome1_strand1.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 1 and strand 1 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G1S1_" << this->genome_path_haplome1_strand1.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi12) {
                        if (this->genome_path_haplome1_strand2.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 1 and strand 2 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G1S2_" << this->genome_path_haplome1_strand2.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi21) {
                        if (this->genome_path_haplome2_strand1.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 2 and strand 1 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G2S1_" << this->genome_path_haplome2_strand1.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi22) {
                        if (this->genome_path_haplome2_strand2.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 2 and strand 2 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G2S2_" << this->genome_path_haplome2_strand2.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi31) {
                        if (this->genome_path_haplome3_strand1.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 3 and strand 1 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G3S1_" << this->genome_path_haplome3_strand1.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi32) {
                        if (this->genome_path_haplome3_strand2.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 3 and strand 2 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G3S2_" << this->genome_path_haplome3_strand2.node2order[start_node][end_node].at(0);
                    }
                }

                assert(edge.sequence.size() - this->k == edge.length);
                assert(edge.sequence.at(this->k) == edge.start_base);

                if (edge.color == 0) {
                    file_dot << "\", color=\"black\"";
                }
                else if (edge.color == 1) {
                    file_dot << "\", color=\"red\"";
                }
                else if (edge.color == 2) {
                    file_dot << "\", color=\"blue\"";
                }
                else if (edge.color == 3) {
                    file_dot << "\", color=\"purple\"";
                }
                else {
                    std::cout << "Error: edge cannot be assigned color (graph output)" << std::endl;
                }
                file_dot << ", penwidth=1]\n";
            }
        }
    }
    file_dot << "}" << std::endl;
    file_dot.close();
}

void Graph::dfs(std::unordered_map<std::string, Node>& graph, std::string node, std::unordered_set<std::string>& component) {
    if (component.find(node) != component.end())
        return;
    component.insert(node);
    for (auto&& n : graph[node].incoming_edges)
        dfs(graph, n.first, component);
    for (auto&& n : graph[node].outgoing_edges)
        dfs(graph, n.first, component);
    return;
}

void Graph::count_duplicaction_graph(const std::string& prefix) {
    std::unordered_map<std::string, Node> duplication_graph = this->graph;
    // remove all 0+1, 1+1 and 1+0 edges
    std::unordered_set<std::string> nodes_to_remove;
    for (auto&& node1 : duplication_graph) {
        std::string n1 = node1.first;
        std::unordered_set<std::string> nodes2_to_remove;
        for (auto&& node2 : node1.second.outgoing_edges) {
            std::string n2 = node2.first;
            for (int i = int(node2.second.size()) - 1; i >= 0; --i) {
                auto& edge = node2.second.at(i);
                if (edge.multi11 + edge.multi12 <= 1 && edge.multi21 + edge.multi22 <= 1) {
                    duplication_graph[n1].outgoing_edges[n2].erase(duplication_graph[n1].outgoing_edges[n2].begin() + i);
                    duplication_graph[n2].incoming_edges[n1].erase(duplication_graph[n2].incoming_edges[n1].begin() + i);
                }
            }
            if (duplication_graph[n1].outgoing_edges[n2].size() == 0)
                nodes2_to_remove.insert(n2);
        }
        for (auto&& n2 : nodes2_to_remove) {
            duplication_graph[n1].outgoing_edges.erase(n2);
            duplication_graph[n2].incoming_edges.erase(n1);
            if (duplication_graph[n1].outgoing_edges.size() == 0 && duplication_graph[n1].incoming_edges.size() == 0) {
                nodes_to_remove.insert(n1);
                nodes_to_remove.insert(reverse_complementary_node(n1));
            }
            if (duplication_graph[n2].outgoing_edges.size() == 0 && duplication_graph[n2].incoming_edges.size() == 0) {
                nodes_to_remove.insert(n2);
                nodes_to_remove.insert(reverse_complementary_node(n2));
            }
        }
    }

    for (auto&& n : nodes_to_remove)
        duplication_graph.erase(n);

    write_graph_abnormal(duplication_graph, prefix + "_duplication1");

    nodes_to_remove.clear();
    for (auto&& n : duplication_graph) {
        if (nodes_to_remove.find(n.first) != nodes_to_remove.end() || nodes_to_remove.find(reverse_complementary_node(n.first)) != nodes_to_remove.end())
            continue;

        collapse_non_branching_path(duplication_graph, n.first, nodes_to_remove);
        collapse_non_branching_path(duplication_graph, reverse_complementary_node(n.first), nodes_to_remove);
    }
    for (auto&& n : nodes_to_remove)
        duplication_graph.erase(n);

    write_graph_abnormal(duplication_graph, prefix + "_duplication2");

    int num_components = 0, circular_components = 0, linear_components = 0, complex_components = 0;
    std::unordered_set<std::string> nodes_scanned;
    std::unordered_set<std::string> current_component;
    for (auto&& n : duplication_graph) {
        if (nodes_scanned.find(n.first) != nodes_scanned.end())
            continue;
        else {
            current_component.clear();
            dfs(duplication_graph, n.first, current_component);
            if (current_component.size() == 1) {
                circular_components += 1;
                // std::cout << "Circular: ";
            }
            else if (current_component.size() == 2) {
                auto it = current_component.begin();
                std::string node1 = *it;
                std::string node2 = *(++it);

                int num_edges = 0;
                for (auto&& n : duplication_graph[node1].outgoing_edges) {
                    for (auto&& e : n.second)
                        num_edges += 1;
                }
                for (auto&& n : duplication_graph[node2].outgoing_edges) {
                    for (auto&& e : n.second)
                        num_edges += 1;
                }

                if (num_edges == 1) {
                    linear_components += 1;
                    // std::cout << "Linear: ";
                }
                else {
                    complex_components += 1;
                    // std::cout << "Complex: ";
                }

            }
            else {
                complex_components += 1;
                // std::cout << "Complex: ";
            }
            for (auto&& c : current_component) {
                nodes_scanned.insert(c);
                // std::cout << c << " ";
            }
            // std::cout << std::endl;
            num_components += 1;
        }
    }

    std::cout << "\nCircular components: " << circular_components << std::endl;
    std::cout << "Linear components: " << linear_components << std::endl;
    std::cout << "Complex components: " << complex_components << std::endl;
    std::cout << "All components: " << num_components << std::endl;
}

void Graph::count_imbalanced_graph(const std::string& prefix) {
    std::unordered_map<std::string, Node> imbalanced_graph = this->graph;
    for (auto&& node1 : imbalanced_graph) {
        std::string n1 = node1.first;
        std::unordered_set<std::string> nodes2_to_remove;
        for (auto&& node2 : node1.second.outgoing_edges) {
            std::string n2 = node2.first;
            for (int i = int(node2.second.size()) - 1; i >= 0; --i) {
                auto& edge = node2.second.at(i);
                if (edge.multi11 + edge.multi12 == edge.multi21 + edge.multi22) {
                    if (this->genome_path_haplome3_strand1.genome_path.size() && edge.multi11 + edge.multi12 != edge.multi31 + edge.multi32)
                        continue;
                    imbalanced_graph[n1].outgoing_edges[n2].at(i).block_id = 1;
                    imbalanced_graph[n2].incoming_edges[n1].at(i).block_id = 1;
                }
            }
        }
    }
    write_graph(imbalanced_graph, prefix + "_imbalanced");

    // remove all 0+1, 1+1 and 1+0 edges
    std::unordered_set<std::string> nodes_to_remove;
    for (auto&& node1 : imbalanced_graph) {
        std::string n1 = node1.first;
        std::unordered_set<std::string> nodes2_to_remove;
        for (auto&& node2 : node1.second.outgoing_edges) {
            std::string n2 = node2.first;
            for (int i = int(node2.second.size()) - 1; i >= 0; --i) {
                auto& edge = node2.second.at(i);
                if (edge.multi11 + edge.multi12 == edge.multi21 + edge.multi22 && edge.multi11 + edge.multi12 == edge.multi31 + edge.multi32) {
                    imbalanced_graph[n1].outgoing_edges[n2].erase(imbalanced_graph[n1].outgoing_edges[n2].begin() + i);
                    imbalanced_graph[n2].incoming_edges[n1].erase(imbalanced_graph[n2].incoming_edges[n1].begin() + i);
                }
            }
            if (imbalanced_graph[n1].outgoing_edges[n2].size() == 0)
                nodes2_to_remove.insert(n2);
        }
        for (auto&& n2 : nodes2_to_remove) {
            imbalanced_graph[n1].outgoing_edges.erase(n2);
            imbalanced_graph[n2].incoming_edges.erase(n1);
            if (imbalanced_graph[n1].outgoing_edges.size() == 0 && imbalanced_graph[n1].incoming_edges.size() == 0) {
                nodes_to_remove.insert(n1);
                nodes_to_remove.insert(reverse_complementary_node(n1));
            }
            if (imbalanced_graph[n2].outgoing_edges.size() == 0 && imbalanced_graph[n2].incoming_edges.size() == 0) {
                nodes_to_remove.insert(n2);
                nodes_to_remove.insert(reverse_complementary_node(n2));
            }
        }
    }

    for (auto&& n : nodes_to_remove)
        imbalanced_graph.erase(n);

    write_graph_abnormal(imbalanced_graph, prefix + "_imbalanced1");

    nodes_to_remove.clear();
    for (auto&& n : imbalanced_graph) {
        if (nodes_to_remove.find(n.first) != nodes_to_remove.end() || nodes_to_remove.find(reverse_complementary_node(n.first)) != nodes_to_remove.end())
            continue;

        collapse_non_branching_path(imbalanced_graph, n.first, nodes_to_remove);
        collapse_non_branching_path(imbalanced_graph, reverse_complementary_node(n.first), nodes_to_remove);
    }
    for (auto&& n : nodes_to_remove)
        imbalanced_graph.erase(n);

    write_graph_abnormal(imbalanced_graph, prefix + "_imbalanced2");

    // int num_components = 0, circular_components = 0, linear_components = 0, complex_components = 0;
    // std::unordered_set<std::string> nodes_scanned;
    // std::unordered_set<std::string> current_component;
    // for (auto&& n : imbalanced_graph) {
    //     if (nodes_scanned.find(n.first) != nodes_scanned.end())
    //         continue;
    //     else {
    //         current_component.clear();
    //         dfs(imbalanced_graph, n.first, current_component);
    //         if (current_component.size() == 1) {
    //             circular_components += 1;
    //             // std::cout << "Circular: ";
    //         }
    //         else if (current_component.size() == 2) {
    //             auto it = current_component.begin();
    //             std::string node1 = *it;
    //             std::string node2 = *(++it);

    //             int num_edges = 0;
    //             for (auto&& n : imbalanced_graph[node1].outgoing_edges) {
    //                 for (auto&& e : n.second)
    //                     num_edges += 1;
    //             }
    //             for (auto&& n : imbalanced_graph[node2].outgoing_edges) {
    //                 for (auto&& e : n.second)
    //                     num_edges += 1;
    //             }

    //             if (num_edges == 1) {
    //                 linear_components += 1;
    //                 // std::cout << "Linear: ";
    //             }
    //             else {
    //                 complex_components += 1;
    //                 // std::cout << "Complex: ";
    //             }

    //         }
    //         else {
    //             complex_components += 1;
    //             // std::cout << "Complex: ";
    //         }
    //         for (auto&& c : current_component) {
    //             nodes_scanned.insert(c);
    //             // std::cout << c << " ";
    //         }
    //         // std::cout << std::endl;
    //         num_components += 1;
    //     }
    // }

    // std::cout << "\nCircular components: " << circular_components << std::endl;
    // std::cout << "Linear components: " << linear_components << std::endl;
    // std::cout << "Complex components: " << complex_components << std::endl;
    // std::cout << "All components: " << num_components << '\n' << std::endl;
}

void Graph::write_graph(const std::string& prefix, int thick, int min_synteny_length) {
    std::string graph_dot = prefix + ".dot";
    std::string graph_fasta = prefix + ".fasta";
    std::string haplome1_strand1 = prefix + ".haplome1.strand1.path";
    std::string haplome1_strand2 = prefix + ".haplome1.strand2.path";
    std::string haplome2_strand1 = prefix + ".haplome2.strand1.path";
    std::string haplome2_strand2 = prefix + ".haplome2.strand2.path";
    std::string haplome3_strand1 = prefix + ".haplome3.strand1.path";
    std::string haplome3_strand2 = prefix + ".haplome3.strand2.path";
    std::string haplome1_sequence = prefix + ".haplome1.fa";
    std::string haplome2_sequence = prefix + ".haplome2.fa";
    std::string haplome3_sequence = prefix + ".haplome3.fa";
    std::string haplome1_locations = prefix + ".haplome1.locus";
    std::string haplome1_blocks = prefix + ".haplome1.blocks";
    std::string haplome2_locations = prefix + ".haplome2.locus";
    std::string haplome2_blocks = prefix + ".haplome2.blocks";
    std::string haplome3_locations = prefix + ".haplome3.locus";
    std::string haplome3_blocks = prefix + ".haplome3.blocks";

    std::cout << "\nWrite graph " << graph_dot << ", fasta " << graph_fasta << std::endl;
    std::ofstream file_dot(graph_dot);
    std::ofstream file_fasta(graph_fasta);

    std::ofstream file_haplome1_strand1(haplome1_strand1);
    std::ofstream file_haplome1_strand2(haplome1_strand2);
    std::ofstream file_haplome1_sequence(haplome1_sequence);
    std::ofstream file_haplome1_locus(haplome1_locations);
    file_haplome1_sequence << ">contig_1" << std::endl;
    this->genome_path_haplome1_strand1.construct_edge_orders();
    this->genome_path_haplome1_strand2.construct_edge_orders();

    auto graph_gp = graph;

    // Haplome 1
    int len_genome1 = 0;
    int previous_multi = -1;
    std::vector<Block_in_Genome_Path> blocks_1;
    for (int i = 0; i + 1 < this->genome_path_haplome1_strand1.genome_path.size(); ++i) {
        file_haplome1_locus << len_genome1 << ",";
        file_haplome1_strand1 << this->genome_path_haplome1_strand1.genome_path[i] << "->" << this->genome_path_haplome1_strand1.genome_path[i + 1] << "\n";
        std::string node1 = this->genome_path_haplome1_strand1.genome_path[i], node2 = this->genome_path_haplome1_strand1.genome_path[i + 1];

        // find edge index
        int edge_index = -1;
        int cnt = 0;
        for (int e = 0; e < graph_gp[node1].outgoing_edges[node2].size(); ++e) {
            if (graph_gp[node1].outgoing_edges[node2].at(e).multi11 > 0) {
                edge_index = e;
                break;
            }
        }
        int min_loc = -1;
        for (int e = 0; e < graph_gp[node1].outgoing_edges[node2].size(); ++e) {
            if (graph_gp[node1].outgoing_edges[node2].at(e).original_edge && graph_gp[node1].outgoing_edges[node2].at(e).multiplicity11.size() > 0) {
                cnt += 1;
                if (min_loc == -1 || graph_gp[node1].outgoing_edges[node2].at(e).multiplicity11.at(0) < min_loc) {
                    min_loc = graph_gp[node1].outgoing_edges[node2].at(e).multiplicity11.at(0);
                    edge_index = e;
                }
            }
        }
        assert(edge_index != -1);

        Edge current_edge = graph[node1].outgoing_edges[node2].at(edge_index);
        if (i == 0)
            file_haplome1_sequence << current_edge.sequence;
        else
            file_haplome1_sequence << current_edge.sequence.substr(k);
        len_genome1 += current_edge.length;

        graph_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity -= 1;
        graph_gp[node1].outgoing_edges[node2].at(edge_index).multi11 -= 1;
        if (min_loc != -1) {
            graph_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity11.erase(graph_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity11.begin());
        }
    }
    for (int i = 0; i + 1 < this->genome_path_haplome1_strand2.genome_path.size(); ++i) {
        file_haplome1_strand2 << this->genome_path_haplome1_strand2.genome_path[i] << "->" << this->genome_path_haplome1_strand2.genome_path[i + 1] << "\n";
    }
    file_haplome1_strand1.close();
    file_haplome1_strand2.close();
    file_haplome1_sequence.close();
    file_haplome1_locus << len_genome1 + k << std::endl;
    file_haplome1_locus.close();

    // Haplome 2
    int len_genome2 = 0;
    if (this->genome_path_haplome2_strand1.genome_path.size()) {
        std::ofstream file_haplome2_strand1(haplome2_strand1);
        std::ofstream file_haplome2_strand2(haplome2_strand2);
        std::ofstream file_haplome2_sequence(haplome2_sequence);
        std::ofstream file_haplome2_locus(haplome2_locations);
        file_haplome2_sequence << ">contig_1" << std::endl;
        this->genome_path_haplome2_strand1.construct_edge_orders();
        this->genome_path_haplome2_strand2.construct_edge_orders();
        previous_multi = -1;
        std::vector<Block_in_Genome_Path> blocks_2;
        for (int i = 0; i + 1 < this->genome_path_haplome2_strand1.genome_path.size(); ++i) {
            file_haplome2_locus << len_genome2 << ",";
            file_haplome2_strand1 << this->genome_path_haplome2_strand1.genome_path[i] << "->" << this->genome_path_haplome2_strand1.genome_path[i + 1] << "\n";
            std::string node1 = this->genome_path_haplome2_strand1.genome_path[i], node2 = this->genome_path_haplome2_strand1.genome_path[i + 1];

            // find edge index
            int edge_index = -1;
            int cnt = 0;
            for (int e = 0; e < graph_gp[node1].outgoing_edges[node2].size(); ++e) {
                if (graph_gp[node1].outgoing_edges[node2].at(e).multi21 > 0) {
                    edge_index = e;
                    break;
                }
            }
            int min_loc = -1;
            for (int e = 0; e < graph_gp[node1].outgoing_edges[node2].size(); ++e) {
                if (graph_gp[node1].outgoing_edges[node2].at(e).original_edge && graph_gp[node1].outgoing_edges[node2].at(e).multiplicity21.size() > 0) {
                    cnt += 1;
                    if (min_loc == -1 || graph_gp[node1].outgoing_edges[node2].at(e).multiplicity21.at(0) < min_loc) {
                        min_loc = graph_gp[node1].outgoing_edges[node2].at(e).multiplicity21.at(0);
                        edge_index = e;
                    }
                }
            }
            assert(edge_index != -1);

            Edge current_edge = graph[node1].outgoing_edges[node2].at(edge_index);
            if (i == 0)
                file_haplome2_sequence << current_edge.sequence;
            else
                file_haplome2_sequence << current_edge.sequence.substr(k);
            len_genome2 += current_edge.length;

            graph_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity -= 1;
            graph_gp[node1].outgoing_edges[node2].at(edge_index).multi21 -= 1;
            if (min_loc != -1) {
                graph_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity21.erase(graph_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity21.begin());
            }
        }
        for (int i = 0; i + 1 < this->genome_path_haplome2_strand2.genome_path.size(); ++i) {
            file_haplome2_strand2 << this->genome_path_haplome2_strand2.genome_path[i] << "->" << this->genome_path_haplome2_strand2.genome_path[i + 1] << "\n";
        }
        file_haplome2_strand1.close();
        file_haplome2_strand2.close();
        file_haplome2_sequence.close();
        file_haplome2_locus << len_genome2 + k << std::endl;
        file_haplome2_locus.close();
    }

    // Haplome 3
    int len_genome3 = 0;
    if (this->genome_path_haplome3_strand1.genome_path.size()) {
        std::ofstream file_haplome3_strand1(haplome3_strand1);
        std::ofstream file_haplome3_strand2(haplome3_strand2);
        std::ofstream file_haplome3_sequence(haplome3_sequence);
        std::ofstream file_haplome3_locus(haplome3_locations);
        file_haplome3_sequence << ">contig_1" << std::endl;
        this->genome_path_haplome3_strand1.construct_edge_orders();
        this->genome_path_haplome3_strand2.construct_edge_orders();
        previous_multi = -1;
        std::vector<Block_in_Genome_Path> blocks_3;
        for (int i = 0; i + 1 < this->genome_path_haplome3_strand1.genome_path.size(); ++i) {
            file_haplome3_locus << len_genome3 << ",";
            file_haplome3_strand1 << this->genome_path_haplome3_strand1.genome_path[i] << "->" << this->genome_path_haplome3_strand1.genome_path[i + 1] << "\n";
            std::string node1 = this->genome_path_haplome3_strand1.genome_path[i], node2 = this->genome_path_haplome3_strand1.genome_path[i + 1];

            // find edge index
            int edge_index = -1;
            int cnt = 0;
            for (int e = 0; e < graph_gp[node1].outgoing_edges[node2].size(); ++e) {
                if (graph_gp[node1].outgoing_edges[node2].at(e).multi31 > 0) {
                    edge_index = e;
                    break;
                }
            }
            int min_loc = -1;
            for (int e = 0; e < graph_gp[node1].outgoing_edges[node2].size(); ++e) {
                if (graph_gp[node1].outgoing_edges[node2].at(e).original_edge && graph_gp[node1].outgoing_edges[node2].at(e).multiplicity31.size() > 0) {
                    cnt += 1;
                    if (min_loc == -1 || graph_gp[node1].outgoing_edges[node2].at(e).multiplicity31.at(0) < min_loc) {
                        min_loc = graph_gp[node1].outgoing_edges[node2].at(e).multiplicity31.at(0);
                        edge_index = e;
                    }
                }
            }
            assert(edge_index != -1);

            Edge current_edge = graph[node1].outgoing_edges[node2].at(edge_index);
            if (i == 0)
                file_haplome3_sequence << current_edge.sequence;
            else
                file_haplome3_sequence << current_edge.sequence.substr(k);
            len_genome3 += current_edge.length;

            graph_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity -= 1;
            graph_gp[node1].outgoing_edges[node2].at(edge_index).multi31 -= 1;
            if (min_loc != -1) {
                graph_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity31.erase(graph_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity31.begin());
            }
        }
        for (int i = 0; i + 1 < this->genome_path_haplome3_strand2.genome_path.size(); ++i) {
            file_haplome3_strand2 << this->genome_path_haplome3_strand2.genome_path[i] << "->" << this->genome_path_haplome3_strand2.genome_path[i + 1] << "\n";
        }
        file_haplome3_strand1.close();
        file_haplome3_strand2.close();
        file_haplome3_sequence.close();
        file_haplome3_locus << len_genome3 + k << std::endl;
        file_haplome3_locus.close();
    }

    unsigned cnt_long = 0, cnt_short = 0, length_long = 0, length_short = 0, cnt_red = 0, length_red = 0, cnt_black = 0, length_black = 0, cnt_blue = 0, length_blue = 0;
    unsigned cnt_multi1 = 0, len_multi1 = 0, cnt_multi2 = 0, len_multi2 = 0, cnt_multi3 = 0, len_multi3 = 0, cnt_more_3 = 0, len_more_3 = 0;
    unsigned cnt_non_short = 0, cnt_red_non_short = 0, cnt_black_non_short = 0, cnt_blue_non_short = 0, cnt_multi1_non_short = 0, cnt_multi2_non_short = 0, cnt_multi3_non_short = 0, cnt_more_3_non_short = 0;
    unsigned cnt_imbalanced = 0, cnt_imbalanced_non_short = 0, len_imbalanced = 0;
    unsigned total_len = 0;
    file_dot << "digraph {\nnodesep = 0.5;\n";
    for (auto&& node : graph) {
        file_dot << node.first << " [style=filled fillcolor=\"white\"]\n";
    }
    for (auto&& node : graph) {
        std::string start_node = node.first;
        for (auto&& edges : node.second.outgoing_edges) {
            std::string end_node = edges.first;
            for (auto&& edge : edges.second) {
                if (this->genome_path_haplome3_strand1.genome_path.size())
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << '+' << edge.multi21 + edge.multi22 << '+' << edge.multi31 + edge.multi32 << ")";
                else if (this->genome_path_haplome2_strand1.genome_path.size())
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << '+' << edge.multi21 + edge.multi22 << ")";
                else
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << ")";
                assert(edge.multiplicity == edge.multi11 + edge.multi12 + edge.multi21 + edge.multi22 + edge.multi31 + edge.multi32);
                if (edge.multiplicity == 1) {
                    if (edge.multi11) {
                        if (this->genome_path_haplome1_strand1.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 1 and strand 1 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G1S1_" << this->genome_path_haplome1_strand1.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi12) {
                        if (this->genome_path_haplome1_strand2.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 1 and strand 2 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G1S2_" << this->genome_path_haplome1_strand2.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi21) {
                        if (this->genome_path_haplome2_strand1.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 2 and strand 1 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G2S1_" << this->genome_path_haplome2_strand1.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi22) {
                        if (this->genome_path_haplome2_strand2.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 2 and strand 2 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G2S2_" << this->genome_path_haplome2_strand2.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi31) {
                        if (this->genome_path_haplome3_strand1.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 3 and strand 1 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G3S1_" << this->genome_path_haplome3_strand1.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi32) {
                        if (this->genome_path_haplome3_strand2.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 3 and strand 2 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G3S2_" << this->genome_path_haplome3_strand2.node2order[start_node][end_node].at(0);
                    }
                }

                std::string color_out;
                if (edge.color == 0) {
                    cnt_black += 1;
                    length_black += edge.length;
                    color_out = "black";
                    // file_dot << "\", color=\"black\"";
                    if (edge.length > 300)
                        cnt_black_non_short += 1;
                }
                else if (edge.color == 1) {
                    cnt_red += 1;
                    length_red += edge.length;
                    color_out = "red";
                    // file_dot << "\", color=\"red\"";
                    if (edge.length > 300)
                        cnt_red_non_short += 1;
                }
                else if (edge.color == 2) {
                    cnt_blue += 1;
                    length_blue += edge.length;
                    color_out = "blue";
                    // file_dot << "\", color=\"blue\"";
                    if (edge.length > 300)
                        cnt_blue_non_short += 1;
                }
                else if (edge.color == 3) {
                    // cnt_blue += 1;
                    // length_blue += edge.length;
                    color_out = "purple";
                    // file_dot << "\", color=\"blue\"";
                    // if (edge.length > 300)
                    //     cnt_blue_non_short += 1;
                }
                else {
                    std::cout << "Error: edge cannot be assigned color (graph output)" << std::endl;
                }

                if (edge.block_id != 0)
                    file_dot << "\", color=\"azure3\"";
                else
                    file_dot << "\", color=\"" + color_out + "\"";

                total_len += edge.length * edge.multiplicity;

                if (edge.multi11 + edge.multi12 != edge.multi21 + edge.multi22) {
                    cnt_imbalanced += 1;
                    len_imbalanced += edge.length;
                    if (edge.length > 300)
                        cnt_imbalanced_non_short += 1;
                }

                if (edge.length >= thick) {
                    cnt_long += 1;
                    length_long += edge.length;
                    file_dot << ", penwidth=5]\n";
                }
                else {
                    cnt_short += 1;
                    length_short += edge.length;
                    file_dot << ", penwidth=1]\n";
                }

                if (edge.length > 300)
                    cnt_non_short += 1;

                if (edge.multiplicity == 1) {
                    cnt_multi1 += 1;
                    len_multi1 += edge.length;
                    if (edge.length > 300)
                        cnt_multi1_non_short += 1;
                }
                else if (edge.multiplicity == 2) {
                    cnt_multi2 += 1;
                    len_multi2 += edge.length;
                    if (edge.length > 300)
                        cnt_multi2_non_short += 1;
                }
                else if (edge.multiplicity == 3) {
                    cnt_multi3 += 1;
                    len_multi3 += edge.length;
                    if (edge.length > 300)
                        cnt_multi3_non_short += 1;
                }
                else {
                    cnt_more_3 += 1;
                    len_more_3 += edge.length;
                    if (edge.length > 300)
                        cnt_more_3_non_short += 1;
                }
                assert(edge.sequence.size() - this->k == edge.length);
                assert(edge.sequence.at(this->k) == edge.start_base);
                file_fasta << ">" << start_node << "_" << end_node << "_" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << '+' << edge.multi21 + edge.multi22 << ")" << edge.color << "\n";
                file_fasta << edge.sequence << "\n";
            }
        }
    }
    file_dot << "}" << std::endl;
    file_dot.close();
    file_fasta.close();

    std::cout << "Total number of nodes: " << this->get_num_nodes() << std::endl;
    std::cout << "Total number of edges (non-short edges): " << cnt_long + cnt_short << " (" << cnt_non_short << ")" << std::endl;
    if (!this->genome_path_haplome2_strand1.genome_path.size()) {
        std::cout << "Length of red edges: " << length_red << std::endl;
    }
    else if (!this->genome_path_haplome1_strand1.genome_path.size()) {
        std::cout << "Length of blue edges: " << length_blue << std::endl;

    }
    else {
        std::cout << "Bicolored edges: " << cnt_black << std::endl;
        std::cout << "Length of bicolored edges: " << length_black << std::endl;
        std::cout << "Red edges: " << cnt_red << std::endl;
        std::cout << "Length of red edges: " << length_red << std::endl;
    }
    std::cout << "Multi 1 (non-short edges): " << cnt_multi1 << " (" << cnt_multi1_non_short << ")" << std::endl;
    std::cout << "Length: " << len_multi1 << std::endl;
    std::cout << "Multi 2 (non-short edges): " << cnt_multi2 << " (" << cnt_multi2_non_short << ")" << std::endl;
    std::cout << "Length: " << len_multi2 << std::endl;
    std::cout << "Multi 3 (non-short edges): " << cnt_multi3 << " (" << cnt_multi3_non_short << ")" << std::endl;
    std::cout << "Length: " << len_multi3 << std::endl;
    std::cout << "Multi >3 (non-short edges): " << cnt_more_3 << " (" << cnt_more_3_non_short << ")" << std::endl;
    std::cout << "Length: " << len_more_3 << std::endl;

    if (this->genome_path_haplome1_strand1.genome_path.size()) {
        std::cout << "Genome 1 edges: " << this->genome_path_haplome1_strand1.genome_path.size() - 1 << std::endl;
        std::cout << "Genome 1 length: " << len_genome1 + k << std::endl;
    }
    if (this->genome_path_haplome2_strand1.genome_path.size()) {
        std::cout << "Genome 2 edges: " << this->genome_path_haplome2_strand1.genome_path.size() - 1 << std::endl;
        std::cout << "Genome 2 length: " << len_genome2 + k << std::endl;
    }
    if (this->genome_path_haplome3_strand1.genome_path.size()) {
        std::cout << "Genome 3 edges: " << this->genome_path_haplome3_strand1.genome_path.size() - 1 << std::endl;
        std::cout << "Genome 3 length: " << len_genome3 + k << std::endl;
    }

    std::cout << "Imbalanced edges (non-short edges): " << cnt_imbalanced << " (" << cnt_imbalanced_non_short << ")" << std::endl;
    std::cout << "Length: " << len_imbalanced << std::endl;
}

void Graph::write_graph(std::unordered_map<std::string, Node>& graph_modified, const std::string& prefix, int thick, int min_synteny_length) {
    std::string graph_modified_dot = prefix + ".dot";
    std::string graph_modified_fasta = prefix + ".fasta";
    std::string haplome1_strand1 = prefix + ".haplome1.strand1.path";
    std::string haplome1_strand2 = prefix + ".haplome1.strand2.path";
    std::string haplome2_strand1 = prefix + ".haplome2.strand1.path";
    std::string haplome2_strand2 = prefix + ".haplome2.strand2.path";
    std::string haplome3_strand1 = prefix + ".haplome3.strand1.path";
    std::string haplome3_strand2 = prefix + ".haplome3.strand2.path";
    std::string haplome1_sequence = prefix + ".haplome1.fa";
    std::string haplome2_sequence = prefix + ".haplome2.fa";
    std::string haplome3_sequence = prefix + ".haplome3.fa";
    std::string haplome1_locations = prefix + ".haplome1.locus";
    std::string haplome1_blocks = prefix + ".haplome1.blocks";
    std::string haplome2_locations = prefix + ".haplome2.locus";
    std::string haplome2_blocks = prefix + ".haplome2.blocks";
    std::string haplome3_locations = prefix + ".haplome3.locus";
    std::string haplome3_blocks = prefix + ".haplome3.blocks";

    std::cout << "\nWrite graph " << graph_modified_dot << ", fasta " << graph_modified_fasta << std::endl;
    std::ofstream file_dot(graph_modified_dot);
    std::ofstream file_fasta(graph_modified_fasta);

    std::ofstream file_haplome1_strand1(haplome1_strand1);
    std::ofstream file_haplome1_strand2(haplome1_strand2);
    std::ofstream file_haplome1_sequence(haplome1_sequence);
    std::ofstream file_haplome1_locus(haplome1_locations);
    file_haplome1_sequence << ">contig_1" << std::endl;
    this->genome_path_haplome1_strand1.construct_edge_orders();
    this->genome_path_haplome1_strand2.construct_edge_orders();

    auto graph_modified_gp = graph_modified;

    // Haplome 1
    int len_genome1 = 0;
    int previous_multi = -1;
    std::vector<Block_in_Genome_Path> blocks_1;
    for (int i = 0; i + 1 < this->genome_path_haplome1_strand1.genome_path.size(); ++i) {
        file_haplome1_locus << len_genome1 << ",";
        file_haplome1_strand1 << this->genome_path_haplome1_strand1.genome_path[i] << "->" << this->genome_path_haplome1_strand1.genome_path[i + 1] << "\n";
        std::string node1 = this->genome_path_haplome1_strand1.genome_path[i], node2 = this->genome_path_haplome1_strand1.genome_path[i + 1];

        // find edge index
        int edge_index = -1;
        int cnt = 0;
        for (int e = 0; e < graph_modified_gp[node1].outgoing_edges[node2].size(); ++e) {
            if (graph_modified_gp[node1].outgoing_edges[node2].at(e).multi11 > 0) {
                edge_index = e;
                break;
            }
        }
        int min_loc = -1;
        for (int e = 0; e < graph_modified_gp[node1].outgoing_edges[node2].size(); ++e) {
            if (graph_modified_gp[node1].outgoing_edges[node2].at(e).original_edge && graph_modified_gp[node1].outgoing_edges[node2].at(e).multiplicity11.size() > 0) {
                cnt += 1;
                if (min_loc == -1 || graph_modified_gp[node1].outgoing_edges[node2].at(e).multiplicity11.at(0) < min_loc) {
                    min_loc = graph_modified_gp[node1].outgoing_edges[node2].at(e).multiplicity11.at(0);
                    edge_index = e;
                }
            }
        }
        assert(edge_index != -1);

        Edge current_edge = graph_modified[node1].outgoing_edges[node2].at(edge_index);
        if (i == 0)
            file_haplome1_sequence << current_edge.sequence;
        else
            file_haplome1_sequence << current_edge.sequence.substr(k);
        len_genome1 += current_edge.length;

        graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity -= 1;
        graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multi11 -= 1;
        if (min_loc != -1) {
            graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity11.erase(graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity11.begin());
        }
    }
    for (int i = 0; i + 1 < this->genome_path_haplome1_strand2.genome_path.size(); ++i) {
        file_haplome1_strand2 << this->genome_path_haplome1_strand2.genome_path[i] << "->" << this->genome_path_haplome1_strand2.genome_path[i + 1] << "\n";
    }
    file_haplome1_strand1.close();
    file_haplome1_strand2.close();
    file_haplome1_sequence.close();
    file_haplome1_locus << len_genome1 + k << std::endl;
    file_haplome1_locus.close();

    // Haplome 2
    int len_genome2 = 0;
    if (this->genome_path_haplome2_strand1.genome_path.size()) {
        std::ofstream file_haplome2_strand1(haplome2_strand1);
        std::ofstream file_haplome2_strand2(haplome2_strand2);
        std::ofstream file_haplome2_sequence(haplome2_sequence);
        std::ofstream file_haplome2_locus(haplome2_locations);
        file_haplome2_sequence << ">contig_1" << std::endl;
        this->genome_path_haplome2_strand1.construct_edge_orders();
        this->genome_path_haplome2_strand2.construct_edge_orders();
        previous_multi = -1;
        std::vector<Block_in_Genome_Path> blocks_2;
        for (int i = 0; i + 1 < this->genome_path_haplome2_strand1.genome_path.size(); ++i) {
            file_haplome2_locus << len_genome2 << ",";
            file_haplome2_strand1 << this->genome_path_haplome2_strand1.genome_path[i] << "->" << this->genome_path_haplome2_strand1.genome_path[i + 1] << "\n";
            std::string node1 = this->genome_path_haplome2_strand1.genome_path[i], node2 = this->genome_path_haplome2_strand1.genome_path[i + 1];

            // find edge index
            int edge_index = -1;
            int cnt = 0;
            for (int e = 0; e < graph_modified_gp[node1].outgoing_edges[node2].size(); ++e) {
                if (graph_modified_gp[node1].outgoing_edges[node2].at(e).multi21 > 0) {
                    edge_index = e;
                    break;
                }
            }
            int min_loc = -1;
            for (int e = 0; e < graph_modified_gp[node1].outgoing_edges[node2].size(); ++e) {
                if (graph_modified_gp[node1].outgoing_edges[node2].at(e).original_edge && graph_modified_gp[node1].outgoing_edges[node2].at(e).multiplicity21.size() > 0) {
                    cnt += 1;
                    if (min_loc == -1 || graph_modified_gp[node1].outgoing_edges[node2].at(e).multiplicity21.at(0) < min_loc) {
                        min_loc = graph_modified_gp[node1].outgoing_edges[node2].at(e).multiplicity21.at(0);
                        edge_index = e;
                    }
                }
            }
            assert(edge_index != -1);

            Edge current_edge = graph_modified[node1].outgoing_edges[node2].at(edge_index);
            if (i == 0)
                file_haplome2_sequence << current_edge.sequence;
            else
                file_haplome2_sequence << current_edge.sequence.substr(k);
            len_genome2 += current_edge.length;

            graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity -= 1;
            graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multi21 -= 1;
            if (min_loc != -1) {
                graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity21.erase(graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity21.begin());
            }
        }
        for (int i = 0; i + 1 < this->genome_path_haplome2_strand2.genome_path.size(); ++i) {
            file_haplome2_strand2 << this->genome_path_haplome2_strand2.genome_path[i] << "->" << this->genome_path_haplome2_strand2.genome_path[i + 1] << "\n";
        }
        file_haplome2_strand1.close();
        file_haplome2_strand2.close();
        file_haplome2_sequence.close();
        file_haplome2_locus << len_genome2 + k << std::endl;
        file_haplome2_locus.close();
    }

    // Haplome 3
    int len_genome3 = 0;
    if (this->genome_path_haplome3_strand1.genome_path.size()) {
        std::ofstream file_haplome3_strand1(haplome3_strand1);
        std::ofstream file_haplome3_strand2(haplome3_strand2);
        std::ofstream file_haplome3_sequence(haplome3_sequence);
        std::ofstream file_haplome3_locus(haplome3_locations);
        file_haplome3_sequence << ">contig_1" << std::endl;
        this->genome_path_haplome3_strand1.construct_edge_orders();
        this->genome_path_haplome3_strand2.construct_edge_orders();
        previous_multi = -1;
        std::vector<Block_in_Genome_Path> blocks_3;
        for (int i = 0; i + 1 < this->genome_path_haplome3_strand1.genome_path.size(); ++i) {
            file_haplome3_locus << len_genome3 << ",";
            file_haplome3_strand1 << this->genome_path_haplome3_strand1.genome_path[i] << "->" << this->genome_path_haplome3_strand1.genome_path[i + 1] << "\n";
            std::string node1 = this->genome_path_haplome3_strand1.genome_path[i], node2 = this->genome_path_haplome3_strand1.genome_path[i + 1];

            // find edge index
            int edge_index = -1;
            int cnt = 0;
            for (int e = 0; e < graph_modified_gp[node1].outgoing_edges[node2].size(); ++e) {
                if (graph_modified_gp[node1].outgoing_edges[node2].at(e).multi31 > 0) {
                    edge_index = e;
                    break;
                }
            }
            int min_loc = -1;
            for (int e = 0; e < graph_modified_gp[node1].outgoing_edges[node2].size(); ++e) {
                if (graph_modified_gp[node1].outgoing_edges[node2].at(e).original_edge && graph_modified_gp[node1].outgoing_edges[node2].at(e).multiplicity31.size() > 0) {
                    cnt += 1;
                    if (min_loc == -1 || graph_modified_gp[node1].outgoing_edges[node2].at(e).multiplicity31.at(0) < min_loc) {
                        min_loc = graph_modified_gp[node1].outgoing_edges[node2].at(e).multiplicity31.at(0);
                        edge_index = e;
                    }
                }
            }
            assert(edge_index != -1);

            Edge current_edge = graph_modified[node1].outgoing_edges[node2].at(edge_index);
            if (i == 0)
                file_haplome3_sequence << current_edge.sequence;
            else
                file_haplome3_sequence << current_edge.sequence.substr(k);
            len_genome3 += current_edge.length;

            graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity -= 1;
            graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multi31 -= 1;
            if (min_loc != -1) {
                graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity31.erase(graph_modified_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity31.begin());
            }
        }
        for (int i = 0; i + 1 < this->genome_path_haplome3_strand2.genome_path.size(); ++i) {
            file_haplome3_strand2 << this->genome_path_haplome3_strand2.genome_path[i] << "->" << this->genome_path_haplome3_strand2.genome_path[i + 1] << "\n";
        }
        file_haplome3_strand1.close();
        file_haplome3_strand2.close();
        file_haplome3_sequence.close();
        file_haplome3_locus << len_genome3 + k << std::endl;
        file_haplome3_locus.close();
    }

    unsigned cnt_long = 0, cnt_short = 0, length_long = 0, length_short = 0, cnt_red = 0, length_red = 0, cnt_black = 0, length_black = 0, cnt_blue = 0, length_blue = 0;
    unsigned cnt_multi1 = 0, len_multi1 = 0, cnt_multi2 = 0, len_multi2 = 0, cnt_multi3 = 0, len_multi3 = 0, cnt_more_3 = 0, len_more_3 = 0;
    unsigned cnt_non_short = 0, cnt_red_non_short = 0, cnt_black_non_short = 0, cnt_blue_non_short = 0, cnt_multi1_non_short = 0, cnt_multi2_non_short = 0, cnt_multi3_non_short = 0, cnt_more_3_non_short = 0;
    unsigned cnt_imbalanced = 0, cnt_imbalanced_non_short = 0, len_imbalanced = 0;
    unsigned total_len = 0;
    file_dot << "digraph {\nnodesep = 0.5;\n";
    for (auto&& node : graph_modified) {
        file_dot << node.first << " [style=filled fillcolor=\"white\"]\n";
    }
    for (auto&& node : graph_modified) {
        std::string start_node = node.first;
        for (auto&& edges : node.second.outgoing_edges) {
            std::string end_node = edges.first;
            for (auto&& edge : edges.second) {
                if (this->genome_path_haplome3_strand1.genome_path.size())
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << '+' << edge.multi21 + edge.multi22 << '+' << edge.multi31 + edge.multi32 << ")";
                else if (this->genome_path_haplome2_strand1.genome_path.size())
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << '+' << edge.multi21 + edge.multi22 << ")";
                else
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << ")";
                assert(edge.multiplicity == edge.multi11 + edge.multi12 + edge.multi21 + edge.multi22 + edge.multi31 + edge.multi32);
                if (edge.multiplicity == 1) {
                    if (edge.multi11) {
                        if (this->genome_path_haplome1_strand1.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 1 and strand 1 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G1S1_" << this->genome_path_haplome1_strand1.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi12) {
                        if (this->genome_path_haplome1_strand2.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 1 and strand 2 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G1S2_" << this->genome_path_haplome1_strand2.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi21) {
                        if (this->genome_path_haplome2_strand1.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 2 and strand 1 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G2S1_" << this->genome_path_haplome2_strand1.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi22) {
                        if (this->genome_path_haplome2_strand2.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 2 and strand 2 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G2S2_" << this->genome_path_haplome2_strand2.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi31) {
                        if (this->genome_path_haplome3_strand1.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 3 and strand 1 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G3S1_" << this->genome_path_haplome3_strand1.node2order[start_node][end_node].at(0);
                    }
                    else if (edge.multi32) {
                        if (this->genome_path_haplome3_strand2.node2order[start_node][end_node].size() != 1)
                            std::cout << "Simple bulge in the same haplome 3 and strand 2 " << start_node << "->" << end_node << std::endl;
                        else
                            file_dot << "G3S2_" << this->genome_path_haplome3_strand2.node2order[start_node][end_node].at(0);
                    }
                }

                std::string color_out;
                if (edge.color == 0) {
                    cnt_black += 1;
                    length_black += edge.length;
                    color_out = "black";
                    // file_dot << "\", color=\"black\"";
                    if (edge.length > 300)
                        cnt_black_non_short += 1;
                }
                else if (edge.color == 1) {
                    cnt_red += 1;
                    length_red += edge.length;
                    color_out = "red";
                    // file_dot << "\", color=\"red\"";
                    if (edge.length > 300)
                        cnt_red_non_short += 1;
                }
                else if (edge.color == 2) {
                    cnt_blue += 1;
                    length_blue += edge.length;
                    color_out = "blue";
                    // file_dot << "\", color=\"blue\"";
                    if (edge.length > 300)
                        cnt_blue_non_short += 1;
                }
                else if (edge.color == 3) {
                    // cnt_blue += 1;
                    // length_blue += edge.length;
                    color_out = "purple";
                    // file_dot << "\", color=\"blue\"";
                    // if (edge.length > 300)
                    //     cnt_blue_non_short += 1;
                }
                else {
                    std::cout << "Error: edge cannot be assigned color (graph output)" << std::endl;
                }

                if (edge.block_id != 0)
                    file_dot << "\", color=\"azure3\"";
                else
                    file_dot << "\", color=\"" + color_out + "\"";

                total_len += edge.length * edge.multiplicity;

                if (edge.multi11 + edge.multi12 != edge.multi21 + edge.multi22) {
                    cnt_imbalanced += 1;
                    len_imbalanced += edge.length;
                    if (edge.length > 300)
                        cnt_imbalanced_non_short += 1;
                }

                if (edge.length >= thick) {
                    cnt_long += 1;
                    length_long += edge.length;
                    file_dot << ", penwidth=5]\n";
                }
                else {
                    cnt_short += 1;
                    length_short += edge.length;
                    file_dot << ", penwidth=1]\n";
                }

                if (edge.length > 300)
                    cnt_non_short += 1;

                if (edge.multiplicity == 1) {
                    cnt_multi1 += 1;
                    len_multi1 += edge.length;
                    if (edge.length > 300)
                        cnt_multi1_non_short += 1;
                }
                else if (edge.multiplicity == 2) {
                    cnt_multi2 += 1;
                    len_multi2 += edge.length;
                    if (edge.length > 300)
                        cnt_multi2_non_short += 1;
                }
                else if (edge.multiplicity == 3) {
                    cnt_multi3 += 1;
                    len_multi3 += edge.length;
                    if (edge.length > 300)
                        cnt_multi3_non_short += 1;
                }
                else {
                    cnt_more_3 += 1;
                    len_more_3 += edge.length;
                    if (edge.length > 300)
                        cnt_more_3_non_short += 1;
                }
                assert(edge.sequence.size() - this->k == edge.length);
                assert(edge.sequence.at(this->k) == edge.start_base);
                file_fasta << ">" << start_node << "_" << end_node << "_" << edge.start_base << edge.length << "(" << edge.multi11 + edge.multi12 << '+' << edge.multi21 + edge.multi22 << ")" << edge.color << "\n";
                file_fasta << edge.sequence << "\n";
            }
        }
    }
    file_dot << "}" << std::endl;
    file_dot.close();
    file_fasta.close();

    std::cout << "Total number of nodes: " << this->get_num_nodes() << std::endl;
    std::cout << "Total number of edges (non-short edges): " << cnt_long + cnt_short << " (" << cnt_non_short << ")" << std::endl;
    if (!this->genome_path_haplome2_strand1.genome_path.size()) {
        std::cout << "Length of red edges: " << length_red << std::endl;
    }
    else if (!this->genome_path_haplome1_strand1.genome_path.size()) {
        std::cout << "Length of blue edges: " << length_blue << std::endl;

    }
    else {
        std::cout << "Bicolored edges: " << cnt_black << std::endl;
        std::cout << "Length of bicolored edges: " << length_black << std::endl;
        std::cout << "Red edges: " << cnt_red << std::endl;
        std::cout << "Length of red edges: " << length_red << std::endl;
    }
    std::cout << "Multi 1 (non-short edges): " << cnt_multi1 << " (" << cnt_multi1_non_short << ")" << std::endl;
    std::cout << "Length: " << len_multi1 << std::endl;
    std::cout << "Multi 2 (non-short edges): " << cnt_multi2 << " (" << cnt_multi2_non_short << ")" << std::endl;
    std::cout << "Length: " << len_multi2 << std::endl;
    std::cout << "Multi 3 (non-short edges): " << cnt_multi3 << " (" << cnt_multi3_non_short << ")" << std::endl;
    std::cout << "Length: " << len_multi3 << std::endl;
    std::cout << "Multi >3 (non-short edges): " << cnt_more_3 << " (" << cnt_more_3_non_short << ")" << std::endl;
    std::cout << "Length: " << len_more_3 << std::endl;

    if (this->genome_path_haplome1_strand1.genome_path.size()) {
        std::cout << "Genome 1 edges: " << this->genome_path_haplome1_strand1.genome_path.size() - 1 << std::endl;
        std::cout << "Genome 1 length: " << len_genome1 + k << std::endl;
    }
    if (this->genome_path_haplome2_strand1.genome_path.size()) {
        std::cout << "Genome 2 edges: " << this->genome_path_haplome2_strand1.genome_path.size() - 1 << std::endl;
        std::cout << "Genome 2 length: " << len_genome2 + k << std::endl;
    }
    if (this->genome_path_haplome3_strand1.genome_path.size()) {
        std::cout << "Genome 3 edges: " << this->genome_path_haplome3_strand1.genome_path.size() - 1 << std::endl;
        std::cout << "Genome 3 length: " << len_genome3 + k << std::endl;
    }

    std::cout << "Imbalanced edges (non-short edges): " << cnt_imbalanced << " (" << cnt_imbalanced_non_short << ")" << std::endl;
    std::cout << "Length: " << len_imbalanced << std::endl;
}

// Get reverse complementary edge_sequence
std::string Graph::reverse_complementary(std::string& seq) {
    std::string out;
    for (int i = seq.size() - 1; i >= 0; --i) {
        if (seq.at(i) == 'A')
            out += 'T';
        else if (seq.at(i) == 'T')
            out += 'A';
        else if (seq.at(i) == 'C')
            out += 'G';
        else if (seq.at(i) == 'G')
            out += 'C';
        else
            out += seq.at(i);
    }
    return std::move(out);
}
