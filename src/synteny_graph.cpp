#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "synteny_graph.hpp"
#include <cassert>
#include <algorithm>
#include <iomanip>

Edge::Edge() {}

Edge::Edge(char& base, const unsigned& length, const std::string& sequence, std::vector<int>& multis) {
    this->start_base = base;
    this->length = length;
    this->sequence = sequence;
    this->multis = multis;
    int sum = 0;
    for (auto&& m : multis)
        sum += m;
    this->multiplicity = sum;
    this->original_edge = false;
}

Node::Node() {}

Genome::Genome() {}

void Genome::add_node_to_end(std::string node) {
    this->genome_path.push_back(node);
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

void Graph::find_haplome(std::vector<std::string>& haplome_seqs, std::vector<std::string>& haplome_seqs_r,
    std::string edge_sequence, std::string& edge, std::string& edge_r, std::vector<int>& haplome_indices,
    std::vector<std::vector<int>>& positions_haplomes
) {
    haplome_indices.clear();
    positions_haplomes.clear();
    for (int i = 0; i < haplome_seqs.size(); ++i) {
        std::vector<int> positions_strand1, positions_strand2;
        this->get_locations(haplome_seqs[i], edge_sequence, positions_strand1);
        this->get_locations(haplome_seqs_r[i], edge_sequence, positions_strand2);
        positions_haplomes.emplace_back(positions_strand1);
        positions_haplomes.emplace_back(positions_strand2);
        haplome_indices.push_back(positions_strand1.size());
        haplome_indices.push_back(positions_strand2.size());
    }
}

void Graph::sort_locations(std::unordered_map<std::string, std::vector<std::vector<int>>>& edge2positions, std::vector<Genome>& genomes) {
    std::vector<std::vector<unsigned>> positions;
    std::vector<std::unordered_map<unsigned, std::string>> pos2edge;
    for (auto&& edge : edge2positions) {
        for (int i = 0; i < edge.second.size(); ++i) {
            std::vector<int>& positions_haplome = edge.second[i];
            std::sort(positions_haplome.begin(), positions_haplome.end());
            if (positions.size() < i + 1)
                positions.resize(i + 1);
            if (pos2edge.size() < i + 1)
                pos2edge.resize(i + 1);

            for (auto&& pos : positions_haplome) {
                positions[i].push_back(pos);
                pos2edge[i][pos] = edge.first;
            }
        }
    }

    // sort positions
    genome_paths.resize(positions.size());
    for (int i = 0; i < positions.size(); ++i) {
        std::unordered_map<unsigned, int> pos2id;
        std::sort(positions[i].begin(), positions[i].end());
        for (int j = 0; j < positions[i].size(); ++j) {
            pos2id[positions[i][j]] = j;
            std::string edge_index = pos2edge[i][positions[i][j]];
            size_t pos1 = edge_index.find(" -> ");
            std::string node1 = edge_index.substr(0, pos1);
            size_t pos2 = edge_index.find(" L ", pos1);
            std::string node2 = edge_index.substr(pos1 + 4, pos2 - pos1 - 4);
            if (!genome_paths[i].genome_path.empty()) {
                assert(genome_paths[i].genome_path.at(genome_paths[i].genome_path.size() - 1) == node1);
                genome_paths[i].add_node_to_end(node2);
            }
            else {
                genome_paths[i].add_node_to_end(node1);
                genome_paths[i].add_node_to_end(node2);
            }
        }
    }
}

// Read graph from DOT file
void Graph::read_from_dot(const std::string& graph_dot, const std::string& graph_fasta, const std::vector<std::string>& haplomes, int k) {
    // load fasta sequence
    std::string line, node1, node2, length, multiplicity, start_base;
    std::unordered_map<std::string, std::string> edge2sequence;
    std::unordered_map<std::string, std::vector<int>> edge2haplome;
    std::unordered_map<std::string, std::vector<std::vector<int>>> edge2positions;

    this->k = k;

    // read input haplomes
    std::vector<std::string> haplome_seqs, haplome_seqs_r;
    for (auto&& haplome : haplomes) {
        std::string haplome_seq;
        std::ifstream haplome_f(haplome);
        while (getline(haplome_f, line)) {
            if (line.at(0) != '>') haplome_seq += line;
        };
        haplome_f.close();
        for (auto&& c : haplome_seq) c = toupper(c);
        haplome_seqs.push_back(haplome_seq);
        haplome_seqs_r.push_back(reverse_complementary(haplome_seq));
    }
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
    std::vector<std::vector<int>> edge_indices_openmp(1000);
    std::vector<std::vector<int>> edge_r_indices_openmp(1000);
    std::vector<std::vector<std::vector<int>>> edge_positions_openmp(1000);
    std::vector<std::vector<std::vector<int>>> edge_r_positions_openmp(1000);

    unsigned cnt_line = 0;
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
                    this->find_haplome(haplome_seqs, haplome_seqs_r, edge_sequences_openmp[i], edges_openmp[i], edges_r_openmp[i],
                        edge_indices_openmp[i], edge_positions_openmp[i]);
                    this->find_haplome(haplome_seqs, haplome_seqs_r, reverse_complementary(edge_sequences_openmp[i]), edges_openmp[i], edges_r_openmp[i],
                        edge_r_indices_openmp[i], edge_r_positions_openmp[i]);
                }

                for (int i = 0; i < omp_size; ++i) {
                    edge2haplome[edges_openmp[i]] = edge_indices_openmp[i];
                    edge2haplome[edges_r_openmp[i]] = edge_r_indices_openmp[i];
                    assert(edge2haplome[edges_openmp[i]] == edge2haplome[edges_r_openmp[i]]);
                    edge2positions[edges_openmp[i]] = edge_positions_openmp[i];
                    edge2positions[edges_r_openmp[i]] = edge_r_positions_openmp[i];
                }
                omp_index = 0;
                edges_openmp.clear();
                edges_openmp.resize(omp_size);
                edges_r_openmp.clear();
                edges_r_openmp.resize(omp_size);
                edge_sequences_openmp.clear();
                edge_sequences_openmp.resize(omp_size);
                edge_indices_openmp.clear();
                edge_indices_openmp.resize(omp_size);
                edge_r_indices_openmp.clear();
                edge_r_indices_openmp.resize(omp_size);
            }

            if (++cnt_line % 1000 == 0)
                std::cout << "Processed " << cnt_line << " edge sequences." << std::endl;
        }
    }
    // dealing with remaining edges in the buffer
#pragma omp parallel for
    for (int i = 0; i < omp_index; i++) {
        this->find_haplome(haplome_seqs, haplome_seqs_r, edge_sequences_openmp[i], edges_openmp[i], edges_r_openmp[i],
            edge_indices_openmp[i], edge_positions_openmp[i]);
        this->find_haplome(haplome_seqs, haplome_seqs_r, reverse_complementary(edge_sequences_openmp[i]), edges_openmp[i], edges_r_openmp[i],
            edge_r_indices_openmp[i], edge_r_positions_openmp[i]);
    }

    for (int i = 0; i < omp_index; ++i) {
        edge2haplome[edges_openmp[i]] = edge_indices_openmp[i];
        edge2haplome[edges_r_openmp[i]] = edge_r_indices_openmp[i];
        edge2positions[edges_openmp[i]] = edge_positions_openmp[i];
        edge2positions[edges_r_openmp[i]] = edge_r_positions_openmp[i];
    }

    this->sort_locations(edge2positions, this->genome_paths);
    std::cout << "Read sequences from " << graph_fasta << " finished (loaded " << edge2sequence.size() << " edges, including reverse complements)." << std::endl;

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
        int multi_from_alignment = 0;
        for (auto&& pos : edge2positions[edge_index])
            multi_from_alignment += pos.size();
        if (multiplicity != multi_from_alignment) {
            std::cout << "Correct multiplicity of " << edge_index << " from " << multiplicity << " to " << multi_from_alignment << std::endl;
            multiplicity = multi_from_alignment;
        }
        Edge edge = Edge(start_base, length, edge2sequence[edge_index], edge2haplome[edge_index]);
        edge.original_edge = true;
        edge.haplome_positions = edge2positions[edge_index];
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
                    flag = true;
                    for (int j = 0; j < edge.multis.size() / 2; ++j) {
                        if (edge.multis[2 * j] + edge.multis[2 * j + 1] > x)
                            flag = false;
                    }
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

    std::vector<std::vector<Block_in_Genome_Path>> blocks_haplomes;
    for (int g = 0; g < this->genome_paths.size() / 2; ++g) {
        int g_i = g * 2;
        Genome& genome = this->genome_paths[g_i];
        int len_genome = 0;
        int previous_multi = -1;
        std::vector<Block_in_Genome_Path> blocks;
        for (int i = 0; i + 1 < genome.genome_path.size(); ++i) {
            std::string node1 = genome.genome_path[i], node2 = genome.genome_path[i + 1];

            // find edge index
            int edge_index = -1;
            int cnt = 0;
            for (int e = 0; e < graph_synteny_gp[node1].outgoing_edges[node2].size(); ++e) {
                if (graph_synteny_gp[node1].outgoing_edges[node2].at(e).multis[g_i] > 0) {
                    edge_index = e;
                    break;
                }
            }
            int min_loc = -1;
            for (int e = 0; e < graph_synteny_gp[node1].outgoing_edges[node2].size(); ++e) {
                if (graph_synteny_gp[node1].outgoing_edges[node2].at(e).original_edge && graph_synteny_gp[node1].outgoing_edges[node2].at(e).haplome_positions[g_i].size() > 0) {
                    cnt += 1;
                    if (min_loc == -1 || graph_synteny_gp[node1].outgoing_edges[node2].at(e).haplome_positions[g_i].at(0) < min_loc) {
                        min_loc = graph_synteny_gp[node1].outgoing_edges[node2].at(e).haplome_positions[g_i].at(0);
                        edge_index = e;
                    }
                }
            }
            assert(edge_index != -1);

            Edge current_edge = graph_synteny[node1].outgoing_edges[node2].at(edge_index);
            if (current_edge.multiplicity >= previous_multi) {
                Block_in_Genome_Path block(current_edge.block_id, this->synteny_blocks[current_edge.block_id].multiplicity, this->synteny_blocks[current_edge.block_id].length, node1, node2, len_genome, len_genome + current_edge.length, true, false, this->synteny_blocks[current_edge.block_id].num_edges);
                blocks.emplace_back(block);
            }
            else {
                Block_in_Genome_Path block(current_edge.block_id, this->synteny_blocks[current_edge.block_id].multiplicity, this->synteny_blocks[current_edge.block_id].length, node1, node2, len_genome + k, len_genome + current_edge.length, false, false, this->synteny_blocks[current_edge.block_id].num_edges);
                blocks[blocks.size() - 1].include_right = true;
                blocks[blocks.size() - 1].end += k;
                blocks.emplace_back(block);
            }
            previous_multi = current_edge.multiplicity;
            len_genome += current_edge.length;
            graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multiplicity -= 1;
            graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).multis[g_i] -= 1;
            if (min_loc != -1) {
                graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).haplome_positions[g_i].erase(graph_synteny_gp[node1].outgoing_edges[node2].at(edge_index).haplome_positions[g_i].begin());
            }
        }
        if (blocks.size()) {
            blocks[blocks.size() - 1].include_right = true;
            blocks[blocks.size() - 1].end += k;
            merge_blocks(blocks);
        }
        blocks_haplomes.emplace_back(blocks);
    }

    std::unordered_map<int, int> id_map;
    int id_new = 1;
    for (auto&& blocks : blocks_haplomes) {
        for (auto&& b : blocks) {
            if (id_map.find(b.id) != id_map.end())
                b.id = id_map.at(b.id);
            else {
                id_map[b.id] = id_new;
                id_map[-b.id] = -id_new;
                b.id = id_new;
                id_new += 1;
            }
        }
    }

    for (int i = 0; i < blocks_haplomes.size(); ++i) {
        if (!blocks_haplomes[i].empty()) {
            std::string haplome_blocks = prefix + ".haplome" + std::to_string(i + 1) + ".blocks";
            std::ofstream file_haplome_blocks(haplome_blocks);
            for (int j = 0; j < blocks_haplomes[i].size(); ++j) {
                file_haplome_blocks << blocks_haplomes[i][j].id << '(' << blocks_haplomes[i][j].multiplicity << ')' << std::setprecision(1) << std::fixed << 1.0 * (blocks_haplomes[i][j].end - blocks_haplomes[i][j].start) / 1000 << '\t' << blocks_haplomes[i][j].node_s << '\t' << blocks_haplomes[i][j].node_e << '\t' << blocks_haplomes[i][j].start << '\t' << blocks_haplomes[i][j].end << '\t' << blocks_haplomes[i][j].include_left << '\t' << blocks_haplomes[i][j].include_right << '\t' << blocks_haplomes[i][j].is_short_version << std::endl;
            }
            file_haplome_blocks.close();
        }
    }
}

void Graph::write_graph(std::unordered_map<std::string, Node>& graph_modified, const std::string& prefix, int thick) {
    std::string graph_modified_dot = prefix + ".dot";
    std::string graph_modified_fasta = prefix + ".fasta";

    std::cout << "\nWrite graph " << graph_modified_dot << ", fasta " << graph_modified_fasta << std::endl;
    std::ofstream file_dot(graph_modified_dot);
    std::ofstream file_fasta(graph_modified_fasta);

    for (auto&& genome_path : this->genome_paths)
        genome_path.construct_edge_orders();
    std::cout << genome_paths[0].node2order.size() << std::endl;
    std::cout << genome_paths[1].node2order.size() << std::endl;
    std::cout << genome_paths[2].node2order.size() << std::endl;
    std::cout << genome_paths[3].node2order.size() << std::endl;
    std::cout << genome_paths[4].node2order.size() << std::endl;
    std::cout << genome_paths[5].node2order.size() << std::endl;


    file_dot << "digraph {\nnodesep = 0.5;\n";
    for (auto&& node : graph_modified) {
        file_dot << node.first << " [style=filled fillcolor=\"white\"]\n";
    }
    for (auto&& node : graph_modified) {
        std::string start_node = node.first;
        for (auto&& edges : node.second.outgoing_edges) {
            std::string end_node = edges.first;
            for (auto&& edge : edges.second) {
                file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.start_base << edge.length << "(" << edge.multiplicity << ")";
                if (edge.multiplicity == 1) {
                    for (int i = 0; i < edge.multis.size(); ++i) {
                        if (edge.multis[i]) {
                            if (this->genome_paths[i].node2order[start_node][end_node].size() != 1)
                                std::cout << "Simple bulge in the same haplome " << start_node << "->" << end_node << "(" << this->genome_paths[i].node2order[start_node][end_node].size() << " " << edge.multis[i] << ")" << std::endl;
                            else
                                file_dot << "G" << (i + 2) / 2 << "S" << i % 2 + 1 << "_" << this->genome_paths[i].node2order[start_node][end_node].at(0);
                        }
                    }
                }

                if (edge.length >= thick) {
                    file_dot << "\", penwidth=5]\n";
                }
                else {
                    file_dot << "\", penwidth=1]\n";
                }
                assert(edge.sequence.size() - this->k == edge.length);
                assert(edge.sequence.at(this->k) == edge.start_base);
                file_fasta << ">" << start_node << "_" << end_node << "_" << edge.start_base << edge.length << "(" << edge.multiplicity << ")" << "\n";
                file_fasta << edge.sequence << "\n";
            }
        }
    }
    file_dot << "}" << std::endl;
    file_dot.close();
    file_fasta.close();
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

std::string Graph::reverse_complementary_node(std::string node) {
    if (node.at(0) == '-')
        return node.substr(1);
    else
        return '-' + node;
}