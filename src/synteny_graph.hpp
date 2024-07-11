#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
class Path;
class Edge;

class Block {
public:
    int multiplicity;
    int length;
    int num_edges;

    Block(int multiplicity, int length) {
        this->multiplicity = multiplicity;
        this->length = length;
    }
    Block() {};
};

class Block_in_Genome_Path {
public:
    int id;
    int multiplicity;
    int length;
    int num_edges;
    int start, end;
    std::string node_s, node_e;
    bool include_left = false;
    bool include_right = false;
    bool is_short_version = false;
    Block_in_Genome_Path(int id, int multiplicity, int length, std::string node_s, std::string node_e, int start, int end, bool include_left, bool include_right, int num_edges) {
        this->id = id;
        this->multiplicity = multiplicity;
        this->length = length;
        this->node_s = node_s;
        this->node_e = node_e;
        this->start = start;
        this->end = end;
        this->include_left = include_left;
        this->include_right = include_right;
        this->num_edges = num_edges;
    };
    Block_in_Genome_Path() {};
};

class Edge {
public:
    // haplome 1: 1, haplome 2: 2, both haplomes: 0, unknown: -1
    char start_base;
    unsigned length;
    std::string sequence;
    int multiplicity = 0;
    int max_multi_in_non_branching = 0;
    std::vector<int> multis;

    // locations only for original edge
    bool original_edge = false;
    std::vector<std::vector<int>> haplome_positions;

    // synteny blocks
    int block_id = 0;

    Edge(char& base, const unsigned& length, const std::string& sequence, std::vector<int>& multis);
    Edge();
};

class Node {
public:
    std::unordered_map<std::string, std::vector<Edge>> incoming_edges;
    std::unordered_map<std::string, std::vector<Edge>> outgoing_edges;
    Node();
};

class Genome {
public:
    std::vector<std::string> genome_path;
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<int>>> node2order;
    void add_node_to_end(std::string node);
    void erase_node(std::string node);
    void construct_edge_orders();
    Genome();
};

class Graph {
public:
    int k;
    std::unordered_map<std::string, Node> graph;
    std::vector<Genome> genome_paths;
    std::unordered_map<int, Block> synteny_blocks;

    void read_from_dot(const std::string& graph_dot, const std::string& graph_fasta, const std::vector<std::string>& haplomes, int k);
    void write_graph(std::unordered_map<std::string, Node>& graph, const std::string& prefix, int thick = 50000);
    int get_num_nodes();

    void merge_blocks(std::vector<Block_in_Genome_Path>& blocks);
    void synteny_block_generation(const std::string& prefix);
    void count_duplicaction_graph(const std::string& prefix);
    void count_imbalanced_graph(const std::string& prefix);

    Graph();
private:
    void get_locations(std::string& ref, std::string& seq, std::vector<int>& locations);
    void find_haplome(std::vector<std::string>& haplome_seqs, std::vector<std::string>& haplome_seqs_r,
        std::string edge_sequence, std::string& edge, std::string& edge_r, std::vector<int>& haplome_indices,
        std::vector<std::vector<int>>& positions_haplome
    );
    void sort_locations(std::unordered_map<std::string, std::vector<std::vector<int>>>& edge2positions, std::vector<Genome>& genomes);

    // synteny block generation
    void mask_edges_with_multi_x(std::unordered_map<std::string, Node>& graph, int x, int& start_id);
    bool mask_non_branching_path(std::unordered_map<std::string, Node>& graph, int x, int& start_id);
    void collapse_non_branching_path(std::unordered_map<std::string, Node>& graph, std::string node, std::unordered_set<std::string>& nodes_to_remove);
    void change_multi_non_branching_path(std::unordered_map<std::string, Node>& graph, std::string node, std::vector<std::string>& non_branching_path, int& max_multi);

    std::string reverse_complementary_node(std::string node);
    std::string reverse_complementary(std::string& s);
};