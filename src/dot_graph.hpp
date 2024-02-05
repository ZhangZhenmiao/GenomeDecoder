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
    int color = -1;
    char start_base;
    unsigned length;
    std::string sequence;
    int multiplicity = 0;
    int max_multi_in_non_branching = 0;
    int multi11 = 0, multi12 = 0, multi21 = 0, multi22 = 0, multi31 = 0, multi32 = 0;

    // locations only for original edge
    bool original_edge = false;
    std::vector<int> multiplicity11;
    std::vector<int> multiplicity12;
    std::vector<int> multiplicity21;
    std::vector<int> multiplicity22;
    std::vector<int> multiplicity31;
    std::vector<int> multiplicity32;

    // synteny blocks
    int block_id = 0;

    void add_multi_from_edge_or_path(Edge& edge_or_path);
    void add_multi_from_edge_or_path(Path& edge_or_path);
    void update_color_by_multi();
    void remove_multi_from_path(Path& edge_or_path);
    bool check_equal_multi(Edge& edge_or_path);
    bool check_equal_multi(Path& edge_or_path);
    bool check_folds_multi(Edge& edge_or_path);
    bool check_contain_multi(Path& edge_or_path);
    Edge(char& base, const unsigned& length, const std::string& sequence, int multi11, int multi12, int multi21, int multi22, int multi31, int multi32);
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
    bool find_path(Path& path, std::vector<int>& pos);
    Genome();
};

class Path {
public:
    std::vector<std::string> nodes;
    std::vector<int> bulge_legs;
    std::unordered_map<int, std::string> index2palindromic_seq;

    std::string sequence;
    unsigned length = 0;
    int color = -1;
    int multiplicity = -1;
    int multi11 = 0, multi12 = 0, multi21 = 0, multi22 = 0, multi31 = 0, multi32 = 0;

    std::vector<int> index_haplome1_strand1;
    std::vector<int> index_haplome1_strand2;
    std::vector<int> index_haplome2_strand1;
    std::vector<int> index_haplome2_strand2;
    std::vector<int> index_haplome3_strand1;
    std::vector<int> index_haplome3_strand2;

    bool update_min_multi(Edge edge_or_path);
    void update_color_by_multi();
    Path();
};

class Bulge {
public:
    Path leg1, leg2;
    double max_identity = 0;
    int num_edges = 0;
    Bulge();
    Bulge(Path p1, Path p2, double identity);
};

class Graph {
public:
    int k;
    double sim_simple;
    std::unordered_map<std::string, Node> graph;
    Genome genome_path_haplome1_strand1;
    Genome genome_path_haplome1_strand2;
    Genome genome_path_haplome2_strand1;
    Genome genome_path_haplome2_strand2;
    Genome genome_path_haplome3_strand1;
    Genome genome_path_haplome3_strand2;
    std::unordered_map<int, Block> synteny_blocks;

    void read_from_dot(const std::string& graph_dot, const std::string& graph_fasta, const std::string& haplome1, const std::string& haplome2, const std::string& haplome3, int k);
    void write_graph(const std::string& prefix, int thick = 50000, int min_synteny_length = 2000);
    void write_graph(std::unordered_map<std::string, Node>& graph, const std::string& prefix, int thick = 50000, int min_synteny_length = 2000);
    void write_graph_abnormal(std::unordered_map<std::string, Node>& graph, const std::string& prefix);
    int get_num_nodes();

    void multi_bulge_removal(unsigned& removed_bulges, bool check_overlap = true, double similarity = 0);
    void merge_non_branching_paths();
    void erase_node_from_genome_paths(std::string node);

    void general_whirl_removal(unsigned& removed_whirls, bool simple_whirl = false);
    void resolving_bulge_with_two_multi_edge_paths(unsigned& removed_paths, int x, double identity, bool verbose = false, std::string output = "");

    void merge_blocks(std::vector<Block_in_Genome_Path>& blocks);
    void synteny_block_generation(const std::string& prefix);
    void count_duplicaction_graph(const std::string& prefix);
    void count_imbalanced_graph(const std::string& prefix);

    Graph();
private:
    void get_locations(std::string& ref, std::string& seq, std::vector<int>& locations);
    void find_haplome(std::string& haplome1_seq, std::string& haplome1_seq_r,
        std::string& haplome2_seq, std::string& haplome2_seq_r, std::string& haplome3_seq, std::string& haplome3_seq_r,
        std::string edge_sequence, std::string& edge, std::string& edge_r, int& haplome_index,
        std::vector<int>& positions_haplome1_strand1, std::vector<int>& positions_haplome1_strand2,
        std::vector<int>& positions_haplome2_strand1, std::vector<int>& positions_haplome2_strand2,
        std::vector<int>& positions_haplome3_strand1, std::vector<int>& positions_haplome3_strand2
    );
    void clear_buffer(auto& buffer, const int omp_size);
    void sort_locations(std::unordered_map<std::string, std::vector<int>>& edge2positions, Genome& genome);

    bool check_reverse_bulge(std::string node1, std::string node2, bool check_overlap = true);
    std::string collapse_bulge(std::string node1, std::string node2, unsigned& removed_bulges, double similarity = 0);

    bool check_non_branching(std::string node);
    void merge_edges(std::string node);
    bool add_node_to_path(Path& path, std::string node, int bulge_leg = 0);

    void remove_whirl(Path& unambiguous_path, std::vector<std::string>& nodes_to_remove);
    bool find_path_from_genome_paths(Path& path);
    bool replace_genome_path_with_path(Path& path, Path& path_replace, bool only_check = false);

    std::string collapse_complex_bulge_two_multi_edge_paths(Path& p1, Path& p2, double max_identity, std::vector<std::string>& nodes_to_remove, bool verbose = false, std::string output = "");

    // synteny block generation
    void mask_edges_with_multi_x(std::unordered_map<std::string, Node>& graph, int x, int& start_id);
    bool mask_non_branching_path(std::unordered_map<std::string, Node>& graph, int x, int& start_id);
    void collapse_non_branching_path(std::unordered_map<std::string, Node>& graph, std::string node, std::unordered_set<std::string>& nodes_to_remove);
    void change_multi_non_branching_path(std::unordered_map<std::string, Node>& graph, std::string node, std::vector<std::string>& non_branching_path, int& max_multi);


    bool get_reverse_path(Path& path, Path& path_reverse);

    double PI_edlib(std::string sequence1, std::string sequence2);
    double percent_identity(std::string cigar);
    std::string reverse_complementary_node(std::string node);
    std::string reverse_complementary(std::string& s);

    void dfs(std::unordered_map<std::string, Node>& graph, std::string node, std::unordered_set<std::string>& component);
};