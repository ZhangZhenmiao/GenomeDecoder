#include <iostream>
#include <sys/stat.h>
#include "dot_graph.hpp"
#include "graph_simplification.hpp"
#include "cmdline/cmdline.h"

int main(int argc, char* argv[]) {
    cmdline::parser argParser;
    argParser.add<std::string>("dot", 'd', "graph.dot file under LJA output", true);
    argParser.add<std::string>("fasta", 'f', "the graph.fasta file under LJA output", true);
    argParser.add<std::string>("haplome1", '1', "the reference haplome sequence 1", false);
    argParser.add<std::string>("haplome2", '2', "the reference haplome sequence 2", false);
    argParser.add<std::string>("haplome3", '3', "the reference haplome sequence 3", false);
    argParser.add<double>("similarity_simple", 's', "the similarity for simple bulges", false, 0);
    argParser.add<double>("similarity_complex", 'c', "the similarity for complex bulges", false, 0.65);
    argParser.add<int>("kmer", 'k', "the k-mer size of LJA", true);
    argParser.add<std::string>("output", 'o', "the output directory (should be new)", true);
    argParser.add("only_simple_bulge", '\0', "only conduct simple bulge collapsing");
    argParser.add("only_synteny", '\0', "only generate synteny blocks");


    argParser.parse_check(argc, argv);
    std::string graph_dot = argParser.get<std::string>("dot");
    std::string graph_fasta = argParser.get<std::string>("fasta");
    std::string haplome1 = argParser.get<std::string>("haplome1");
    std::string haplome2 = argParser.get<std::string>("haplome2");
    std::string haplome3 = argParser.get<std::string>("haplome3");
    double sim = argParser.get<double>("similarity_complex");
    double sim_simple = argParser.get<double>("similarity_simple");
    int k = argParser.get<int>("kmer");
    std::string output = argParser.get<std::string>("output");
    bool only_simple_bulge = argParser.exist("only_simple_bulge");
    bool only_synteny = argParser.exist("only_synteny");

    Graph graph;
    graph.sim_simple = sim_simple;

    if (mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) {
        if (errno == EEXIST) {
            std::cout << "Output directory already exist." << std::endl;
            return 0;
        }
        else {
            std::cout << "Create output directory faliled." << std::endl;
            return 0;
        }
    }

    // read graph from file
    graph.read_from_dot(graph_dot, graph_fasta, haplome1, haplome2, haplome3, k);
    graph.write_graph(output + "/graph.with_haplome");
    graph.synteny_block_generation(output + "/graph.with_haplome");
    graph.count_duplicaction_graph(output + "/graph.with_haplome");
    graph.count_imbalanced_graph(output + "/graph.with_haplome");

    if (only_synteny)
        return 0;

    unsigned removed_bulges = 1, removed_whirls = 1, total_whirls = 0;
    std::cout << "Original graph size: " << graph.get_num_nodes() << std::endl;
    std::cout << std::endl;

    if (only_simple_bulge) {
        graph.multi_bulge_removal(removed_bulges, true, sim_simple);
        std::cout << "Graph size after simple bulge removal: " << graph.get_num_nodes() << ", removed_bulges: " << removed_bulges << std::endl;
    }
    else {
        std::cout << "Bulge and whirl removal:" << std::endl;
        removed_bulges = 1;
        int cnt_rounds = 0;
        while (removed_bulges != 0) {
            cnt_rounds++;
            graph.multi_bulge_removal(removed_bulges, true, sim_simple);
            std::cout << "Graph size after simple bulge removal: " << graph.get_num_nodes() << ", removed_bulges: " << removed_bulges << std::endl;
            // if (removed_bulges != 0)
            //     graph.write_graph(output + "/graph.simple_bulge_round" + std::to_string(cnt_rounds));
            removed_whirls = 1;
            while (removed_whirls != 0) {
                graph.general_whirl_removal(removed_whirls, true);
                total_whirls += removed_whirls;
                std::cout << "Graph size after whirl removal: " << graph.get_num_nodes() << " removed_whirls: " << removed_whirls << std::endl;
            }
            std::cout << "Round " << cnt_rounds << " finished" << std::endl;
        }
        // graph.write_graph(output + "/graph.before_complex_bulge");
        std::cout << "Removed whirls: " << total_whirls << std::endl;
        std::cout << "Before complex bulge finished" << std::endl;
        std::cout << std::endl;

        std::cout << "Complex bulges 2+2" << std::endl;
        unsigned removed_paths = 1;
        cnt_rounds = 0;
        int total_removed = 0;
        total_whirls = 0;
        while (removed_paths != 0) {
            cnt_rounds++;
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 2, sim, false, output);
            std::cout << "Graph size after remove bulges with two multi-edge paths: " << graph.get_num_nodes() << ", removed_paths: " << removed_paths << std::endl;

            // graph.write_graph("debug_before_whirl");
            removed_whirls = 1;
            while (removed_whirls != 0) {
                graph.general_whirl_removal(removed_whirls, true);
                total_whirls += removed_whirls;
                std::cout << "Graph size after whirl removal: " << graph.get_num_nodes() << " removed_whirls: " << removed_whirls << std::endl;
                // graph.write_graph("debug_after_whirl");
            }
            total_removed += removed_paths;
            std::cout << "Round " << cnt_rounds << " finished" << std::endl;
        }
        // graph.write_graph(output + "/graph.complex_bulge_x2");
        std::cout << "Removed whirls: " << total_whirls << std::endl;
        std::cout << "Removed complex bulges: " << total_removed << std::endl;
        std::cout << std::endl;

        std::cout << "Complex bulges 3+3" << std::endl;
        removed_paths = 1;
        cnt_rounds = 0;
        total_removed = 0;
        total_whirls = 0;
        while (removed_paths != 0) {
            cnt_rounds++;
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 3, sim, false, output);
            std::cout << "Graph size after remove bulges with two multi-edge paths: " << graph.get_num_nodes() << ", removed_paths: " << removed_paths << std::endl;

            removed_whirls = 1;
            while (removed_whirls != 0) {
                graph.general_whirl_removal(removed_whirls, true);
                total_whirls += removed_whirls;
                std::cout << "Graph size after whirl removal: " << graph.get_num_nodes() << " removed_whirls: " << removed_whirls << std::endl;
            }

            total_removed += removed_paths;
            std::cout << "Round " << cnt_rounds << " finished" << std::endl;
        }
        // graph.write_graph(output + "/graph.complex_bulge_x3");
        std::cout << "Removed whirls: " << total_whirls << std::endl;
        std::cout << "Removed complex bulges: " << total_removed << std::endl;
        std::cout << std::endl;

        std::cout << "Complex bulges 4+4" << std::endl;
        removed_paths = 1;
        cnt_rounds = 0;
        total_removed = 0;
        total_whirls = 0;
        while (removed_paths != 0) {
            cnt_rounds++;
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 4, sim, false, output);
            std::cout << "Graph size after remove bulges with two multi-edge paths: " << graph.get_num_nodes() << ", removed_paths: " << removed_paths << std::endl;

            removed_whirls = 1;
            while (removed_whirls != 0) {
                graph.general_whirl_removal(removed_whirls, true);
                total_whirls += removed_whirls;
                std::cout << "Graph size after whirl removal: " << graph.get_num_nodes() << " removed_whirls: " << removed_whirls << std::endl;
            }

            total_removed += removed_paths;
            std::cout << "Round " << cnt_rounds << " finished" << std::endl;
        }
        // graph.write_graph(output + "/graph.complex_bulge_x4");
        std::cout << "Removed whirls: " << total_whirls << std::endl;
        std::cout << "Removed complex bulges: " << total_removed << std::endl;
    }
    graph.write_graph(output + "/graph.final");
    return 0;
}
