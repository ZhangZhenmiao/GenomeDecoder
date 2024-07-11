#include <iostream>
#include <sys/stat.h>
#include "synteny_graph.hpp"
#include "cmdline/cmdline.h"

std::vector<std::string> splitStringByComma(const std::string& str) {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string item;

    while (std::getline(ss, item, ',')) {
        result.push_back(item);
    }

    return result;
}

int main(int argc, char* argv[]) {
    cmdline::parser argParser;
    argParser.add<std::string>("dot", 'd', "graph.dot file under LJA output", true);
    argParser.add<std::string>("fasta", 'f', "the graph.fasta file under LJA output", true);
    argParser.add<std::string>("haplomes", 'h', "the reference haplome sequences", true);
    argParser.add<int>("kmer", 'k', "the k-mer size of LJA", true);
    argParser.add<std::string>("output", 'o', "the output directory (should be new)", true);

    argParser.parse_check(argc, argv);
    std::string graph_dot = argParser.get<std::string>("dot");
    std::string graph_fasta = argParser.get<std::string>("fasta");
    std::vector<std::string> haplomes = splitStringByComma(argParser.get<std::string>("haplomes"));
    int k = argParser.get<int>("kmer");
    std::string output = argParser.get<std::string>("output");

    Graph graph;

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
    graph.read_from_dot(graph_dot, graph_fasta, haplomes, k);
    graph.write_graph(graph.graph, output + "/graph.with_haplome");
    graph.synteny_block_generation(output + "/graph.with_haplome");

    return 0;
}
