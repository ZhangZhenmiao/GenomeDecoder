#include "edlib.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

std::string align_by_edlib(std::string sequence1, std::string sequence2) {
    EdlibAlignResult result = edlibAlign(sequence1.c_str(), sequence1.size(), sequence2.c_str(), sequence2.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
        std::string cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
        for (auto&& letter : cigar) {
            if (letter == '=')
                letter = 'M';
        }
        return cigar;
    }
    else {
        std::cout << "edlib failed" << std::endl;
        return "";
    }
}

int main(int argc, char* argv[]) {
    std::string file1 = argv[1];
    std::string file2 = argv[2];
    std::string out = argv[3];
    std::string line;
    std::string sequence1, sequence2;

    std::ifstream seq1_f(file1);
    while (getline(seq1_f, line)) {
        if (line.at(0) == '>')
            continue;
        else
            sequence1 += line;
    }
    seq1_f.close();
    std::cout << sequence1.size() << std::endl;

    std::ifstream seq2_f(file2);
    while (getline(seq2_f, line)) {
        if (line.at(0) == '>')
            continue;
        else
            sequence2 += line;
    }
    seq2_f.close();
    std::cout << sequence2.size() << std::endl;

    std::string cigar = align_by_edlib(sequence2, sequence1);

    if (mkdir(out.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) {
        if (errno == EEXIST) {}
        else {
            std::cout << "Create output directory faliled." << std::endl;
            exit(1);
        }
    }

    std::ofstream out_f(out + "/cigar.txt");
    out_f << cigar;
    out_f.close();

    return 0;
}