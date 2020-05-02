#include <iostream>
#include <vector>
#include <string.h>

#include "ugraph.cpp"


UGraph ParseFileToAdjacency(const bool directed, const std::string filename) {
    FILE* pFile = NULL;
    char* line = NULL;
    char* c_buff = NULL;
    bool skip_me = true;
    size_t no_lines = 1;
    size_t line_len = 0;
    size_t i = 0;
    size_t _t;

    std::vector<std::vector<size_t>> adj_lists;
    std::vector<size_t> t_list;

    pFile = fopen(filename.c_str(), "r");
    if (pFile == NULL) { return UGraph(); }

    std::cout << " > reading file: " << filename << "...\n";

    while (getline(&line, &line_len, pFile) != -1)
    {
        skip_me = true;
        c_buff = strtok(line, " \t");

        if (c_buff == NULL) { 
            std::cout << " c_buff is empty! file reading issue...\n";
            // @todo close filestream
        }

        while (c_buff != NULL) {
            if (skip_me) 
            { 
                skip_me = false;
            } else {
                _t = (size_t) atoi(c_buff);
                if (_t > 0) {
                    t_list.push_back(_t - 1);
                }
            }
            c_buff = strtok(NULL, "  \t");
        }

        adj_lists.push_back(t_list);
        t_list.clear();
        no_lines ++;
    }

    std::cout << "Read " << no_lines << "\n";
    UGraph g = UGraph(adj_lists);
    return g;
}


int main() {
    std::string filename; 

    std::cout << "ShitTY Graph ProgrAm! \n";
    std::cout << "-> Filename: ";
    std::cin >> filename;

    auto x = ParseFileToAdjacency(false, filename);

    size_t min_cut = 0;
    std::vector<size_t> vec = std::vector<size_t>(1500);

    for (int i = 0 ; i < 1500 ; i++) {
        UGraph::KargerMinCut(&x, &min_cut);
        vec[i] = min_cut;
    }

    size_t min = vec[0];

    for (int i = 0 ; i < 1500 ; i++) {
        std::cout << vec[i] << " \n";
        if (vec[i] < min) {
            min = vec[i];
        }
    }
    std::cout << " The minimum of a 100 is: "<< min << "\n";
    return 0;
}