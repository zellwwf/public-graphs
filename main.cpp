#include <iostream>
#include <vector>
#include <string.h>

#include "ugraph.cpp"


void ParseFileToAdjacency(const bool directed, const std::string filename) {
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
    if (pFile == NULL) { return; }

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

    UGraph g = UGraph(false, adj_lists);
    g.to_s();

    g.RemoveVertex(2);
    g.to_s();


    g.RemoveVertex(2);
    g.to_s();

    g.RemoveVertex(3);
    g.to_s();
    
    g.RemoveVertex(1);
    g.to_s();

    g.RemoveVertex(0);
    g.to_s();

    UGraph::Contract(&g, (size_t)1, (size_t)0);
    g.to_s();
    
    UGraph::Contract(&g, (size_t)1, (size_t)0);
    g.to_s();

    char _tft_ = getchar();

    size_t x;
    /*
    for (auto l = 0; l < 15; l++) {
        g.KargerMinCut(&g, &x);
        std::cout << "Minimum cut size of graph is: " << x << "\n";
    }
    */
}


int main() {
    std::string filename; 

    std::cout << "ShitTY Graph ProgrAm! \n";
    std::cout << "-> Filename: ";
    std::cin >> filename;

    ParseFileToAdjacency(false, filename);
    return 0;
}