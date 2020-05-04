#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>

#include "digraph.cpp"

void getcommand() {
}

int main(int argc, char *argv[]) {
    std::string filename; 

    std::cout << "Shitty Graph Program! \n";
    std::cout << "-> Filename: ";
    std::cin >> filename;

    auto x = Digraph::ParseFileToAdjacency(filename);

    Digraph::Kosaraju_SCC(&x);
    //for (auto i = 0; i < x.visitation_order.size(); i++) {
    //     std::cout << i << " - " << x.visitation_order[i] << " - "<< x.scc_sets[i] << "\n";
    //}

    
    std::ofstream MyFile("colorings.txt");
    std::cout << "Writing out to file \n";

    for (auto i = 0; i < x.scc_sets.size(); i++) {
         MyFile << x.scc_sets[i] << "\n";
    }

    MyFile.close();
    std::cout << "Dump Done! Bye \n";

    return 0;
}

void mincutsss() {

    std::string filename; 

    std::cout << "ShitTY Graph ProgrAm! \n";
    std::cout << "-> Filename: ";
    std::cin >> filename;

    auto x = Digraph::ParseFileToAdjacency(filename);
    size_t min_cut = 0;
    std::vector<size_t> vec = std::vector<size_t>(100);

    for (int i = 0 ; i < 100 ; i++) {
        Digraph::KargerMinCut(&x, &min_cut);
        vec[i] = min_cut;
    }

    size_t min = vec[0];

    for (int i = 0 ; i < 100 ; i++) {
        std::cout << vec[i] << " \n";
        if (vec[i] < min) {
            min = vec[i];
        }
    }
    std::cout << " The minimum of a 100 is: "<< min << "\n";
}