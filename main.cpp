#include <iostream>
#include <vector>
#include <string.h>

#include "digraph.cpp"

int main(int argc, char *argv[]) {
    std::string filename; 

    std::cout << "ShitTY Graph ProgrAm! \n";
    std::cout << "-> Filename: ";
    std::cin >> filename;

    auto x = Digraph::ParseFileToAdjacency(filename);
    x.to_s();
    x.ReverseEdges();
    x.to_s();
    x.ReverseEdges();
    x.to_s();
    Digraph::StartDFS(&x, 3);

    for (auto i = 0; i < x.visitation_order.size(); i++) {
        std::cout << i + 1 << " --> " << x.visitation_order[i] << "\n";
    }

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