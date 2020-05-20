#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>

#include "digraph.cpp"

//cli.h
void command_loop() {

}

//cli.h 
void parse_args() {
    // @todo TBIL
    command_loop();
}

//cli.h
void say_hi() {
    std::cout << "LesS ShiTty Graph Program! \n";
}

int main(int argc, char *argv[]) {

    say_hi();
    
    if (argc > 0) {
        // parse args
        parse_args();
    } else {
        // show commands
        command_loop();
    }

    return 0;
}

void command_kosaraju() {
    std::string filename; 

    std::cout << "-> Filename: ";
    std::cin >> filename;

    auto x = Digraph::ParseFileToAdjacency(filename);

    Digraph::Kosaraju_SCC(&x);

    std::ofstream MyFile("colorings.txt");
    std::cout << "Writing out to file \n";

    for (auto i = 0; i < x.scc_sets.size(); i++) {
         MyFile << x.scc_sets[i] << "\n";
    }

    MyFile.close();
    std::cout << "Dump Done! Bye \n";
}

void command_mincuts() {
    size_t min_cut = 0;
    size_t karger_max_iterations = 0;
    std::string filename;

    std::cout << "Adjacency List File: ";
    std::cin >> filename;

    std::cout << "Number of Iteration";
    std::cin >> karger_max_iterations;

    Digraph x = Digraph::ParseFileToAdjacency(filename);
    std::vector<size_t> cuts_sample = std::vector<size_t>(karger_max_iterations);

    for (int i = 0 ; i < karger_max_iterations; i++) {
        Digraph::KargerMinCut(&x, &min_cut);
        cuts_sample[i] = min_cut;
    }

    size_t min = cuts_sample[0];

    for (int i = 0 ; i < 100 ; i++) {
        std::cout << cuts_sample[i] << " \n";
        if (cuts_sample[i] < min)   min = cuts_sample[i];
    }
    std::cout << " The minimum of a 100 is: "<< min << "\n";
}