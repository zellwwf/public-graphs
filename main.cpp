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

    std::string filename; 

    std::cout << "-> Filename: ";
    std::getline (std::cin,filename);

    auto x = Digraph::ParseFileToAdjacencyWeighted(filename);

    x.to_s();
    std::pair<std::vector<size_t>, std::vector<size_t>> res = Digraph::Dijkstra(&x, 0);
    std::cout << "============== Dijkistra Done =============== \n";

    auto a1 = res.first;
    auto a2 = res.second;
    
    for (auto i = 0; i < a1.size(); i++) {
         std::cout << a1[i] << " - " << a2[i] << "\n";
    }

    size_t sample_out[10];
    size_t sample_idx[10] = {7,37,59,82,99,115,133,165,188,197};

    std::cout << std::endl;
    for (auto i = 0; i < 10; i++) {
        std::cout << sample_idx[i] << ',';
    }

    std::cout << std::endl;
    for (auto i = 0; i < 10; i++) {
        std::cout << a1[sample_idx[i] - 1] << ',';
    }

    std::cout << std::endl;
    for (auto i = 0; i < 10; i++) {
        std::cout << a2[sample_idx[i] - 1] << ',';
    }

    std::cout << "\n\nDump Done! Bye \n";
    /*    
    if (argc > 0) {
        // parse args
        parse_args();
    } else {
        // show commands
        command_loop();
    }
    */

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