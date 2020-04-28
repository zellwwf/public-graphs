// define linked lists ?? isn't adjacency lists just link lists? or just arrays?
// define Graph::ReadFileAsList(int n);
// define Graph::
//
// define Graph::ReadFileAsMatrix(int n);
//
// define graph Graph.Contract(int i, int j);
//		return -1 if i || j > Graph.n_nodes
//		sort i, j -- assume from now i < j;
//		loop on Graph._AL[j] as k 
//			Graph._AL[i] << k unless k == i;
//		delete Graph._AL[j]
//		return :suceess
//
// define int Graph::KargerMinCut(Graph g);
// 		g->n_vertices <= 2 return the edges between them
//      else:
//		select a random edge: Edge E.
//		g = Contract(E.v, E.u);
//		return KargerMinCut(g);
//
// April 18 - 12:13 PM -- Start on notes
// 			  12:45 PM -- Code
// 	          1:45  PM -- Break -- restart pc 
//            3:12  PM -- fucking i lost the file shits...
//            5:45  PM -- Break..

// April 19 - 1:30  PM -- Starting again

// April 21 - 10:27 AM -- Make it compile i hope.
//          - 11:46 AM -- Compiles! -- break
#include <cstdlib>
#include <vector>
#include "utils.cpp"

class UGraph 
{
	
	private:

	size_t n;
	size_t m;
    bool directed;
    std::vector<std::vector<size_t>> adj_lists;
	
    // This moves the edges from a removee and places them to a replacee
	void RenameEdges(size_t removee, size_t replacee) {
		if (removee == replacee) {return;}

		for (int x=0; x < n_vertices(); x++) {
			for(int y = 0; y < adj_lists[x].size(); y++) {
				if (adj_lists[x][y] == removee) {
					adj_lists[x][y] = replacee;
				}
			}
		}
		return;
	}

	public:
	size_t n_vertices() {return n;}
	size_t n_edges() {return m;}
	/* Constructors */
	UGraph() {
		this->n = 0;
		this->m = 0;
        this->directed = false;
	}
	
	UGraph(UGraph *pGraph) {
		// @todo copy constructor;
		this->n = pGraph->n;
		this->m = pGraph->m;
		this->adj_lists = pGraph->adj_lists;
		// memcopy all the vectors in adj list
	}

	UGraph(bool is_dag, const std::vector<std::vector<size_t>> adj_lists) {
		size_t _m = 0;

		for (auto i = 0; i < adj_lists.size(); i ++) {
			_m += adj_lists[i].size();
		}

		this->n = adj_lists.size();
		this->m = _m;
		this->directed = is_dag;
		this->adj_lists = adj_lists;
	}
	/* Operations */
	void AddEdge(size_t i, size_t j) {
        adj_lists[i].push_back(j);
        if (!this->directed) {adj_lists[j].push_back(i); this->m++;}
		this->m++;
	}

	void RemoveSelfLoops(size_t v) {
        auto s = adj_lists[v].size();
		if (adj_lists[v].size() == 0) { return; }
		for(int i = 0; i < s; i++) {
			if (adj_lists[v][i] == v) {
                adj_lists[v].erase(adj_lists[v].begin() + i);
				this->m--;
			}
		}
	}
	
	void RemoveVertex(size_t i) {
		adj_lists.erase(adj_lists.begin() + i);
		this->n--;
	}
	
	void Contract(UGraph* g, size_t i, size_t j) {
		if (g->n_vertices() <= i || g->n_vertices() <= j) { return;}

		auto list_size = adj_lists[j].size();
		for (auto k = 0; k < list_size; k++) {
			g->AddEdge(i, adj_lists[j][k]);
		}
		
		RenameEdges(i,j);
		RemoveSelfLoops(i);
		g->RemoveVertex(j);
	}

	// @for fucks sake please note that this assumes:
	// 1. The UGRAPH is dia
	void KargerMinCut(UGraph* g, size_t* n_min_cut) {
		UGraph _graph = new UGraph(g);
		size_t v0, v1;
		size_t j=0, c=0, x;
		while (g->n_vertices() > 2) {
			size_t r = rand() % g->n_edges();
			while (r > x) {
				x = adj_lists[j].size();
				c = r - x; 
				j++;
			}
			v0 = adj_lists[j][c];
			v1 = j;
			_graph.Contract(&_graph, v0, v1);
		}
	}
};
