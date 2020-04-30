// define Graph::ReadFileAsMatrix(int n);
//
// April 18 - 12:13 PM -- Start on notes
// 			  12:45 PM -- Code
// 	          1:45  PM -- Break -- restart pc 
//            3:12  PM -- fucking i lost the file shits...
//            5:45  PM -- Break..

// April 19 - 1:30  PM -- Starting again

// April 21 - 10:27 AM -- Make it compile i hope.
//          - 11:46 AM -- Compiles! -- break
// After some depression... On April 29 we managed to read file and do some calcs!
// still a long way to go.

#include <cstdlib>
#include <vector>
#include "utils.cpp"
#define DEBUG_MODE true
#define DebugPrint(str ...) if(DEBUG_MODE) printf(str);
class UGraph 
{
	private:

	size_t n;
	size_t m;
    bool directed;
    std::vector<std::vector<size_t>> adj_lists;
	
    // This moves the edges from a removee and places them to a replacee
	// @todo check if works for Ugraph and Dgraph
	//
	// Assert that if e is an edge in removee --> e is an edge in replacee
	void RenameEdges(size_t replacee, size_t removee) {
		DebugPrint("Renaming Edges\n\tREMOVE= %d - REPLACEE= %d\n", removee, replacee)
		if (removee == replacee) {return;}

		// loop over ALL vertices -- count = n
		for (int x=0; x < n_vertices(); x++) {
			// loop over adj list of a vertex -- count ~ max |n|, min |0|.. the SUM of all this loop is constant, it is 2M
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
	size_t n_edges() {if (directed) return m; return m/2;}
	/* Constructors */
	UGraph() {
		this->n = 0;
		this->m = 0;
        this->directed = false;
	}
	
	// Copy Constructor
	UGraph(UGraph *pGraph) {
		// @todo copy constructor;
		this->directed = pGraph->directed;
		this->n = pGraph->n;
		this->m = pGraph->m;
		this->adj_lists = pGraph->adj_lists; // is this copying or is it moving pointers?
	}

	// Adjacency List Constructor
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

	/* Desctructors */
	

	/* Operations */
	void AddEdge(size_t i, size_t j) {
		DebugPrint("Adding Edge From %d to %d\n", i, j)
        adj_lists[i].push_back(j);
        if (!this->directed) {adj_lists[j].push_back(i); this->m++;}
		this->m++;
	}

	void RemoveSelfLoops(size_t v) {
		DebugPrint("Call to RemoveSelfLoops\n")
        auto s = adj_lists[v].size();
		if (adj_lists[v].size() == 0) { return; }
		DebugPrint("\ton vertex #%d\n", v)
		for(int i = 0; i < s; i++) {
			DebugPrint("\t Looking at #%d ...", adj_lists[v][i])
			if (adj_lists[v][i] == v) {
                adj_lists[v].erase(adj_lists[v].begin() + i);

				auto _t = this->m;
				this->m--;
				DebugPrint("\t found!\n\t-- m before = %d, m after = %d, fuck = %d :D\n", _t, this->m, _t -1)
			}
			DebugPrint(" Not found :(\n")
		}
	}
	
	void RemoveVertex(size_t i) {
		DebugPrint("\t Removing Vertex %d", i)
		adj_lists.erase(adj_lists.begin() + i);
		this->n--;
	}
	
	/* Operations on a UGraph */

	static void Contract(UGraph* g, size_t i, size_t j) {
		if (g->directed) {return; }
		if (g->n_vertices() <= i || g->n_vertices() <= j) { return;}

		g->RenameEdges(i,j);
		g->RemoveSelfLoops(i);
		g->RemoveVertex(j);
	}

	// This function copies internally g, (@todo then deletes it).
	// and is not an in_place function, @see KargerMinCut_inplace
	static void KargerMinCut(UGraph* g, size_t* p_cuts) {
		UGraph _graph = new UGraph(g);
		size_t v0, v1;
		size_t j=0, c=0, x;
		while (_graph.n_vertices() > 2) {
			size_t r = rand() % g->n_edges();
			while (r > x) {
				x = _graph.adj_lists[j].size();
				c = r - x; 
				j++;
			}
			v0 = _graph.adj_lists[j][c];
			v1 = j;
			Contract(&_graph, v0, v1);
		}
		size_t t = _graph.n_edges();
		memcpy(p_cuts, &t , sizeof(size_t));
	}
};
