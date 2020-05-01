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
#include <algorithm> 
#include "utils.cpp"
#define DEBUG_MODE true
#define DebugPrint(str ...) if(DEBUG_MODE) printf(str);
#define WarnPrint(str ...) printf(str);

class UGraph 
{
	private:

	size_t n;
	size_t m;
    bool directed;
    std::vector<std::vector<size_t>> adj_lists;
	std::vector<bool> is_v_deleted;

    // This moves the edges from a removee and places them to a replacee
	// @todo check if works for Ugraph and Dgraph
	//
	// Assert that if e is an edge in removee --> e is an edge in replacee
	void RenameEdges(size_t replacee, size_t removee) {
		if (removee == replacee) {return;}
		if (is_v_deleted[removee] || is_v_deleted[replacee]) {return;}

		DebugPrint("Renaming Edges\n\tsource: %d - target: %d\n", removee+1, replacee+1)
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
	size_t n_edges() {if (directed) return m; return (m/2 - m%2);}
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
		this->is_v_deleted = pGraph->is_v_deleted;
	}

	// =Operator

	// new operator

	// Adjacency List Constructor
	UGraph(bool is_dag, const std::vector<std::vector<size_t>> adj_lists) {
		size_t _m = 0;

		for (auto i = 0; i < adj_lists.size(); i ++) {
			_m += adj_lists[i].size();
		}

		this->n = adj_lists.size();
		this->is_v_deleted = std::vector<bool>(adj_lists.size(), false);
		this->m = _m;
		this->directed = is_dag;
		this->adj_lists = adj_lists;
	}

	/* Desctructors */
	
	/* Displayers */

	//in honor of ruby
	void to_s() {
		//   0: 1 2 3 
		size_t nvertices = adj_lists.size();
		size_t i,j;

		printf("\n== G has %d vertices and %d edges ==\n------------------------------\n", n_vertices(), n_edges());
		for (i = 0;i < nvertices; i++) {
			size_t m_edges = adj_lists[i].size(); 
			printf("\t%d ", i+1);
			if (is_v_deleted[i]) {
				printf(" -- deleted -- ");
			}
			else 
			{
				for (j = 0; j < m_edges; j++) {
					printf(" %5d", adj_lists[i][j] + 1);
				}
			}
			printf("\n");
		}
	}

	/* Operations */

	// Adds a single edge between two vertices.
	void AddEdge(size_t i, size_t j) {
		DebugPrint("Adding Edge (%d, %d)\n", i+1, j+1)
		if (is_v_deleted[i] || is_v_deleted[j]) {return;} // #add warns messages
        adj_lists[i].push_back(j);
        if (!this->directed) {adj_lists[j].push_back(i); this->m++;}
		this->m++;
	}

	void CopyEdges(size_t src, size_t target) {
		if (is_v_deleted[src] || is_v_deleted[target]) {return;} // #add warns messages
		DebugPrint("\t Copying Edges of %d to %d\n", src+1, target+1)
		adj_lists[target].insert(adj_lists[target].end(), adj_lists[src].begin(), adj_lists[src].end());
		this->m += adj_lists[src].size();
	}

	void RemoveSelfLoops(size_t v) {
		if (is_v_deleted[v]) return;
		DebugPrint("Call to RemoveSelfLoops\n")
        auto s = adj_lists[v].size();
		if (adj_lists[v].size() == 0) return;
		DebugPrint("\ton vertex #%d\n", v+1)
		for(int i = 0; i < s; i++) {
			DebugPrint("\t Looking at #%d ...", adj_lists[v][i]+1)
			if (adj_lists[v][i] == v) {
                adj_lists[v].erase(adj_lists[v].begin() + i);
				this->m--;
			}
			DebugPrint(" Not found :(\n")
		}
	}
	
	void RemoveVertex(size_t i) {
		if (is_v_deleted[i]) return;
		DebugPrint("\t Removing Vertex %d\n", i+1)
		size_t s = adj_lists[i].size();

		// remove references to the removed vertex from other lists.
		// in effect this deletes edges
		for (int x = 0; x < s; x++) {
			auto ref = adj_lists[i][x];
			adj_lists[ref].erase(std::remove(adj_lists[ref].begin(), adj_lists[ref].end(), i), adj_lists[ref].end());
		}


		// set of edges where s is source.
		// set of edges where s is target.
		if (this->directed) {
			// @todo
		} else {
			this->m -= s;
		}
		adj_lists[i] = std::vector<size_t>(0);
		is_v_deleted[i] = true;
		this->n--;
		this->m -= s;
	}
	
	/* Operations on a UGraph */

	static void Contract(UGraph* g, size_t i, size_t j) {
		if (g->directed)return;
		if (i == j) return;
		if (g->is_v_deleted[i] || g->is_v_deleted[j]) return;
		if (g->n_vertices() <= i || g->n_vertices() <= j) { return;}

		DebugPrint("\t Contracting on Edge (%d, %d) (current edges#: %d)\n", i+1, j+1, g->n_edges())

		g->CopyEdges(j,i);
		g->RenameEdges(i,j); 
		g->RemoveSelfLoops(i); 
		g->RemoveVertex(j);
	}

	// This function copies internally g, (@todo then deletes it).
	// and is not an in_place function, @see KargerMinCut_inplace
	// it writes back the size of minimum cut in the the p_cuts var. 
	static void KargerMinCut(UGraph* g, size_t* p_cuts) {
		DebugPrint("Call: KargerMinCut\n")

		UGraph _graph = new UGraph(g); // @note this is used as a keyword here, not the operator new.
		size_t v0, v1; // holds the values of the randomly selected edge (v0, v1)
		size_t j=0, c, x;
		DebugPrint("  m  |  n  |  r  |  j  |  x  |  c  |  v0  |  v1  |\n")
		DebugPrint("--------------------------------------------------\n")
		while (_graph.n_vertices() > 2) {
			size_t r = rand() % g->n_edges();
			c = r;
			x = _graph.adj_lists[j].size();
			while (c >= x) {				
				c -= x; 
				j++;
				x = _graph.adj_lists[j].size();
			}
			v0 = _graph.adj_lists[j][c];
			v1 = j;
			DebugPrint("%5d|%5d|%5d|%5d|%5d|%5d|%5d|%5d\n", g->n_edges(), g->n_vertices(), r, j, x, c, v0+1, v1+1)

			Contract(&_graph, v0, v1);
		}
		size_t t = _graph.n_edges();
		memcpy(p_cuts, &t , sizeof(size_t));
	}
};
