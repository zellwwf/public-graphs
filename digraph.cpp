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
#include <iostream>
#include <algorithm> 
#include <string>

#define DEBUG_MODE 0
#define DebugPrint(str ...) if(DEBUG_MODE) printf(str);
#define WarnPrint(str ...) printf(str);

class Digraph 
{
	private:

	size_t n;
	size_t m;

    std::vector<std::vector<size_t>> adj_lists;
	std::vector<bool> is_v_deleted;

	public:

	size_t n_vertices() {return n;}
	size_t n_edges() {return m;}
	size_t get_m() {return m;}

	/* Constructors */
	Digraph() {
		this->n = 0;
		this->m = 0;
		// @todo add the pointers / arrays... they will be problematic.
	}
	
	// Copy Constructor
	Digraph(Digraph *pGraph) {
		this->n = pGraph->n;
		this->m = pGraph->m;
		// @todo probably a bug: is this copying or is it moving pointers?
		this->adj_lists = pGraph->adj_lists;
		this->is_v_deleted = pGraph->is_v_deleted;
	}

	// Adjacency List Constructor
	Digraph(const std::vector<std::vector<size_t>> adj_lists) {
		DebugPrint("Call: Digraph Adj. List Constructor\n")
		size_t _m = 0;
		size_t _n = adj_lists.size();

		for (auto i = 0; i < adj_lists.size(); i ++) {
			_m += adj_lists[i].size();
		}

		this->n = _n;
		this->is_v_deleted = std::vector<bool>(_n, false);
		this->m = _m;
		this->adj_lists = adj_lists; 
	}

	/* Desctructors */
	
	/* Displayers */

	//in honor of ruby
	void to_s() {
		//   0: 1 2 3 
		size_t nvertices = adj_lists.size();
		size_t i,j;

		printf("\n* G has %d vertices and %d edges \n------------------------------\n", n_vertices(), n_edges());
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

    // This moves the edges from a removee and places them to a replacee
	// @todo check if works for Ugraph and Dgraph
	//
	// Assert that if e is an edge in removee --> e is an edge in replacee
	void RenameEdges(const size_t replacee, size_t removee) {
		DebugPrint("Call: RenameEdges\n\tsource: %d - target: %d\n", removee+1, replacee+1)
		if (removee == replacee) {
			WarnPrint("ERROR: Same Vertex\n")
			return;
			}
		if (is_v_deleted[removee] || is_v_deleted[replacee]) {
			WarnPrint("ERROR: Operating On a Deleted Vertex\n")
			return;
		}

		for (int x=0; x < n_vertices(); x++) {
			for(int y = 0; y < adj_lists[x].size(); y++) {
				if (adj_lists[x][y] == removee) {
					adj_lists[x][y] = replacee;
				}
			}
		}
		return;
	}

	// Adds a single edge between two vertices.
	void AddEdge(size_t i, size_t j) {
		DebugPrint("Adding Edge (%d, %d)\n", i+1, j+1)
		if (is_v_deleted[i] || is_v_deleted[j]) {
			WarnPrint("ERROR: Operating On a Deleted Vertex\n")
			return;
		}

        adj_lists[i].push_back(j);
		this->m++;
        adj_lists[j].push_back(i); 
		this->m++;
	}

	void CopyEdges(size_t src, size_t target) {
		DebugPrint("Call: Copying Edges of %d to %d\n", src+1, target+1)
		if (is_v_deleted[src] || is_v_deleted[target]) {
			WarnPrint("ERROR: Operating On a Deleted Vertex\n")
			return;
		}

		auto n_refs = adj_lists[src].size();
		for (auto i = 0; i < n_refs; i++) {
			auto ref = adj_lists[src][i];
			this->adj_lists[ref].push_back(target);
			this->m++;
		}

		adj_lists[target].insert(adj_lists[target].end(), adj_lists[src].begin(), adj_lists[src].end());
		this->m += adj_lists[src].size();
	}

	void RemoveSelfLoops(size_t v) {
		DebugPrint("Call: RemoveSelfLoops\n")
		if (is_v_deleted[v]) {
			WarnPrint("ERROR: Operating On a Deleted Vertex\n")
			return;
		}
        auto old_size = adj_lists[v].size();

		if (adj_lists[v].size() == 0) {
			WarnPrint("ERROR: Operating On a Disconnected\n")
			return;
		}
		adj_lists[v].erase(std::remove(adj_lists[v].begin(), adj_lists[v].end(), v), adj_lists[v].end());
		auto diff = old_size - adj_lists[v].size();
		this-> m -= diff;
	}
	
	void RemoveVertex(size_t i) {
		DebugPrint("Call: Removing Vertex %d\n", i+1)
		if (is_v_deleted[i]) {
			WarnPrint("ERROR: Attempting To Deleted a Deleted Vertex\n")
			return;
		}
		size_t s = adj_lists[i].size();

		// remove references to the removed vertex from other lists.
		// in effect this deletes edges
		for (int x = 0; x < s; x++) {
			auto ref = adj_lists[i][x];
			adj_lists[ref].erase(std::remove(adj_lists[ref].begin(), adj_lists[ref].end(), i), adj_lists[ref].end());
		}

		// set of edges where s is source.
		// set of edges where s is target.
		this->m -= s;
		adj_lists[i] = std::vector<size_t>(0);
		is_v_deleted[i] = true;
		this->n--;
		this->m -= s;
	}
	
	static void Contract(Digraph* g, size_t i, size_t j) {
		DebugPrint("Call: Contract (remainder: %d, removed: %d)\n", i+1, j+1)
		if (i == j) {
			WarnPrint("ERROR: Same Vertex\n")
			return;
		}
		if (g->is_v_deleted[i] || g->is_v_deleted[j]) {
			WarnPrint("ERROR: Contracting on a Deleted Vertex\n")
			return;
		}

		DebugPrint("\t Contracting on Edge (%d, %d) (current edges: %d)\n", i+1, j+1, g->n_edges())

		g->CopyEdges(j,i);
		g->RemoveSelfLoops(i); 
		g->RemoveVertex(j);
	}

	static void KargerMinCut(Digraph* g, size_t* p_cuts) {
		DebugPrint("Call: KargerMinCut\n")

		Digraph _graph = Digraph(g);
		size_t v0, v1; // holds the values of the randomly selected edge (v0, v1)
		size_t i, j, x;

		while (_graph.n_vertices() > 2) {
			size_t r = rand() % _graph.n_edges();

			i = 0; // Adj Row Index
			j = r; // Adj Column Index
			x = _graph.adj_lists[i].size();

			while (j >= x) {	
				i++;	// point to the next row		
				j -= x; // decrement the column counter
				x = _graph.adj_lists[i].size(); // fetch new size of row
			}
			v0 = _graph.adj_lists[i][j];
			v1 = i;

			Contract(&_graph, v1, v0);
		}

		size_t t = _graph.n_edges();
		memcpy(p_cuts, &t , sizeof(size_t));
	}

	// Dfs Vars
	std::vector<bool> unexplored;	// Maps vertices into explored or not.
	std::vector<size_t> visitation_order;
	size_t visit_iterator;
	bool is_acyclic;
	bool dfs_running;

	static void DFS_Loop(Digraph* g) {

	}

	static void StartDFS(Digraph* g, size_t vertex) {
		WarnPrint("Starting DFS...\n");
		auto n = g->n_vertices();
		g->dfs_running = true;
		g->unexplored = std::vector<bool>(n,true);
		g->visitation_order = std::vector<size_t>(n,0);
		g->visit_iterator = n;
		g->DFS(vertex);
		g->visitation_order[vertex] = g->visit_iterator;
		g->visit_iterator--;
	}

	void DFS(size_t vertex) {
		WarnPrint("DFS on %d\n", vertex);
		// reserve mem for stack
		unexplored[vertex] = false;

		size_t n_edges = adj_lists[vertex].size(); // the number of edges of in v.
		for (auto i = 0; i < n_edges; i++) {
			auto new_v = adj_lists[vertex][i];
			// solve this find x such that adj_lists[new_v][x] = vertex;
			auto nerww_v = adj_lists[i][vertex];
			if (unexplored[new_v]) { 
				DFS(new_v);
				visitation_order[new_v] = visit_iterator;
				visit_iterator--;
			}
		}
	}

	void ReverseEdges() {
		// if adj[x][y] = z --> adj[y][x] = z
		auto newEdges = std::vector<std::vector<size_t>> (this->adj_lists.size());

		// visit each vertex
		for (auto i = 0; i < adj_lists.size(); i++) {
			
			auto n_verts = adj_lists[i].size(); // 2    [1: 2 3]

			// 
			for (auto j = 0; j < n_verts; j++) {
				// put the reverse arrow 
				auto target = adj_lists[i][j];
				newEdges[target].push_back(i); 
			}

		}

		this->adj_lists = newEdges;
	}

	// G[x] --> edges out of x
	// G[x][y] --> edges out of x of [y] --> vertex
	// G[y][x] --> 
	static Digraph ParseFileToAdjacency(const std::string filename) {
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
		if (pFile == NULL) { return Digraph(); }

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
		Digraph g = Digraph(adj_lists);
		return g;
	}
};
