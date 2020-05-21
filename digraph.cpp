#include <cstdlib>
#include <vector>
#include <iostream>
#include <algorithm> 
#include <string.h>
#include <map>
#include <numeric> // USED ONLY FOR IOTA! in Dijsktra

#define DEBUG
#define VERBOSITY 2 // 0, 1, 2 -- (errors, <-- + warns, <-- + info )

#ifdef DEBUG
#define DebugPrint(str ...) printf(str)
#else
#define DebugPrint(str ...) do { } while(0)
#endif
#define InfoPrint(str ...) if(VERBOSITY > 1) printf(str)
#define WarnPrint(str ...) if(VERBOSITY > 0) printf(str)
#define ErrorPrint(str ...) printf(str);

class Digraph
{
	private:

	size_t n;	// |V|
	size_t m;	// |E|

    std::vector<std::vector<size_t>> adj_lists;
	std::vector<std::vector<size_t>> weights;

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

	// Assignment Constructor Intentionally Left Blank (Why? wats the default behavior?)

	// Adjacency List Constructor
	Digraph(const std::vector<std::vector<size_t>> adj_lists) {
		DebugPrint("Call: Digraph Adj. List Constructor\n");
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

	Digraph(const std::vector<std::vector<size_t>> adj_lists, const std::vector<std::vector<size_t>> weights) {
		DebugPrint("Call: Digraph Adj. List Constructor\n");
		size_t _m = 0;
		size_t _n = adj_lists.size();

		for (auto i = 0; i < adj_lists.size(); i ++) {
			_m += adj_lists[i].size();
		}

		this->n = _n;
		this->is_v_deleted = std::vector<bool>(_n, false);
		this->m = _m;
		this->adj_lists = adj_lists;
		this->weights = weights;
	}
	
	/* Desctructors */
	
	/* Displayers */

	// in honor of ruby
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

	void SetEdgeWeights(std::vector<std::vector<size_t>> new_weights) {
		if (new_weights.size() == get_m()) {
			this->weights = new_weights;
		} else {
			WarnPrint("You passed an array of size: %d and we have %d edges only!", new_weights.size(), get_m());
		}
	}

    // This moves the edges from a removee and places them to a replacee
	// @todo check if works for Ugraph and Dgraph
	//
	// Assert that if e is an edge in removee --> e is an edge in replacee
	void RenameEdges(const size_t replacee, size_t removee) {
		DebugPrint("Call: RenameEdges\n\tsource: %d - target: %d\n", removee+1, replacee+1);
		if (removee == replacee) {
			WarnPrint("ERROR: Same Vertex\n");
			return;
			}
		if (is_v_deleted[removee] || is_v_deleted[replacee]) {
			WarnPrint("ERROR: Operating On a Deleted Vertex\n");
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
		DebugPrint("Adding Edge (%d, %d)\n", i+1, j+1);
		if (is_v_deleted[i] || is_v_deleted[j]) {
			WarnPrint("ERROR: Operating On a Deleted Vertex\n");
			return;
		}

        adj_lists[i].push_back(j);
		this->m++;
        adj_lists[j].push_back(i); 
		this->m++;
	}

	/*
	@deprecated
	*/
	void CopyEdges(size_t src, size_t target) {
		DebugPrint("Call: Copying Edges of %d to %d\n", src+1, target+1);
		if (is_v_deleted[src] || is_v_deleted[target]) {
			WarnPrint("ERROR: Operating On a Deleted Vertex\n");
			return;
		}

		// this block is responsible for making sure the set of edges out of a source, will be copied to target
		// currently its broken for digraphs since it reverse the direction
		auto n_refs = adj_lists[src].size();
		for (auto i = 0; i < n_refs; i++) {
			auto ref = adj_lists[src][i];
			this->adj_lists[target].push_back(ref);
			this->adj_lists[ref].push_back(target); // this makes this method not workj for digraphs @todo @deprecated
			this->m++;
		}

		adj_lists[target].insert(adj_lists[target].end(), adj_lists[src].begin(), adj_lists[src].end());
		this->m += adj_lists[src].size();
	}

	void RemoveSelfLoops(size_t v) {
		DebugPrint("Call: RemoveSelfLoops\n");
		if (is_v_deleted[v]) {
			WarnPrint("ERROR: Operating On a Deleted Vertex\n");
			return;
		}
        auto old_size = adj_lists[v].size();

		if (adj_lists[v].size() == 0) {
			WarnPrint("ERROR: Operating On a Disconnected\n");
			return;
		}
		adj_lists[v].erase(std::remove(adj_lists[v].begin(), adj_lists[v].end(), v), adj_lists[v].end());
		auto diff = old_size - adj_lists[v].size();
		this-> m -= diff;
	}
	
	// @deprecated -- needs rework on digraphs
	void RemoveVertex(size_t i) {
		DebugPrint("Call: Removing Vertex %d\n", i+1);
		if (is_v_deleted[i]) {
			return;
		}
		size_t n_out_degree = adj_lists[i].size();

		// remove references to the removed vertex from other lists.
		// in effect this deletes edges of the vertex i before deleting it.
		// Erase (Remove) line t == n_edges_of_i
		// For loops with the erase_remove line is going to be linear in n of edges, t(m)
		// the loop iterates exactly in_degree of (i) times.
		for (int x = 0; x < n_out_degree; x++) {
			auto ref = adj_lists[i][x];
			adj_lists[ref].erase(std::remove(adj_lists[ref].begin(), adj_lists[ref].end(), i), adj_lists[ref].end());
		}

		// set of edges where s is source.
		// set of edges where s is target.
		this->m -= n_out_degree;
		adj_lists[i] = std::vector<size_t>(0);	// this is to perseve the ordering (indexing) of adj_lists, for now.
		is_v_deleted[i] = true;
		this->n--;
		this->m -= n_out_degree;
	}
	
	// @deprecated -- needs rework on digraphs
	static void Contract(Digraph* g, size_t i, size_t j) {
		DebugPrint("Call: Contract (remainder: %d, removed: %d)\n", i+1, j+1);
		if (i == j) {
			WarnPrint("ERROR: Same Vertex\n");
			return;
		}
		if (g->is_v_deleted[i] || g->is_v_deleted[j]) {
			WarnPrint("ERROR: Contracting on a Deleted Vertex\n");
			return;
		}

		DebugPrint("\t Contracting on Edge (%d, %d) (current edges: %d)\n", i+1, j+1, g->n_edges());

		g->CopyEdges(j,i);
		g->RemoveSelfLoops(i); 
		g->RemoveVertex(j);
	}

	static void KargerMinCut(Digraph* g, size_t* p_cuts) {
		DebugPrint("Call: KargerMinCut\n");

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
	std::vector<size_t> scc_sets;
	size_t scc_current_leader;
	std::vector<size_t> shorted_path_weights;

	static void Kosaraju_SCC(Digraph* g) {
		DebugPrint("Call: Kosaraju...\n");
		g->ReverseEdges();
		InfoPrint("Going to First DFS Loop");
		DFS_Loop(g, true);
		g->ReverseEdges();
		InfoPrint("Going to 2nd DFS Loop");
		DFS_Loop(g, false);
	}

	/*
	
	*/
	static void DFS_Loop(Digraph* g, bool initial_pass) {
		DebugPrint("Call: DFS Loop...\n as pass: %d\n", initial_pass ? 1 : 2);
		auto n = g->n_vertices();
		std::vector<size_t> reorder = std::vector<size_t>(n);
		g->dfs_running = true;
		g->unexplored = std::vector<bool>(n,true);

		// if your on the first pass, scc set will store the 'magical' ordering u will need on the reverse graph.
		g->visit_iterator = n;	// please note that visits and current leader are starting by 1
		g->scc_current_leader = n;
		if (initial_pass) {
			g->scc_sets = std::vector<size_t>(n,0);
			g->visitation_order = std::vector<size_t>(n,0);
		}

		if (!initial_pass) {
			// some O(n) preprocessing to make shit faster
			for (auto k = 0; k < n; k++) {
				reorder[g->visitation_order[k]-1] = k; 
			}
		}

		for (auto i=0; i < n; i++) {
			auto vertex = initial_pass ? i : reorder[i];
			if (g->unexplored[vertex]) {

				g->scc_current_leader = i+1;
				g->DFS(vertex, initial_pass);

				if (initial_pass) {
					g->visitation_order[vertex] = g->visit_iterator;
					g->visit_iterator--;
				}
			}
		}

		if (!initial_pass) {
			g->dfs_running = false;
		}
	}

	static void StartDFS(Digraph* g, size_t vertex) {
		DebugPrint("Call: DFS Start...\n");
		auto n = g->n_vertices();
		g->dfs_running = true;
		g->unexplored = std::vector<bool>(n,true);
		g->visitation_order = std::vector<size_t>(n,0);
		g->visit_iterator = n;
		g->DFS(vertex, true);
		g->visitation_order[vertex] = g->visit_iterator;
		g->visit_iterator--;
	}

	/*
	DFS
	---

	T: 			
	S_max: 		ω(2*(n+m)) 
	
	n : colorings + explorings
				coloring + explorings + visitation ordering
	*/
	void DFS(size_t vertex, bool initial_pass) {
		DebugPrint("DFS on %d\n", vertex);
		scc_sets[vertex] = scc_current_leader;	// coloring strongly connected components with the same color.
		// reserve mem for stack
		unexplored[vertex] = false;

		size_t n_edges = adj_lists[vertex].size(); // the number of edges of in v.
		for (auto i = 0; i < n_edges; i++) {
			auto new_v = adj_lists[vertex][i];
			// solve this: find x such that adj_lists[new_v][x] = vertex;
			if (unexplored[new_v]) { 
				DFS(new_v, initial_pass);
				if (initial_pass) {
					visitation_order[new_v] = visit_iterator;
					visit_iterator--;
				}
			}
		}
	}

	/*
	Reverses The Edges of the Graph
	-------------------------------
	Applies this rule on the adj list array: if adj[x][y] = z --> adj[y][x] = z
	
	visit each vertex (1st loop) and find its outgoing edges. Then iterate over the found arrows
	(2nd for loop) and reverse their direction.

	T: 			θ(m)
	S:			θ(n+m)
	S_max: 		ω(2*(n+m)) 

	we can also do the reversing in place, but due to the adjacency list construction, we will increase time complextiy
	in favor or memory. 
	*/
	void ReverseEdges() {
		auto n_vertices = this->adj_lists.size();
		auto new_edges = std::vector<std::vector<size_t>> (n_vertices);
		size_t n_edges = 0;

		for (auto i = 0; i < n_vertices; i++) {

			n_edges = adj_lists[i].size();

			for (auto j = 0; j < n_edges; j++) {
				auto target = adj_lists[i][j];
				new_edges[target].push_back(i); 
			}
		}

		this->adj_lists = new_edges;
	}


	/*
	@return std::vector<size_t> of all the 'destination' / 'adjacent' wathave u vertices.
	*/
	std::vector<size_t> GetAdjList_Subgraph(const std::vector<size_t> vertices) {
		size_t adj_size;
		std::vector<size_t> c_adj_list;
		std::vector<size_t> v;	
		std::map<size_t, bool> map;

		for (auto i=0; i<vertices.size() ; i++) {
			c_adj_list = adj_lists[vertices[i]];
			adj_size = c_adj_list.size();
			for (auto j=0; j< adj_size; j++) {
				map.insert(std::pair<size_t,bool>(c_adj_list[j], true));
			}
		}

		printf("hi");
		for(std::map<size_t,bool>::iterator it = map.begin(); it != map.end(); ++it) {
  			v.push_back(it->first);
		}

		return v;
	}
	/*

	*/
    std::vector<size_t> DijkstrasShortestPath(size_t v0) {
		DebugPrint("DijkstrasShortestPath on %d\n", v0);

		// Define: vector of shorted_paths_vals, initialized to infinity (or flip all bits to 1s)
		size_t infinity = std::numeric_limits<size_t>::max();
		std::vector<size_t> shortest_paths_vals = std::vector<size_t>(n, infinity);

		// Define: the Abyss (Where sp is not defined, infinite)
		std::vector<size_t> abyss = std::vector<size_t>(n);
		size_t start = 0;
		std::iota(abyss.begin(), abyss.end(), start); // [1,2,3,4,5...n-1]
		abyss.erase(abyss.begin() + v0); // Thanks for counting from 0 to n.

		// Define: Dijsktra's Subtree
		std::vector<size_t> subtree = std::vector<size_t>();
		subtree.reserve(n);

		shortest_paths_vals[v0] = 0;
		subtree.push_back(v0);

		size_t vp;	// the vertex we came from, the parent.
		size_t v = v0;
		size_t c_n_edges, c=0;
		std::vector<size_t> c_outgoing;
		std::vector<size_t> c_weights;
		//std::priority_queue <size_t, std::vector<size_t>, std::greater<size_t>> c_weights;
	
		size_t lowest_weight, v_dest, v_greedy, new_s_path, idx_lowest_weight;
		bool c_i_in_abyss;

		// WARNING this loop assumes that v0 is connected to each (other) vertex. 
		while (subtree.size() != n) {
			c_outgoing = GetAdjList_Subgraph(subtree);
			c_n_edges = c_outgoing.size();
			c_weights.clear();
			c_weights.reserve(c_n_edges);

			for (auto i = 0; i < c_n_edges; i++) {
				// pick the one with the lowest weight.
				v_dest = c_outgoing[i];

				// Go through the outgoing edges. Note that getting the 'crossing' set of edges
				// (def as: { (x,y) | x in subtree, y in abyss }) is not θ(1).
				// @review: I think if work on finding the crossing set of edges i can make it 
				// faster, without priority queue even? I think i can make it θ(log(n)) if we sort adj list.

				// STD::find assumes UNSORTED (thus fucking slower linear search)
				c_i_in_abyss = std::find(abyss.begin(), abyss.end(),v_dest) != abyss.end();

				// if (v, v_dest) is a crossing edge, place it in the heap.
				if (c_i_in_abyss) { 
					new_s_path = shortest_paths_vals[v] + weights[v][v_dest];
					if (shortest_paths_vals[v_dest] > new_s_path) {shortest_paths_vals[v_dest] = new_s_path;}
					
					c_weights[i] = shortest_paths_vals[v_dest];
					//c_weights.push(shorted_paths_vals[v_dest]); 
				} // else just skip the rest of the loop and check the next outgoing edge.
			}
			lowest_weight = infinity;
			idx_lowest_weight;
			for (int j = 0; j < c_n_edges; j++) {
				if (lowest_weight > c_weights[j]) {
					lowest_weight = c_weights[j];
					idx_lowest_weight = j;
				}
			}
			v_greedy = idx_lowest_weight;

			// Apply the greedy logic and follow the lowest path:
			vp = v;
			v = v_greedy;
			subtree.push_back(v);
		}

		//shortest_paths_vals 
		return shortest_paths_vals;
	}
	/*

	*/
	static Digraph ParseFileToAdjacency(const std::string filename) {
		bool skip_first_entry = true;
		size_t no_lines = 1;
		size_t line_len = 0;
		size_t i = 0;
		size_t _t;
		FILE* pFile = NULL;
		char* line = NULL;
		char* input_char_buffer = NULL;
		std::vector<std::vector<size_t>> adj_lists;
		std::vector<size_t> t_list;

		pFile = fopen(filename.c_str(), "r");
		if (pFile == NULL) { return Digraph(); }

		printf("* reading file: %s\n");

		while (getline(&line, &line_len, pFile) != -1)
		{
			skip_first_entry = true;
			input_char_buffer = strtok(line, " \t");

			if (input_char_buffer == NULL) { 
				std::cout << " input_char_buffer is empty! file reading issue...\n";
				fclose(pFile);
			}

			while (input_char_buffer != NULL) {
				if (skip_first_entry) { skip_first_entry = false;} 
				else {
					_t = (size_t) atoi(input_char_buffer);
					if (_t > 0) {
						t_list.push_back(_t - 1);
					}
				}
				input_char_buffer = strtok(NULL, "  \t");
			}

			adj_lists.push_back(t_list);
			t_list.clear();
			no_lines ++;
		}

		std::cout << "Read " << no_lines << "\n";
		Digraph g = Digraph(adj_lists);
		return g;
	}

	static Digraph ParseFileToAdjacencyWeighted(const std::string filename) {
		printf("Call: ParseFileWeighted");
		bool skip_first_entry = true;
		size_t no_lines = 1;
		size_t line_len = 0;
		size_t i = 0;
		size_t _t;
		FILE* pFile = NULL;
		char* line = NULL;
		char* input_char_buffer = NULL;
		char *low = NULL;
		char *high = NULL;
		size_t bsize;
		size_t csize = sizeof(char);

		unsigned short idx_comma;
		std::vector<std::vector<size_t>> adj_lists;
		std::vector<std::vector<size_t>> adj_lists_weights;
		std::vector<size_t> t_list;
		std::vector<size_t> t_weights;

		pFile = fopen(filename.c_str(), "r");
		if (pFile == NULL) { return Digraph(); }

		printf("* reading file: %s\n");

		while (getline(&line, &line_len, pFile) != -1)
		{
			skip_first_entry = true;
			input_char_buffer = strtok(line, " \t"); // "1 2,3\t"
			if (input_char_buffer == NULL) { 
				std::cout << " input_char_buffer is empty! file reading issue...\n";
				fclose(pFile);
			}

			while (input_char_buffer != NULL) {
				if (skip_first_entry) { 
					skip_first_entry = false;
				} 
				else 
				{
					bsize = sizeof(input_char_buffer);
					for (auto _k = 0; _k<bsize; _k++) {
						if (input_char_buffer[_k]==',') {
							idx_comma = _k;
							break;
						}
					} // ur doing valgrind and thinking about rounding in max/2
					
					low = (char*)malloc(csize * (1 + (bsize/2)));
					high = (char*)malloc(csize * (1 + (bsize/2)));
					
					strncpy(low, input_char_buffer, idx_comma-1);
					strcpy(high, &input_char_buffer[idx_comma+1]);

					_t = (size_t) atoi(low);
					if (_t > 0) {
						t_list.push_back(_t - 1);
						_t = (size_t) atoi(high);
						t_weights.push_back(_t);
					}
				}
				input_char_buffer = strtok(NULL, "  \t");
			}

			adj_lists.push_back(t_list);
			adj_lists_weights.push_back(t_weights);
			t_list.clear();
			t_weights.clear();
			no_lines ++;
		}
		
		int ss = 0;

		for (auto i = 0; i < 10; i++) {
			printf("- %d: ", i);
			ss = adj_lists[i].size();
			for (auto j = 0; j < ss; j++) {
				
				printf("%d (%d),",adj_lists[i][j], adj_lists_weights[i][j]);
			}
			printf("\n");
		}
		
		printf("* Read (%d) lines\n", no_lines);
		Digraph g = Digraph(adj_lists, adj_lists_weights);
		return g;
	}
};
