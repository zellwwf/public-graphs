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

		for (size_t i = 0; i < adj_lists.size(); i ++) {
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

		for (size_t i = 0; i < adj_lists.size(); i ++) {
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
		printf("\n* Graph G has %ld vertices and %ld edges \n  -------------------------------------------\n", n_vertices(), n_edges());
		for (i = 0;i < nvertices; i++) {
			size_t m_edges = adj_lists[i].size(); 
			printf("\t%ld ", i+1);
			if (is_v_deleted[i]) {
				printf(" -- deleted -- ");
			}
			else 
			{
				for (j = 0; j < m_edges; j++) {
					printf(" %5ld", adj_lists[i][j] + 1);
					if (weights.size() > 0) {
						printf(" (%5ld)", weights[i][j]);
					}
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
			WarnPrint("You passed an array of size: %ld and we have %ld edges only!", new_weights.size(), get_m());
		}
	}

    // This moves the edges from a removee and places them to a replacee
	// @todo check if works for Ugraph and Dgraph
	//
	// Assert that if e is an edge in removee --> e is an edge in replacee
	void RenameEdges(const size_t replacee, size_t removee) {
		DebugPrint("Call: RenameEdges\n\tsource: %ld - target: %ld\n", removee+1, replacee+1);
		if (removee == replacee) {
			WarnPrint("ERROR: Same Vertex\n");
			return;
			}
		if (is_v_deleted[removee] || is_v_deleted[replacee]) {
			WarnPrint("ERROR: Operating On a Deleted Vertex\n");
			return;
		}

		for (size_t x=0; x < n_vertices(); x++) {
			for(size_t y = 0; y < adj_lists[x].size(); y++) {
				if (adj_lists[x][y] == removee) {
					adj_lists[x][y] = replacee;
				}
			}
		}
		return;
	}

	// Adds a single edge between two vertices.
	void AddEdge(size_t i, size_t j) {
		DebugPrint("Adding Edge (%ld, %ld)\n", i+1, j+1);
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
		DebugPrint("Call: Copying Edges of %ld to %ld\n", src+1, target+1);
		if (is_v_deleted[src] || is_v_deleted[target]) {
			WarnPrint("ERROR: Operating On a Deleted Vertex\n");
			return;
		}

		// this block is responsible for making sure the set of edges out of a source, will be copied to target
		// currently its broken for digraphs since it reverse the direction
		size_t n_refs = adj_lists[src].size();
		for (size_t i = 0; i < n_refs; i++) {
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
		DebugPrint("Call: Removing Vertex %ld\n", i+1);
		if (is_v_deleted[i]) {
			return;
		}
		size_t n_out_degree = adj_lists[i].size();

		// remove references to the removed vertex from other lists.
		// in effect this deletes edges of the vertex i before deleting it.
		// Erase (Remove) line t == n_edges_of_i
		// For loops with the erase_remove line is going to be linear in n of edges, t(m)
		// the loop iterates exactly in_degree of (i) times.
		for (size_t x = 0; x < n_out_degree; x++) {
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
		DebugPrint("Call: Contract (remainder: %ld, removed: %ld)\n", i+1, j+1);
		if (i == j) {
			WarnPrint("ERROR: Same Vertex\n");
			return;
		}
		if (g->is_v_deleted[i] || g->is_v_deleted[j]) {
			WarnPrint("ERROR: Contracting on a Deleted Vertex\n");
			return;
		}

		DebugPrint("\t Contracting on Edge (%ld, %ld) (current edges: %ld)\n", i+1, j+1, g->n_edges());

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
			for (size_t k = 0; k < n; k++) {
				reorder[g->visitation_order[k]-1] = k; 
			}
		}

		for (size_t i=0; i < n; i++) {
			size_t vertex = initial_pass ? i : reorder[i];
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
		DebugPrint("DFS on %ld\n", vertex);
		scc_sets[vertex] = scc_current_leader;	// coloring strongly connected components with the same color.
		// reserve mem for stack
		unexplored[vertex] = false;

		size_t n_edges = adj_lists[vertex].size(); // the number of edges of in v.
		for (size_t i = 0; i < n_edges; i++) {
			size_t new_v = adj_lists[vertex][i];
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

		for (size_t i = 0; i < n_vertices; i++) {

			n_edges = adj_lists[i].size();

			for (size_t j = 0; j < n_edges; j++) {
				auto target = adj_lists[i][j];
				new_edges[target].push_back(i); 
			}
		}

		this->adj_lists = new_edges;
	}


	/*
	Should return a Map sorted by weights (keys), with 
	@return std::vector<size_t> of all the 'destination' / 'adjacent' wathave u vertices.
	*/
	static std::vector<std::tuple<size_t,std::vector<size_t>, std::vector<size_t>>> GetAdjList_Subgraph(const Digraph *g, const std::vector<size_t> vertices) {
		size_t adj_size, vs;
		std::vector<std::tuple<size_t,std::vector<size_t>, std::vector<size_t>>> results;
		std::vector<size_t> c_sources;
		std::vector<size_t> c_adj_list;
		std::vector<size_t> c_weights;
		std::vector<size_t> v, w;

		results.reserve(vertices.size());

		for (size_t i=0; i<vertices.size() ; i++) {
			vs = vertices[i];
			c_adj_list = g->adj_lists[vs];
			c_weights = g->weights[vs];

			adj_size = c_adj_list.size();
			for (size_t j=0; j< adj_size; j++) {
				v.push_back(c_adj_list[j]);
				w.push_back(c_weights[j]);
			}
			results[i] = std::make_tuple(vs, v, w);
		}

		return results;
	}

	static std::pair<size_t,size_t> GetClosestVertex(const Digraph *g, std::vector<bool> *visited, size_t *vp, std::vector<size_t> *_distances) {
		size_t min_v;
		size_t min = std::numeric_limits<size_t>::max();
		size_t dist;
		size_t adj_v;
		for (size_t j = 0; j < g->n; j++) {
			if (visited->at(j) == true) {
				DebugPrint("GetClosestVertex of V%lu \n", j + 1);
				for (size_t i = 0; i < g->adj_lists[j].size(); i++) {
					dist = g->weights[j][i] + _distances->at(j);
					adj_v = g->adj_lists[j][i];
					DebugPrint("\tChecking %lu \n", adj_v + 1);
					if (visited->at(adj_v) == false && dist < min) {
						min_v = adj_v;
						min = dist;
						*vp = j;
						DebugPrint("\t\tVmin=%lu | w=%lu | vp=%lu \n", min_v + 1, dist, *vp + 1);
					} else {
						DebugPrint("\t\tSkipping this I %lu\n", i+1);
					}
				}
			} else {
				DebugPrint("\tSkipping this J %lu\n", j+1);
			}
		}

		return std::make_pair(min_v, min);
	}
	
	static std::pair<std::vector<size_t>, std::vector<size_t>> Dijkstra(const Digraph *g, const size_t v0) {
		size_t n = g->n;
		size_t m = g->m;
		size_t infinity = std::numeric_limits<size_t>::max();
		size_t nearest_v, nearest_w, i, j, ith_number_edges, vj, vp, min_dist, new_dist;

		std::vector<size_t> ith_edges;
		std::vector<size_t> ith_weights;

		// Initialize shorted distances ary _distances, and Parent Ary
		std::vector<size_t> _distances = std::vector<size_t>(n, infinity);
		std::vector<size_t> _parents = std::vector<size_t>(n);
		std::vector<bool> _visited = std::vector<bool>(n, false);
		std::pair<size_t, size_t> nearest;
		// define tree (of vertices with computed dist_min)
		// define abyss/unexplored/fringe This ideally should be a priority queue, so it can give u the smallest edge.
		size_t c=0;
		std::fill(_parents.begin(), _parents.end(), c++);
		
		_distances[v0] = 0;
		_parents[v0] = 0; // it should be NULL, but NULL is zero, and zero is a valid vertex. but infinity is not.
		_visited[v0] = true;
		vp = v0;
		nearest_v = v0;
		while (std::find(_visited.begin(), _visited.end(), false) != _visited.end()) {
			min_dist = infinity;
			// Look at all the possible paths from v, and update those whose newer cost is lower
			// if _distances[vj] < weight[v][vj] + _)distance
			
			// Find closest vertex such that it is unvisited
			nearest = GetClosestVertex(g, &_visited, &vp, &_distances);
			nearest_v = nearest.first;
			nearest_w = nearest.second;
			_visited[nearest_v] = true;

			DebugPrint("Standing at V%ld, nearest V is V%ld %ld km away\n", vp, nearest_v, nearest_w);
			DebugPrint("Looking at Adjacent vertex V%ld\n", vp + 1);
			// new_dist = distance_to_here + distance_to_new_v;
			_distances[nearest_v] = nearest_w;
			_parents[nearest_v] = vp;
			DebugPrint("==> SPD to V%lu =%u | from: V%lu  \n", nearest_v + 1, _distances[nearest_v], vp + 1);
		}
		// Return
		return std::make_pair(_distances, _parents);
	}

	static Digraph ParseFileToAdjacency(const std::string filename) {
		bool skip_first_entry = true;
		size_t no_lines = 1;
		size_t line_len = 0;
		size_t _t;
		FILE* pFile = NULL;
		char* line = NULL;
		char* input_char_buffer = NULL;
		std::vector<std::vector<size_t>> adj_lists;
		std::vector<size_t> t_list;

		pFile = fopen(filename.c_str(), "r");
		if (pFile == NULL) { return Digraph(); }

		printf("* reading file: %s\n", filename.c_str());

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

		fclose(pFile);
		free(line);
		free(input_char_buffer);
		Digraph g = Digraph(adj_lists);
		return g;
	}

	static Digraph ParseFileToAdjacencyWeighted(const std::string filename) {
		printf("Call: ParseFileWeighted");
		bool skip_first_entry = true;
		size_t no_lines = 1;
		size_t line_len = 0;
		size_t bsize = 0;
		size_t i = 0;
		size_t _t = 0;
		size_t idx_comma = 0;

		FILE* pFile = NULL;
		char* line = NULL;
		char* input_char_buffer = NULL;
		std::string fly;

		std::vector<std::vector<size_t>> adj_lists;
		std::vector<std::vector<size_t>> adj_lists_weights;
		std::vector<size_t> t_list;
		std::vector<size_t> t_weights;

		pFile = fopen(filename.c_str(), "r");
		if (pFile == NULL) { return Digraph(); }

		printf("* reading file: %s\n",filename.c_str());

		while (getline(&line, &line_len, pFile) != -1)
		{
			skip_first_entry = true;
			input_char_buffer = strtok(line, " \t");
			if (input_char_buffer == NULL) { 
				ErrorPrint("Cannot read your file!\n");
				fclose(pFile);
			}

			while (input_char_buffer != NULL) {
				if ((char)*input_char_buffer == '\r') {
					//printf(" FUCK U DOCTOR STOP MAKINBG WEIRD FILES!\n");
					break;
				}
				fly = std::string(input_char_buffer);
				//printf(" SUBLOOP %d (Buffer = %s | %#016x):\n",i, input_char_buffer, (char) *input_char_buffer );
				if (skip_first_entry) { 
					skip_first_entry = false;
				}
				else 
				{

					bsize = sizeof(input_char_buffer);

					for (size_t _k = 0; _k<bsize; _k++) {
						if (input_char_buffer[_k]==',') {
							idx_comma = _k;
							break;
						}
					}

					_t = (size_t) std::stoi(fly.substr(0,idx_comma));
					t_list.push_back(_t - 1);

					_t = (size_t) std::stoi(fly.substr(idx_comma+1));
					t_weights.push_back(_t);
				}
				input_char_buffer = strtok(NULL, "  \t");
				i++;
			}

			adj_lists.push_back(t_list);
			adj_lists_weights.push_back(t_weights);
			t_list.clear();
			t_weights.clear();
			no_lines ++;
		}

		fclose(pFile);
		free(line);
		
		printf("* Read (%ld) lines\n", no_lines-1);
		Digraph g = Digraph(adj_lists, adj_lists_weights);

		return g;
	}
};
