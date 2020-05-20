# Read Me

The LeSs ShiTty GraPh Program is a fun program!

To build the small software just light the candle with the following command:
`./candle.sh`
If your lighter is broken, and it gives you a permission denied then easily

`chmod 555 candle.sh` before lighting it.


The Code is immflammable and flammable, use at your own risk.
## Implemented Algorithms
### Undirected Graphs:
*

### Directed Graphs:
* 

## Todos

### Algorithms
#### Undirected Graphs:
* Breadth First Search

#### Directed Graphs:
* Breadth First Search
* Dijkstra Shortest Path

### Others
* Review ugraph and digraph and mark obsolete methods and methods untested on the class. use @todo and @deprecated
    - If the fix requires more than 10 mins do it later.
    - If the thing makes no sense (not well defined mathematically) delete it.
* Measure the performance of the subroutines in time and memory. The How is not clear, plan it.
* (optional) Graph said measurements 
* select a trial easy method to parallelize using openmpi or mpich
    1. the loop over karger mincut is the easier (start many kmincuts on different cpus)
* abstract away (if its efficient enough) from digraph and ugraph to graph

## Building the code
* I will write the make file a bit later.
but for now, lets see all the steps


```gcc -E bin/prprocessor_out.e -S bin/asm.a -C bin/comp.c -o bin/out.o```



graph
directed graph
tree
weighted directed graph
weighted graph..

// assume this is the abstract class
set<vertices> inbound_vertices(v) {
    adjlist[v];
}

set<vertices> outbound_vertices(v) {
    reverse look up
}

auto unknown_graph_type = ReadGraph(file.txt);

Tree = ReadGraph(preprocessing)