#ifndef __CPP_GRAPH_H__
#define __cpp_graph_H__

#include <assert.h>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <unordered_set>
#include <map>
#include <bzlib.h>



class Vertex;
class Edge;
class Graph;



// I/O functions for raeding a graph
Graph *ReadGraph(const char input_filename[4096]);
Graph *ReadBZ2Graph(const char input_filename[4096]);



class Vertex {
public:
    // constructors/destructors
    Vertex(Graph *graph, long index, long enumeration_index, long community, long color);
    ~Vertex();

    // modifying functions
    void AddEdge(Edge *edge);

    // instance variables
    Graph *graph;
    long index;
    long enumeration_index;
    long community;
    long color;


    // extra instance variables keep track of the ingoing and outgoing edges from the vertex
    std::vector<Edge *> incoming_edges;
    std::vector<Edge *> outgoing_edges;
    // keep track of incoming and outgoing neighbors
    std::unordered_set<long> incoming_neighbors;
    std::unordered_set<long> outgoing_neighbors;
    std::unordered_set<long> neighbors;
};



class Edge {
public:
    // constructors/destructors
    Edge(Graph *graph, long source_index, long destination_index, double weight);
    ~Edge();

    // instance variables
    Graph *graph;
    long source_index;
    long destination_index;
    double weight;
};



struct edge_pair_hash
{	std::size_t operator() (const std::pair<long, long> &pair) const
	{
		return std::hash<long>()(pair.first) ^ std::hash<long>()(pair.second);
	}
};



class Graph {
public:
    // constructors/destructors
    Graph(const char input_prefix[128], bool input_directed, bool input_colored);
    ~Graph();

    // modifying functions
    void AddVertex(long index, long enumeration_index, long community = -1, long color = -1);
    void AddEdge(long source_index, long destination_index, double weight = -1);

    // attribute functions
    long NVertices(void);
    long NEdges(void);

    // instance variables
    char prefix[128];
    bool directed;
    bool colored;
    std::map<long, Vertex *> vertices;
    std::vector<Edge *> edges;
    std::unordered_set<std::pair<long, long>, edge_pair_hash> edge_set;
};



#endif
