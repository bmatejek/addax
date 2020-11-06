#ifndef __CPP_GRAPH_H__
#define __cpp_graph_H__

#include <assert.h>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <unordered_set>
#include <unordered_map>



class Vertex;
class Edge;
class Graph;



class Vertex {
public:
    // constructors/destructors
    Vertex(Graph *graph, long index, int community);
    ~Vertex();

    // modifying functions
    void AddEdge(Edge *edge);

private:
    // instance variables
    Graph *graph;
    long index;
    int community;

    std::vector<Edge *> incoming_edges;
    std::vector<Edge *> outgoing_edges;
    std::unordered_set<Edge *> incoming_neighbors;
    std::unordered_set<Edge *> outgoing_neighbors;
    std::unordered_set<Edge *> neighbors;
};



class Edge {
public:
    // constructors/destructors
    Edge(Graph *graph, long source_index, long destination_index, float weight);
    ~Edge();

private:
    // instance variables
    Graph *graph;
    long source_index;
    long destination_index;
    float weight;
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
    Graph(const char input_prefix[128], bool input_directed);
    ~Graph();

    // modifying functions
    void AddVertex(long index, int community = -1);
    void AddEdge(long source_index, long destination_index, float weight = -1);

private:
    char prefix[128];
    bool directed;
    std::unordered_map<long, Vertex *> vertices;
    std::vector<Edge *> edges;
    std::unordered_set<std::pair<long, long>, edge_pair_hash> edge_set;
};






#endif
