#include "cpp-graph.h"



Vertex::Vertex(Graph *graph, long index, int community) :
graph(graph),
index(index),
community(community)
{
    incoming_edges = std::vector<Edge *>();
    outgoing_edges = std::vector<Edge *>();
    incoming_neighbors = std::unordered_set<Edge *>();
    outgoing_neighbors = std::unordered_set<Edge *>();
    neighbors = std::unordered_set<Edge *>();
}



Vertex::~Vertex(void)
{
}



void Vertex::AddEdge(Edge *edge)
{
    assert (false);
}



Edge::Edge(Graph *graph, long source_index, long destination_index, float weight) :
graph(graph),
source_index(source_index),
destination_index(destination_index),
weight(weight)
{
}


Edge::~Edge(void)
{
}


Graph::Graph(const char input_prefix[128], bool input_directed)
{
    /*
    Graph class defines the basic graph structure for addax used for clustering commmunities, motif discovery,
    and generating random examples

    @param prefix: a string to reference this graph by
    @param directed: indicates if the graph is directed or undirected
    */
    strncpy(prefix, input_prefix, 128);
    directed = input_directed;

    vertices = std::unordered_map<long, Vertex *>();
    edges = std::vector<Edge *>();
    edge_set = std::unordered_set<std::pair<long, long>, edge_pair_hash>();
}

Graph::~Graph(void)
{
    /*
    Destructor for eliminating graph object
    */
}

void Graph::AddVertex(long index, int community)
{
    /*
    Add a vertex to the graph

    @param index: the index for the vertex
    @param community: the community that the vertex belongs to (default = -1)
    */

    // vertices must have unique indices
    assert (vertices.find(index) == vertices.end());

    // create the vertex and add it to the mapping
    Vertex *vertex = new Vertex(this, index, community);
    vertices[index] = vertex;
}

void Graph::AddEdge(long source_index, long destination_index, float weight)
{
    /*
    Add an edge to the graph

    @param source_index: the integer of the source index in the graph
    @param destination_index: the integer of the destination index in the graph
    @param weight: the weight of this edge where higher values indicate greater strength (default = 1)
    */

    // the source and destination indices must actually belong to vertices
    assert (vertices.find(source_index) != vertices.end());
    assert (vertices.find(destination_index) != vertices.end());

    // do not allow self loops
    assert (source_index != destination_index);

    // if the graph is undirected, make the source destination the smaller of the two indices
    if (!directed && destination_index < source_index) {
        long tmp = destination_index;
        destination_index = source_index;
        source_index = tmp;
    }

    // create the edge and add it to the list of edges
    Edge *edge = new Edge(this, source_index, destination_index, weight);
    edges.push_back(edge);

    // add to the set of edges in the graph for easier look up
    if (directed) {
        edge_set.insert(std::pair<long, long>(source_index, destination_index));
    }
    else {
        // directed edges go in both directions
        edge_set.insert(std::pair<long, long>(source_index, destination_index));
        edge_set.insert(std::pair<long, long>(destination_index, source_index));
    }

    // add the edge to both vertices
    vertices[source_index]->AddEdge(edge);
    vertices[destination_index]->AddEdge(edge);
}
