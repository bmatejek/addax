#include "cpp-graph.h"



Vertex::Vertex(Graph *input_graph, long input_index, int input_community) :
graph(graph),
index(index),
community(community)
{
    graph = input_graph;
    index = input_index;
    community = input_community;

    // extra instance variables keep track of the ingoing and outgoing edges from the vertex
    incoming_edges = std::vector<Edge *>();
    outgoing_edges = std::vector<Edge *>();
    // keep track of incoming and outgoing neighbors
    incoming_neighbors = std::unordered_set<long>();
    outgoing_neighbors = std::unordered_set<long>();
    neighbors = std::unordered_set<long>();
}

Vertex::~Vertex(void)
{
    /*
    Destructor for eliminating vertex object
    */
}

void Vertex::AddEdge(Edge *edge)
{
    /*
    Add this edge to the set of edges for this vertex and ensure no edge parallelism

    @param edge: the edge that connects this vertex to another
    */

    // assert that this is a valid edge for this vertex
    assert (edge->source_index == index || edge->destination_index == index);

    // if the graph is directed, add the incoming or outgoing edge
    if (graph->directed) {
        if (edge->source_index == index) {
            outgoing_edges.push_back(edge);
            assert (outgoing_neighbors.find(edge->destination_index) == outgoing_neighbors.end());
            outgoing_neighbors.insert(edge->destination_index);
            neighbors.insert(edge->destination_index);
        }
        else {
            incoming_edges.push_back(edge);
            assert (incoming_neighbors.find(edge->source_index) == incoming_neighbors.end());
            incoming_neighbors.insert(edge->source_index);
            neighbors.insert(edge->source_index);
        }
    }
    // if the graph is not directed, add the edge to both incoming and outgoing
    else {
        incoming_edges.push_back(edge);
        outgoing_edges.push_back(edge);

        if (edge->source_index == index) {
            assert (incoming_neighbors.find(edge->destination_index) == incoming_neighbors.end());
            assert (outgoing_neighbors.find(edge->destination_index) == outgoing_neighbors.end());
            incoming_neighbors.insert(edge->destination_index);
            outgoing_neighbors.insert(edge->destination_index);
            neighbors.insert(edge->destination_index);
        }
        else {
            assert (incoming_neighbors.find(edge->source_index) == incoming_neighbors.end());
            assert (outgoing_neighbors.find(edge->source_index) == outgoing_neighbors.end());
            incoming_neighbors.insert(edge->source_index);
            outgoing_neighbors.insert(edge->source_index);
            neighbors.insert(edge->source_index);
        }
    }
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

    vertices = std::map<long, Vertex *>();
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

long Graph::NVertices(void)
{
    /*
    Return the number of vertices in this graph
    */
    return vertices.size();
}

long Graph::NEdges(void)
{
    /*
    Return the number of edges in this graph
    */
    return edges.size();
}




Graph *ReadGraph(const char input_filename[4096])
{
    FILE *fp = fopen(input_filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to open %s\n", input_filename); return NULL; }

    // read the basic attributes for the graph
    long nvertices, nedges;
    bool directed;
    if (fread(&nvertices, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return NULL; }
    if (fread(&nedges, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return NULL; }
    if (fread(&directed, sizeof(bool), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return NULL; }

    // read the prefix
    char prefix[128];
    if (fread(&prefix, sizeof(char), 128, fp) != 128) { fprintf(stderr, "Failed to read %s\n", input_filename); return NULL; }

    Graph *graph = new Graph(prefix, directed);

    // read all the vertices and add them to the graph
    for (long iv = 0; iv < nvertices; ++iv) {
        long index, community;
        if (fread(&index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return NULL; }
        if (fread(&community, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return NULL; }

        graph->AddVertex(index, community);
    }

    // read all of the edges and add them to the graph
    for (long ie = 0; ie < nedges; ++ie) {
        long source_index, destination_index;
        float weight;
        if (fread(&source_index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return NULL; }
        if (fread(&destination_index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return NULL; }
        if (fread(&weight, sizeof(double), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", input_filename); return NULL; }

        graph->AddEdge(source_index, destination_index, weight);
    }

    return graph;
}
