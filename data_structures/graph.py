import struct



import networkx as nx



import community as community_louvain



class Graph(object):
    def __init__(self, prefix, directed):
        """
        Graph class defines the basic graph structure for addax used for clustering commmunities, motif discovery,
        and generating random examples

        @param directed: indicates if the graph is directed or undirected
        """
        self.prefix = prefix
        self.directed = directed

        # vertices is a mapping from the vertex index to the vertex object
        self.vertices = {}
        # edges is a list of edges with sources, destinations, and weights
        self.edges = []
        # the edge set contains a list of (source, destination) indices
        self.edge_set = set()

    def AddVertex(self, index, community = -1):
        """
        Add a vertex to the graph

        @param index: the index for the vertex
        @param community: the community that the vertex belongs to (default = -1)
        """
        # vertices must have unique indices
        assert (not index in self.vertices)

        # create the vertex and add it to the mapping
        vertex = self.Vertex(self, index, community)
        self.vertices[index] = vertex

    def AddEdge(self, source_index, destination_index, weight = 1):
        """
        Add an edge to the graph

        @param source_index: the integer of the source index in the graph
        @param destination_index: the integer of the destination index in the graph
        @param weight: the weight of this edge where higher values indicate greater strength (default = 1)
        """
        # the source and destination indices must actually belong to vertices
        assert (source_index in self.vertices)
        assert (destination_index in self.vertices)

        # do not allow self loops
        assert (not source_index == destination_index)

        # if the graph is undirected, make the source destination the smaller of the two indices
        if not self.directed and destination_index < source_index:
            tmp = destination_index
            destination_index = source_index
            source_index = tmp

        # create the edge and add it to the list of edges
        edge = self.Edge(self, source_index, destination_index, weight)
        self.edges.append(edge)

        # add to the set of edges in the graph for easier look up
        if self.directed:
            self.edge_set.add((source_index, destination_index))
        else:
            # directed edges go in both directions
            self.edge_set.add((source_index, destination_index))
            self.edge_set.add((destination_index, source_index))

        # add the edge to both vertices
        self.vertices[source_index].AddEdge(edge)
        self.vertices[destination_index].AddEdge(edge)

    def NVertices(self):
        """
        Return the number of vertices in this graph
        """
        return len(self.vertices.keys())

    def NEdges(self):
        """
        Return the number of edges in this graph
        """
        return len(self.edges)

    def DetectCommunities(self, output_filename = None):
        """
        Returns a list of communities based on the Louvain algorithm
        """
        # initialize a networkx graph
        G = nx.Graph()

        # add all vertices
        for vertex in self.vertices.values():
            G.add_node(vertex.index)

        # add all edges to the networkx graph
        undirected_edges = {}
        for edge in self.edges:
            # get the min and max edge
            edge_one = min(edge.source_index, edge.destination_index)
            edge_two = max(edge.destination_index, edge.source_index)

            if not (edge_one, edge_two) in undirected_edges:
                undirected_edges[(edge_one, edge_two)] = edge.weight
            else:
                undirected_edges[(edge_one, edge_two)] += edge.weight

        # add the undirected edge to the graph
        for (edge_one, edge_two) in undirected_edges:
            G.add_edge(edge_one, edge_two, weight=undirected_edges[(edge_one, edge_two)])

        # determine communities in the graph
        partition = community_louvain.best_partition(G)

        # write the partition to file
        if not output_filename == None:
            with open(output_filename, 'wb') as fd:
                fd.write(struct.pack('q', self.NVertices()))
                for (neuron_id, community) in partition.items():
                    fd.write(struct.pack('qq', neuron_id, community))

        # get a list of communities using the Louvain algorithm
        return partition

    def Communities(self):
        """
        Return a mapping from vertex indices to communities
        """
        communities = {}

        for vertex in self.vertices.values():
            communities[vertex.index] = vertex.community

        return communities

    class Vertex(object):
        def __init__(self, graph, index, community = -1):
            """
            Vertex class defines the vertices in a graph that are labeled by the index

            @param graph: the larger graph that contains this vertex
            @param index: the integer index that corresponds to this vertex
            @param community: the community that the vertex belongs to (default = -1)
            """
            self.graph = graph
            self.index = index
            self.community = community

            # extra instance variables keep track of the ingoing and outgoing edges from the vertex
            self.incoming_edges = []
            self.outgoing_edges = []
            # keep track of incoming and outgoing neighbors
            self.incoming_neighbors = set()
            self.outgoing_neighbors = set()
            self.neighbors = set()

        def AddEdge(self, edge):
            """
            Add this edge to the set of edges for this vertex and ensure no edge parallelism

            @param edge: the edge that connects this vertex to another
            """
            # ensure that this is a valid edge for this vertex
            assert (edge.source_index == self.index or edge.destination_index == self.index)

            # if the graph is directed, add the incoming or outgoing edge
            if self.graph.directed:
                if edge.source_index == self.index:
                    self.outgoing_edges.append(edge)
                    assert (not edge.destination_index in self.outgoing_neighbors)
                    self.outgoing_neighbors.add(edge.destination_index)
                    self.neighbors.add(edge.destination_index)
                else:
                    self.incoming_edges.append(edge)
                    assert (not edge.source_index in self.incoming_neighbors)
                    self.incoming_neighbors.add(edge.source_index)
                    self.neighbors.add(edge.source_index)
            # if the graph is not directed, add the edge to both incoming and outgoing
            else:
                self.incoming_edges.append(edge)
                self.outgoing_edges.append(edge)

                if edge.source_index == self.index:
                    assert (not edge.destination_index in self.incoming_neighbors and not edge.destination_index in self.outgoing_neighbors)
                    self.incoming_neighbors.add(edge.destination_index)
                    self.outgoing_neighbors.add(edge.destination_index)
                    self.neighbors.add(edge.destination_index)
                else:
                    assert (not edge.source_index in self.incoming_neighbors and not edge.source_index in self.outgoing_neighbors)
                    self.incoming_neighbors.add(edge.source_index)
                    self.outgoing_neighbors.add(edge.source_index)
                    self.neighbors.add(edge.source_index)

        def IncomingNeighborIndices(self):
            """
            Returns the neighbors with edges going from
            """
            return self.incoming_neighbors

        def OutgoingNeighborIndices(self):
            """
            Returns the neighbors with an edge from this vertex to that neighbor
            """
            return self.outgoing_neighbors

        def NeighborIndices(self):
            """
            Return all neighbors from this vertex regardless of incoming and outgoing status
            """
            return self.neighbors


    class Edge(object):
        def __init__(self, graph, source_index, destiantion_index, weight = 1):
            """
            Edge class defines the edges in a graph that connect the vertices

            @param graph: the larger graph that contains this edge
            @param source_index: the integer of the source index in the graph
            @param destination_index: the integer of the destination index in the graph
            @param weight: the weight of this edge where higher values indicate greater strength (default = 1)
            """
            self.graph = graph
            self.source_index = source_index
            self.destination_index = destiantion_index
            self.weight = weight
