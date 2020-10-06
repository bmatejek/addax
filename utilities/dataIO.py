import struct



from addax.data_structures.graph import Graph



def ReadGraph(input_filename):
    """
    Read a graph data structure from disk

    @param input_filename: the filename where the graph data is stored
    """
    with open(input_filename, 'rb') as fd:
        # read the basic attributes for the graph
        nvertices, nedges, directed, = struct.unpack('qq?', fd.read(17))

        # read the prefix
        prefix, = struct.unpack('128s', fd.read(128))
        prefix = prefix.decode().strip('\0')

        graph = Graph(prefix, directed)

        # read all the vertices and add them to the graph
        for _ in range(nvertices):
            index, community, = struct.unpack('qq', fd.read(16))

            graph.AddVertex(index, community)

        # read all of the edges and add them to the graph
        for _ in range(nedges):
            source_index, destination_index, weight, = struct.unpack('qqd', fd.read(24))

            graph.AddEdge(source_index, destination_index, weight)

        return graph



def WriteGraph(graph, output_filename):
    """
    Write a graph to disk for later I/O access

    @param graph: the graph data structure to save to disk
    @param output_filename: the location to save the graph data structure
    """
    with open(output_filename, 'wb') as fd:
        # write the basic attributes for the graph to disk
        nvertices = graph.NVertices()
        nedges = graph.NEdges()
        directed = graph.directed
        prefix = graph.prefix

        fd.write(struct.pack('qq?', nvertices, nedges, directed))
        fd.write(struct.pack('128s', prefix.encode()))

        # write all of the vertices and their attributes
        for vertex in graph.vertices.values():
            fd.write(struct.pack('qq', vertex.index, vertex.community))

        # write all of the edges and their attributes
        for edge in graph.edges:
            fd.write(struct.pack('qqd', edge.source_index, edge.destination_index, edge.weight))
