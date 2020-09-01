import struct



from addax.data_structures.graph import Graph



def ReadGraph(filename):
    # create empty arrays for nodes, edges, and directed orientation
    nodes = []
    edges = []
    directed = False

    with open(filename, 'rb') as fd:
        # read all of the nodes
        nnodes, = struct.unpack('q', fd.read(8))
        for _ in range(nnodes):
            node, = struct.unpack('q', fd.read(8))
            nodes.append(node)

        # read all of the edge tuples
        nedges, = struct.unpack('q', fd.read(8))
        for _ in range(nedges):
            node_one, node_two,  = struct.unpack('qq', fd.read(16))
            edges.append((node_one, node_two))

        directed, = struct.unpack('q', fd.read(8))

    graph = Graph(nodes, edges, directed)

    return graph



def WriteGraph(graph, filename):
    nodes = graph.nodes
    edges = graph.edges
    directed = graph.directed

    nnodes = len(nodes)
    nedges = len(edges)

    with open(filename, 'wb') as fd:
        # write all of the nodes
        fd.write(struct.pack('q', nnodes))
        fd.write(struct.pack('%sq' % nnodes, *nodes))

        # write all of the edge tuples
        fd.write(struct.pack('q', nedges))
        for (node_one, node_two) in edges:
            fd.write(struct.pack('qq', node_one, node_two))

        # write the directed orientation
        fd.write(struct.pack('q', directed))
