import random



from addax.data_structures.graph import Graph
from addax.utilities.dataIO import WriteGraph
from addax.visualize.graph import VisualizeGraph



def GenerateSimpleRandomGraphs(N, nnodes, nedges, directed):
    """
    Generate a simple random graph that randomly chooses edges to connect vertices
    @param N: the number of random graphs to generate
    @param nnodes: the number of nodes in the random graph
    @param nedges: the number of edges in the random graph
    @param directed: are the edges directed or not
    """

    for ii in range(N):
        # create an array of vertex indices (connectomes have non consecutive indices)
        vertex_indices = random.sample(range(2**31, 2**32), nnodes)

        # for directed graphs, all combinations of two nodes, outside of self loops, are valid edges
        if directed: edges = [(vertex_one, vertex_two) for vertex_one in vertex_indices for vertex_two in vertex_indices if not vertex_one == vertex_two]
        # for undirected graphs, unique combinations of two nodes are valid (non-parallel) edges
        else: edges = [(vertex_indices[iv1], vertex_indices[iv2]) for iv1 in range(nnodes) for iv2 in range(iv1 + 1, nnodes)]

        # take the first nedges of the edges array to add to the graph
        random.shuffle(edges)
        edges = edges[:nedges]

        # create a new graph object
        graph = Graph(directed)

        # add the vertex indices (there are no communities in this simple graph)
        for vertex_index in vertex_indices:
            graph.AddVertex(vertex_index)

        # add the edges for this graph
        for edge in edges:
            # there are no weights in this simple graph
            graph.AddEdge(edge[0], edge[1])

        # write the simple graph to file
        if directed: output_filename = 'random/simple/graphs/nnodes-{}-nedges-{}-directed-{:05d}.graph'.format(nnodes, nedges, ii)
        else: output_filename = 'random/simple/graphs/nnodes-{}-nedges-{}-undirected-{:05d}.graph'.format(nnodes, nedges, ii)

        WriteGraph(graph, output_filename)

        # visualize the graph
        if directed: output_filename = 'random/simple/visuals/nnodes-{}-nedges-{}-directed-{:05d}.dot'.format(nnodes, nedges, ii)
        else: output_filename = 'random/simple/visuals/nnodes-{}-nedges-{}-undirected-{:05d}.dot'.format(nnodes, nedges, ii)

        VisualizeGraph(graph, output_filename)
