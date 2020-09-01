import random



from addax.data_structures.graph import Graph
from addax.utilities.dataIO import ReadGraph, WriteGraph
from addax.visualize.graph import VisualizeGraph


def GenerateSimpleRandomGraph(nnodes, nedges, directed, prefix):
    # get a list of node integers (0 to n - 1)
    nodes = range(nnodes)

    edge_list = set()

    while len(edge_list) < nedges:
        node_one = random.choice(nodes)
        node_two = random.choice(nodes)

        # don't allow self loops
        if node_one == node_two: continue

        # don't allow parallel loops
        if (node_one, node_two) in edge_list: continue
        if not directed and (node_two, node_one) in edge_list: continue

        edge_list.add((node_one, node_two))

    # create a graph object
    graph = Graph(nodes, edge_list, directed)

    # output the graph to file
    output_filename = 'graphs/{}.graph'.format(prefix)
    WriteGraph(graph, output_filename)

    # visualize the output graph
    output_filename = 'visuals/{}.dot'.format(prefix)
    VisualizeGraph(graph, output_filename)
