import networkx as nx



def VisualizeGraph(graph, output_filename):
    nodes = graph.nodes
    edges = graph.edges
    directed = graph.directed

    if directed:
        viz_graph = nx.DiGraph()
    else:
        viz_graph = nx.Graph()

    node_mapping = {}
    for index, node in enumerate(nodes):
        viz_graph.add_node(index, label=str(node))

        # create a mapping from the node number to the index
        node_mapping[node] = index

    for (node_one, node_two) in edges:
        node_one_index = node_mapping[node_one]
        node_two_index = node_mapping[node_two]

        viz_graph.add_edge(node_one_index, node_two_index)

    A = nx.nx_agraph.to_agraph(viz_graph)
    A.layout(prog='dot')
    A.draw(output_filename)
