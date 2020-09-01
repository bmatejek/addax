import itertools



from addax.data_structures.unionfind import UnionFindElement, Find, Union
from addax.kavosh.enumerate import EnumerateSubgraphsSequentially
from addax.utilities.dataIO import ReadGraph



def TestSubgraphEnumeration(filename):
    """
    filename: input file to both enumerate via Kavosh and brute force algorithm

    Test subgraph enumeration of the kavosh algorithm against brute force
    """

    # read the graph and number of nodes
    graph = ReadGraph(filename)
    nodes = graph.nodes
    nnodes = len(nodes)

    # consider ks from 2 to the number of nodes
    ks = range(2, nnodes)

    for k in ks:
        subgraphs = EnumerateSubgraphsSequentially(graph, k)

        # make sure all subgraphs are unique
        unique_subgraphs = set()
        for subgraph in subgraphs:
            assert (not subgraph in unique_subgraphs)
            unique_subgraphs.add(subgraph)

        # consider all combinations of subgraphs
        nbrute_force_combinations = 0
        for combination in itertools.combinations(nodes, k):
            vertices = {}

            # create a union find element for all nodes
            for node in combination:
                vertices[node] = UnionFindElement(node)

            # connect nodes that share edges within the subgraph
            for node in combination:
                for neighbor_node in graph.EdgesAdjacentToKthNode(node):
                    # do not consider edges that leave the combination subgraph
                    if not neighbor_node in combination: continue

                    Union(vertices[node], vertices[neighbor_node])

            # get the parent nodes of the set of vertices
            parent_nodes = set()
            for node in combination:
                parent = Find(vertices[node])
                parent_nodes.add(parent.label)

            # this is only a valid subgraph is there is one parent
            if len(parent_nodes) == 1:
                nbrute_force_combinations += 1

        assert (len(subgraphs) == nbrute_force_combinations)
