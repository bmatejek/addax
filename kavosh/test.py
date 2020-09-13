import itertools



from addax.data_structures.unionfind import UnionFindElement, Union, Find
from addax.kavosh.enumerate import EnumerateSubgraphsSequentially
from addax.utilities.dataIO import ReadGraph



def TestSubgraphEnumeration(filename):
    """
    Test the subgraph enumeration algorithm defined by the Kavosh algorithm and implemented in
    addax.kavosh.enumerate.py

    @param filename: input file to enumerate all subgraphs via Kavosh and a brute force strategy.
    """

    # read the graph and number of nodes
    graph = ReadGraph(filename)
    nodes = graph.vertices
    nnodes = len(nodes)

    # consider ks from 2 to the number of nodes
    ks = range(2, nnodes)

    for k in ks:
        subgraphs = set()

        for subgraph in EnumerateSubgraphsSequentially(graph, k):
            assert (not subgraph in subgraphs)
            subgraphs.add(subgraph)

        # consider all combinations of subgraphs
        nbrute_force_combinations = 0
        for combination in itertools.combinations(nodes, k):
            vertices = {}

            # create a union find element for all nodes
            for node in combination:
                vertices[node] = UnionFindElement(node)

            # connect nodes that share edges within the subgraph
            for node in combination:
                for neighbor_node in graph.vertices[node].Neighbors():
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
