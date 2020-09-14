import time
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

    # read the graph into memory and get the vertex indices and the number of vertices
    graph = ReadGraph(filename)
    vertices = graph.vertices
    nvertices = graph.NVertices()

    # consider ks from 2 to the number of vertices
    ks = range(2, nvertices)

    kavosh_total_time = 0
    brute_force_total_time = 0

    # iterate over all considered motif sizes
    for k in ks:
        kavosh_time = time.time()

        # create a set of subgraphs discovered by the kavosh enumeration strategy
        subgraphs = set()

        # retrieve all subgraphs in the graph of size k
        for subgraph in EnumerateSubgraphsSequentially(graph, k):
            assert (not subgraph in subgraphs)
            subgraphs.add(subgraph)

        kavosh_total_time += (time.time() - kavosh_time)

        brute_force_time = time.time()

        # consider all combinations of subgraphs
        nbrute_force_subgraphs = 0
        for vertex_combination in itertools.combinations(vertices.keys(), k):
            # create a union find data structure to see if these vertices are connected in the subgraph
            union_find_elements = {}
            for vertex in vertex_combination:
                # create a union find element with this label
                union_find_elements[vertex] = UnionFindElement(vertex)

            # connect vertices that share edges within the subgraph
            for vertex in vertex_combination:
                for neighbor_vertex in graph.vertices[vertex].NeighborIndices():
                    # do not consider any edges that remain outside the considered subgraph
                    if not neighbor_vertex in vertex_combination: continue

                    # link together the vertices that share an edge
                    Union(union_find_elements[vertex], union_find_elements[neighbor_vertex])

            # if all vertcies share the same parent, they belong to the same tree and the subgraph is valid
            union_find_parents = set()
            for vertex in vertex_combination:
                parent = Find(union_find_elements[vertex])
                union_find_parents.add(parent)

            # if there is not one parent, the set of vertices is not connected
            if not len(union_find_parents) == 1: continue

            # this is only a valid subgraph is there is one parent
            nbrute_force_subgraphs += 1

        brute_force_total_time += (time.time() - brute_force_time)

        assert (len(subgraphs) == nbrute_force_subgraphs)

    print ('Verified {}'.format(filename))
    print ('  Kavosh Enumeration Time: {:0.2f} seconds'.format(kavosh_total_time))
    print ('  Brute Force Enumeration Time: {:0.2f} seconds'.format(brute_force_total_time))
