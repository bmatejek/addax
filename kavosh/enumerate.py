import sys
import time
import itertools



def Validate(G, parents, u, visited):
    """
    Find the valid vertices for the next recursive level. Excludes elements with indices smaller
    than the root vertex and those already visited at a previous recursive level

    @param G: graph
    @param parents: selected vertices of last layer
    @param u: root vertex
    visited: a list of vertices already visited

    Returns a sorted list of valid vertices for the next level of analysis
    """
    # create a set of valid vertices at this level
    valid_vertices = set()

    # iterate over all the parents of this layer
    for v in parents:
        # iterate over all the neighbors of this parent
        for w in G.vertices[v].NeighborIndices():
            # if the root vertex is less than the neighbor and the neighbor has not been visited
            if u < w and not w in visited:
                visited.add(w)

                valid_vertices.add(w)

    return sorted(list(valid_vertices))



def EnumerateVertex(G, u, S, rem, i, visited):
    """
    Enumerate all subgraphs of size rem that contain vertices in S[0] ... S[i - 1]

    @param G: graph
    @param u: root vertex
    @param S: selection (S = {S_0, S_i, ... S_{k - 1}}) is an array of the set of all S_i
    @param rem: number of remaining vertices to be selected
    @param i: current depth of the tree
    @param visited: a list of vertices already visited

    Returns a generator that continually gives the next subgraph that contains S[0] ... S[i - 1] of
    the appropriate size (k)
    """
    # if there are no remaining vertices to add, subgraph size limit reached
    if not rem:
        # the set of vertices in S[0] to S[i - 1] contain the subgraph
        # note that the sets in S[i] ... S[k] do not belong to the subgraph but a previous iteration
        enumerated_subgraph = set()
        for index in range(i):
            enumerated_subgraph = enumerated_subgraph.union(S[index])

        # return a sorted tuple for this subgraph
        yield tuple(sorted(enumerated_subgraph))

    # get the valid vertices from the other
    valid_vertices = Validate(G, S[i - 1], u, visited)

    # the max number of vertices for this layer is the minimum of the number of children or the remaining
    n_i = min(len(valid_vertices), rem)

    # iterate over the total number of possible vertices at this layer (follows number compositions)
    for k_i in range(1, n_i + 1):
        # get all combinations of k for the valid vertices
        combinations = itertools.combinations(valid_vertices, k_i)

        # try all different combinations for this size k
        for combination in combinations:
            S[i] = set(combination)

            # enumerate all subgraphs of size rem - k_i that contain all vertices in S[0] ... S[i]
            yield from EnumerateVertex(G, u, S, rem - k_i, i + 1, visited)

    # remove all of the valid vertices from the list of visited
    # finished all yields (subgraphs) for this level and proceed back to level above (closer to root)
    for v in valid_vertices:
        visited.remove(v)



def EnumerateSubgraphsFromNode(G, k, u):
    """
    Enumerate all subgraphs of a given motif size rooted at a given vertex

    @param G: graph
    @param k: motif size
    @param u: root vertex index

    Returns a generator that continually gives the next subgraph in the graph rooted as this vertex
    """
    # make sure this vertex appears in the graph
    assert (u in G.vertices.keys())

    # keep track globally (through parameter passing) the vertices visited at higher enumeration steps
    visited = set()

    # add the root vertex to the visited set
    visited.add(u)

    # create the selection (first layer only ever has the root vertex)
    S = {}
    S[0] = set()
    S[0].add(u)

    # enumerate all subgraphs of size k - 1 that contain the root u
    yield from EnumerateVertex(G, u, S, k - 1, 1, visited)



def EnumerateSubgraphsSequentially(G, k):
    """
    Enumerate all subgraphs of a given size in the graph

    @param G: graph
    @param k: motif size

    Returns a generator that continually gives the next subgraph in the graph
    """

    # iterate over all vertices in the graph
    for u in G.vertices.keys():
        # enumerate all subgraphs of size k that contain the root u
        yield from EnumerateSubgraphsFromNode(G, k, u)
