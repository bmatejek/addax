import sys
import time
import itertools



def Validate(G, parents, u, visited):
    """
    G: graph
    parents: selected vertices of last layer
    u: root vertex
    visited: a list of nodes already visited

    Returns a sorted list of valid vertices for the next level
    """
    # create a set of valid vertices at this level
    valid_vertices = set()

    # iterate over all the parents of this layer
    for v in parents:
        # iterate over all the neighbors of this parent
        for w in G.EdgesAdjacentToKthNode(v):
            # if the root vertex is less than the neighbor and the neighbor has not been visited
            if u < w and not w in visited:
                visited.add(w)

                valid_vertices.add(w)

    return sorted(list(valid_vertices))



def EnumerateVertex(G, u, S, rem, i, visited, subgraphs):
    """
    G: graph
    u: root vertex
    S: selection (S = {S_0, S_i, ... S_{k - 1}}) is an array of the set of all S_i
    rem: number of remaining vertices to be selected
    i: current depth of the tree
    visited: a list of nodes already visited
    """
    # if there are no remaining vertices to add, subgraph size limit reached
    if not rem:
        # the set of vertices in S[0] to S[i - 1] contain the subgraph
        # note that the sets in S[i] ... S[k] do not belong to the subgraph but a previous iteration
        enumerated_subgraph = set()
        for index in range(i):
            enumerated_subgraph = enumerated_subgraph.union(S[index])

        subgraphs.append(tuple(sorted(enumerated_subgraph)))

        return

    # get the valid vertices from the other
    valid_vertices = Validate(G, S[i - 1], u, visited)

    # the max number of nodes for this layer is the minimum of the number of children or the remaining
    n_i = min(len(valid_vertices), rem)

    # iterate over the total number of possible nodes at this layer
    for k_i in range(1, n_i + 1):
        # get all combinations of k for the valid vertices
        combinations = itertools.combinations(valid_vertices, k_i)

        # try all different combinations for this size k
        for combination in combinations:
            S[i] = set(combination)

            # enumerate the vertices from
            EnumerateVertex(G, u, S, rem - k_i, i + 1, visited, subgraphs)

    # remove all of the valid vertices from the list of visited
    for v in valid_vertices:
        visited.remove(v)



def EnumerateSubgraphsSequentially(G, k):
    """
    G: graph
    k: motif size
    """
    nodes = G.nodes

    subgraphs = []

    # iterate over all nodes
    for u in nodes:
        # keep track globally (through parameter passing) of the vertices visited in enumeration
        visited = set()

        visited.add(u)

        # create the selection (first layer only has the root)
        S = {}
        S[0] = set()
        S[0].add(u)

        EnumerateVertex(G, u, S, k - 1, 1, visited, subgraphs)

        visited.remove(u)

    return subgraphs
