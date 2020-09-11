class UnionFindElement:
    def __init__(self, label):
        self.label = label
        self.parent = self
        self.rank = 0

    def Label(self):
        return self.label

    def Parent(self):
        return self.parent

    def Rank(self):
        return self.rank

def Find(element):
    if (element.parent != element):
        element.parent = Find(element.parent)
    return element.parent

def Union(element_one, element_two):
    root_one = Find(element_one)
    root_two = Find(element_two)

    if (root_one == root_two): return

    if (root_one.rank < root_two.rank):
        root_one.parent = root_two
    elif (root_one.rank > root_two.rank):
        root_two.parent = root_one
    else:
        root_two.parent = root_one
        root_one.rank = root_one.rank + 1





import itertools



from addax.kavosh.enumerate import EnumerateSubgraphsSequentially
from addax.utilities.dataIO import ReadGraph



def TestSubgraphEnumeration(filename):
    """
    filename: input file to both enumerate via Kavosh and brute force algorithm

    Test subgraph enumeration of the kavosh algorithm against brute force
    """

    # read the graph and number of nodes
    graph = ReadGraph(filename)
    nodes = graph.vertices
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
