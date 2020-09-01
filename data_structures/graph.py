class Graph(object):
    def __init__(self, nodes, edges, directed):
        # for directed graphs nodes (node_one, node_two) in edges correspond to an edge
        # from node one to node two
        self.nodes = sorted(nodes)
        self.edges = edges
        self.directed = directed
        self.adjacent_edges = {}
        self.outgoing_directed_edges = {}
        self.ingoing_directed_edges = {}

        # create adjacency lists for all nodes
        for node in self.nodes:
            self.adjacent_edges[node] = set()
            self.outgoing_directed_edges[node] = set()
            self.ingoing_directed_edges[node] = set()

        # edges correspond to a tuple
        for (node_one, node_two) in edges:
            # make sure that both nodes exist in the nodes list
            assert (node_one in self.nodes)
            assert (node_two in self.nodes)

            # add the opposite node to the adjacency lists
            self.adjacent_edges[node_one].add(node_two)
            self.adjacent_edges[node_two].add(node_one)

            if directed:
                self.outgoing_directed_edges[node_one].add(node_two)
                self.ingoing_directed_edges[node_two].add(node_one)
            else:
                # if the graph is not directed, then node one and node two have both in and outgoing edges
                self.outgoing_directed_edges[node_one].add(node_two)
                self.outgoing_directed_edges[node_two].add(node_one)
                self.ingoing_directed_edges[node_one].add(node_two)
                self.ingoing_directed_edges[node_two].add(node_one)

        # create ordered lists for the adjacent edges dictionary
        for node in self.adjacent_edges.keys():
            self.adjacent_edges[node] = sorted(list(self.adjacent_edges[node]))
            self.outgoing_directed_edges[node] = sorted(list(self.outgoing_directed_edges[node]))
            self.ingoing_directed_edges[node] = sorted(list(self.ingoing_directed_edges[node]))



    def NEdges(self):
        return len(self.edges)



    def KthNode(self, k):
        return self.nodes[k]



    def NNodes(self):
        return len(self.nodes)



    def EdgesAdjacentToKthNode(self, k):
        return self.adjacent_edges[k]



    def EdgesFromKthNode(self, k):
        return self.outgoing_directed_edges[k]



    def EdgesToKthNode(self, k):
        return self.ingoing_directed_edges[k]



    def DirectedEdgesFromK(self, k):
        return self.outward_
