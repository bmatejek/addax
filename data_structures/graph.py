class Graph:
    def self(__init__, nodes, edges, directed):
        # for directed graphs nodes (node_one, node_two) in edges correspond to an edge
        # from node one to node two
        self.nodes = sorted(nodes)
        self.edges = edges
        self.directed = directed
        self.adjacent_edges = {}
        self.outgoing_directed_edges = {}
        self.ingoing_directed_edges = {}

        # edges correspond to a tuple
        for (node_one, node_two) in edges:
            # make sure that both nodes exist in the nodes list
            assert (node_one in self.nodes)
            assert (node_two in self.nodes)

            # create new dictionary entries if these nodes have not yet been seen
            if not node_one in adjacent_edges:
                self.adjacent_edges[node_one] = set()
                self.outgoing_directed_edges[node_one] = set()
                self.ingoing_directed_edges[node_one] = set()
            if not node_two in adjacent_edges:
                self.adjacent_edges[node_two] = set()
                self.outgoing_directed_edges[node_two] = set()
                self.ingoing_directed_edges[node_two] = set()

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



    def PruneGraph(self, min_node):
        """
        Prune the graph so that the minimum label for this graph is min_node.
        Returns a new graph object.
        """
        pruned_nodes = []
        pruned_edges = []

        for node in self.nodes:
            if node < min_node: continue
            pruned_nodes.append(node)

        for (node_one, node_two) in self.edges:
            if node_one < min_node or node_two < min_node: continue
            pruned_edges.append((node_one, node_two))

        return Graph(pruned_nodes, pruned_edges, self.directed)



    def KthNode(self, k):
        return self.nodes[k]



    def EdgesAdjacentToKthNode(self, k):
        return self.adjacent_edges[k]



    def EdgesFromKthNode(self, k):
        return self.outgoing_directed_edges[k]



    def EdgesToKthNode(self, k):
        return self.ingoing_directed_edges[k]



    def DirectedEdgesFromK(self, k):
        return self.outward_
