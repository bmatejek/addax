import sys
import pynauty



from addax.data_structures.graph import Graph
from addax.kavosh.enumerate import EnumerateSubgraphsFromNode, EnumerateSubgraphsSequentially
from addax.utilities.dataIO import PickleData, ReadPickledData
from addax.visualize.graph import VisualizeGraph



def ClassifySubgraph(G, subgraph, k):
    """
    Given a subgraph, determine the unique certificate corresponding to the canonical labeling of that graph

    @param G: the graph from which the subgraph comes
    @param subgraph: the subgraph extracted via Kavosh for which to determine the canonical label
    @param k: the size of the subgraph

    Returns a certificate for the subgraph corresponding to the canonical labeling of that graph
    """
    # create a mapping from the vertex indices to numpy indices
    vertex_mapping = {}
    for iv, vertex in enumerate(subgraph):
        vertex_mapping[vertex] = iv

    # construct a dictionary of neighbors in the subgraph
    adj_neighbors = {}

    # iterate over all of the indices in the subgraph
    for vertex_index in subgraph:
        # get the index between 0 and k - 1 for this vertex
        index_one = vertex_mapping[vertex_index]

        # keep track of the adjacent neighbors
        adj_neighbors[index_one] = []

        # find all the neighbors to this vertex in the graph
        for neighbor_index in G.vertices[vertex_index].OutgoingNeighborIndices():
            # skip neighbors not in the subgraph
            if not neighbor_index in vertex_mapping: continue

            # get the index between 0 and k - 1 for this vertex
            index_two = vertex_mapping[neighbor_index]

            assert (not index_one == index_two)

            adj_neighbors[index_one].append(index_two)

    # create a pynauty graph of size k for the subgraph with directedness dependent on G
    py_subgraph = pynauty.Graph(k, directed = G.directed, adjacency_dict = adj_neighbors)

    # get the certificate of this subgraph as a byte string
    certificate = pynauty.certificate(py_subgraph)

    return certificate



def ParseCertificate(k, certificate, directed):
    """
    Parse a certificate given by pynauty and determine the canonical motif

    @param k: motif size
    @param certificate: the certificate from Nauty to parse
    @param directed: is the motif a directed or undirected motif

    Returns a canonical graph for this motif as a graph object
    """
    # each vertex gets a certain number of bits in the certificate that contain the adjacency
    # matrix for that vertex. The number of bits is determined by the no_setwords in pynauty. This variable
    # matches the number of 64-bit words (on a 64-bit system) needed to fit all vertices.
    assert (not (len(certificate) % k))
    bytes_per_vertex = len(certificate) // k

    # create an empty graph object
    canonical_graph = Graph(directed)
    for iv in range(k):
        canonical_graph.AddVertex(iv)

    for vertex in range(k):
        certificate_bytes = certificate[vertex * bytes_per_vertex:(vertex + 1) * bytes_per_vertex]

        # the certificate starts with the "farthest right" bit
        for byte_offset, byte in enumerate(certificate_bytes[::-1]):
            assert (byte < 256)

            # 8 bits per byte
            for bit_offset in range(8):
                bit = byte & (2 ** 7)
                byte = byte << 1

                if bit:
                    # determine the neighbor by the byte and bit offsets
                    neighbor_vertex = 8 * byte_offset + (bit_offset)
                    canonical_graph.AddEdge(vertex, neighbor_vertex)

    return canonical_graph



def ClassifySubgraphsFromNode(G, k, u):
    """
    Classify all the subgraphs of size k from vertex u using the Nauty framework

    @param G: graph
    @param k: motif size
    @param u: root vertex index

    Returns a dictionary of certificates
    """

    # create a dictionary for the certificates
    certificates = {}

    # enumerate all subgraphs rooted at vertext u using the Kavosh method
    for subgraph in EnumerateSubgraphsFromNode(G, k, u):
        # get the classification for this subgraph
        certificate = ClassifySubgraph(G, subgraph, k)

        # create a dictionary of certificates
        if not certificate in certificates:
            certificates[certificate] = 1
        else:
            certificates[certificate] += 1

    # save the certificates to disk
    output_directory = 'temp/{}'.format(G.prefix)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok = True)
    output_filename = '{}/motif-size-{:03d}-node-{:016d}-certificates.pickle'.format(output_directory, k, u)
    PickleData(certificates, output_filename)

    return certificates



def ClassifySubgraphsSequentially(G, k):
    """
    Classify all subgraphs of size k using the Nauty framework

    @param G: graph
    @param k: motif size

    Returns a dictionary of certificates
    """

    # create a dictionary for the certificates
    certificates = {}

    # enumerate all subgraphs using the Kavosh method
    for subgraph in EnumerateSubgraphsSequentially(G, k):
        # get the classification for this subgraph
        certificate = ClassifySubgraph(G, subgraph, k)

        # create a dictionary of certificates
        if not certificate in certificates:
            certificates[certificate] = 1
        else:
            certificates[certificate] += 1

    # save the certificates to disk
    output_directory = 'temp/{}'.format(G.prefix)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok = True)
    output_filename = '{}/motif-size-{:03d}-certificates.pickle'.format(output_directory, k)
    PickleData(certificates, output_filename)

    return certificates



def CombineSubgraphCertificates(G, k):
    """
    Combine the generated certificates for all nodes in the graph G

    @param G: graph
    @param k: motif size

    Returns a dictionary of certificates
    """

    certificates = {}

    # get the input/output directory
    directory = 'temp/{}'.format(G.prefix)
    if not os.path.exists(directory):
        sys.stderr.write('Failed to find any input certificates...\n')
        sys.exit(-1)

    for u in G.vertices.keys():
        # make sure that the input file exists
        input_filename = '{}/motif-size-{:03d}-node-{:016d}-certificates.pickle'.format(directory, k, u)
        if not os.path.exists(input_filename):
            sys.stderr.write('Missing certificates file: {}...\n'.format(input_filename))
            sys.exit(-1)

        # read the input certifiicates from disk
        certificates_from_node = ReadPickledData(input_filename)

        # add the certificates to the global array
        for certificate in certificates_from_node:
            if not certificate in certificates:
                certificates[certificate] = certificates_from_node[certificate]
            else:
                certificates[certificate] += certificates_from_node[certificate]

    # save the certificates to disk
    output_directory = 'temp/{}'.format(G.prefix)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok = True)
    output_filename = '{}/motif-size-{:03d}-certificates.pickle'.format(output_directory, k)
    PickleData(certificates, output_filename)

    return certificates
