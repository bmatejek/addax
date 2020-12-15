import networkx as nx



from addax.utilities.dataIO import ReadGraph



def ParseCertificate(graph, k, certificate, color_mapping = {}):
    """
    Produce a canonical graph from a certificate

    @param graph: the input graph from which the certificate was generated
    @param k: motif size
    @param certificate: the certificate from Nauty to parse
    @param color_mapping: a dictionary that converts color indices to strings
    """
    # this will require substantial debugging and validation for k > 256
    # currently skipped since motifs that size are computationally infeasible
    # nauty almost certainly can never run on graphs that size
    assert (k < 256)

    # each vertex gets a certain number of bits in the certificate that contain the adjacency
    # matrix for that vertex. The number of bits is determined by the no_setwords in pynauty. This variable
    # matches the number of 64-bit words (on a 64-bit system) needed to fit all vertices.
    # for colored graphs we also add on 64 bits for each node in the subgraph for coloring
    if graph.colored:
        # if the graph is colored, there are 16 characters per vertex
        # rest of the string represents the canonical labeling
        coloring = certificate[-16 * k:]
        certificate = certificate[:-16 * k]

        # convert the vertex bytes (in hex) to base 10 integer
        vertex_colors = [int(coloring[16 * iv: 16 * (iv + 1)], 16) for iv in range(k)]

    # create a new networkx graph object
    if graph.directed:
        nx_graph = nx.DiGraph()
    else:
        nx_graph = nx.Graph()

    # add a vertex for each of the k nodes in the subgraph
    for index in range(k):
        # colored graphs have labels for their colors
        if graph.colored:
            if vertex_colors[index] in color_mapping: nx_graph.add_node(index, label = color_mapping[vertex_colors[index]])
            else:nx_graph.add_node(index, label = 'Color: {}'.format(vertex_colors[index]))
        else:
            nx_graph.add_node(index)

    # get the number of words per veretx (must be a multiple of 2 * k)
    assert (not (len(certificate) % (2 * k)))
    # based on our output of using hexademical in cpp-enumerate.cpp (i.e., %02x),
    # each byte is written with two letters in our input string (e.g., aa, 02, 10, etc.)
    bytes_per_vertex = len(certificate) // k // 2

    # iterate over every vertex and extract the corresponding adjacency matrix
    for vertex in range(k):
        # multiple by two since we are using two characters in the string per byte (wrriten in hexademical)
        certificate_bytes = certificate[vertex * bytes_per_vertex * 2:(vertex + 1) * bytes_per_vertex * 2]

        for byte_offset, iv in enumerate(range(0, len(certificate_bytes), 2)):
            # get the byte as an integer (using hexadecimal currently)
            byte = int(certificate_bytes[iv:iv+2], 16)
            assert (byte < 256)

            # 8 bits per byte (little endian)
            for bit_offset in range(8):
                bit = byte & (2 ** 7)
                byte = byte << 1

                # if this bit is 1, there is an edge from vertex to this location
                if bit:
                    # determine the neighbor by the byte and bit offset
                    neighbor_vertex = 8 * (bytes_per_vertex - byte_offset - 1) + bit_offset
                    nx_graph.add_edge(vertex, neighbor_vertex)

    return nx_graph



def ParseCertificates(input_filename, k, community_based = False, color_mapping = {}):
    """
    Parse the certificates generated for this subgraph

    @param input_filename: location for the graph to enumerate
    @parak k: the motif subgraph size to find
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    @param color_mapping: a dictionary that converts color indices to strings
    """
    # read the graph
    graph = ReadGraph(input_filename, vertices_only = True)

    # read the combined enumerated subgraphs file
    if community_based: subgraphs_filename = 'subgraphs-community-based/{}/motif-size-{:03d}-certificates.txt'.format(graph.prefix, k)
    else: subgraphs_filename = 'subgraphs/{}/motif-size-{:03d}-certificates.txt'.format(graph.prefix, k)

    with open(subgraphs_filename, 'r') as fd:
        # read all of the certificates and enumerated subgraphs
        for subgraph_index, certificate_info in enumerate(fd.readlines()[:-1]):
            certificate = certificate_info.split(':')[0].strip()
            nsubgraphs = int(certificate_info.split(':')[1].strip())

            # parse the certificate for this graph
            nx_graph = ParseCertificate(graph, k, certificate, color_mapping)

            # create the graph drawing structure
            A = nx.nx_agraph.to_agraph(nx_graph)
            A.layout(prog='dot')

            if community_based: output_filename = 'subgraphs-community-based/{}/motif-size-{:03d}-motif-{}-found-{}.dot'.format(graph.prefix, k, subgraph_index, nsubgraphs)
            else: output_filename = 'subgraphs/{}/motif-size-{:03d}-motif-{}-found-{}.dot'.format(graph.prefix, k, subgraph_index, nsubgraphs)

            A.draw(output_filename)
