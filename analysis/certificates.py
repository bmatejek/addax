from addax.utilities.dataIO import ReadGraph
from addax.kavosh.enumerate import CreateDirectoryStructure



def ReadCertificates(input_filename, k, vertex_colored, edge_colored, community_based):
    """
    Read the certificates for this graph from the subgraphs directory

    @param input_filename: location for the graph to enumerate
    @parak k: the motif subgraph size to find
    @param vertex_colored: a boolean flag to allow for vertex colors
    @param edge_colored: a boolean flag to allow for edge colors
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    """
    # read the graph
    graph = ReadGraph(input_filename, vertices_only = True)

    # get the temp directory
    temp_directory = CreateDirectoryStructure(input_filename, vertex_colored, edge_colored, community_based, False)

    # get the input directory
    input_directory = 'subgraphs/{}'.format('/'.join(temp_directory.split('/')[1:]))

    # read the combined enumerated subgraphs file
    subgraphs_filename = '{}/motif-size-{:03d}-certificates.txt'.format(input_directory, k)

    certificates = {}

    # keep track of the total subgraphs
    total_subgraphs = 0

    with open(subgraphs_filename, 'r') as fd:
        # read all of the certificates and enumerated subgraphs
        for certificate_info in fd:
            if certificate_info.startswith('Found'): continue
            elif certificate_info.startswith('Enumerated'):
                total_subgraphs = int(certificate_info.split()[1])
                total_time = float(certificate_info.split()[-2])
            else:
                certificate = certificate_info.split(':')[0].strip()
                nsubgraphs = int(certificate_info.split(':')[1].strip())

                certificates[certificate] = nsubgraphs

    # return the certificates and the number of subgraphs
    return certificates, total_subgraphs, total_time
