import os
import glob



import numpy as np



cimport cython
cimport numpy as np
import ctypes
from libcpp cimport bool



from addax.utilities.dataIO import ReadGraph, ReadPrefix



cdef extern from 'cpp-enumerate.h':
    void CppSetCommunityBased(bool community_based)
    void CppSetWriteSubgraphs(bool write_subgraphs)
    void CppEnumerateSubgraphsSequentially(const char *prefix, short k)
    void CppEnumerateSubgraphsFromNodes(const char *prefix, short k, long *nodes, long nnodes, long output_suffix)



def CreateDirectoryStructure(input_filename, community_based, write_subgraphs):
    """
    Create the directory structure for enumeration

    @param input_filename: location for the graph to enumerate
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    @param write_subgraphs: a boolean flag to write all enumerated subgraphs to disk
    """
    # is this directory community based
    if community_based: tmp_directory = 'temp-community-based'
    else: tmp_directory = 'temp'

    # get the prefix for the dataset
    prefix = ReadPrefix(input_filename)

    # create the certificate and subgraph directory
    certificate_directory = '{}/{}/certificates'.format(tmp_directory, prefix)
    if not os.path.exists(certificate_directory):
        os.makedirs(certificate_directory, exist_ok = True)

    if write_subgraphs:
        subgraphs_directory = '{}/{}/subgraphs'.format(tmp_directory, prefix)
        if not os.path.exists(subgraphs_directory):
            os.makedirs(subgraphs_directory, exist_ok = True)



def EnumerateSubgraphsSequentially(input_filename, k, community_based = False, write_subgraphs = False):
    """
    Enumerate all subgraphs in the graph specified by input_filename

    @param input_filename: location for the graph to enumerate
    @parak k: the motif subgraph size to find
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    @param write_subgraphs: a boolean flag to write all enumerated subgraphs to disk
    """
    import time
    start_time = time.time()
    # create the temp directory if it does not exist
    CreateDirectoryStructure(input_filename, community_based, write_subgraphs)

    # set the community based flag
    CppSetCommunityBased(community_based)
    # set the write subgraphs flag
    CppSetWriteSubgraphs(write_subgraphs)

    # enumerate the subgraph, cast the string into a character array
    CppEnumerateSubgraphsSequentially(input_filename.encode('utf-8'), k)
    print ('Enumerated all subgraphs in {:0.2f} seconds'.format(time.time() - start_time))



def EnumerateSubgraphsFromNodes(input_filename, k, nodes, output_suffix, community_based = False, write_subgraphs = False):
    """
    Enumerate all subgraphs in the graph starting at the nodes array

    @param input_filename: location for the graph to enumerate
    @parak k: the motif subgraph size to find
    @param nodes: an array of nodes to enumerate starting at
    @param output_suffix: a integer identifying a unique file to which to save the results
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    @param write_subgraphs: a boolean flag to write all enumerated subgraphs to disk
    """
    # create the temp directory if it does not exist
    CreateDirectoryStructure(input_filename, community_based, write_subgraphs)

    # set the community based flag
    CppSetCommunityBased(community_based)
    # set the write subgraphs flag
    CppSetWriteSubgraphs(write_subgraphs)

    # convert the array of nodes into a c array
    nnodes = len(nodes)
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_nodes = np.ascontiguousarray(nodes, dtype=ctypes.c_int64)

    # enumerate the subgraph, cast the string into a character array
    CppEnumerateSubgraphsFromNodes(input_filename.encode('utf-8'), k, &(cpp_nodes[0]), nnodes, output_suffix)

    # free memory
    del cpp_nodes



def ParseCertificateFile(input_filename):
    """
    Read the input file and validate that it is correct

    @param input_filename: location of the certificates file saved in enumeration
    """
    # there are two modes, certificate mode and summary mode
    certificate_mode = True

    certificates = {}
    total_nsubgraphs, total_time = 0, 0







def CombineEnumeratedSubgraphs(input_filename, k, community_based = False):
    """
    Combine all of the enumerated subgraphs for a given file and motif size.

    @param input_filename: location for the graph to enumerate
    @parak k: the motif subgraph size to find
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    """
    # read the graph (only vertices)
    graph = ReadGraph(input_filename, vertices_only = True)

    # create the list of vertices
    vertices = set(list(graph.vertices.keys()))

    # get the temporary directory
    if community_based: tmp_directory = 'temp-community-based/{}/certificates'.format(graph.prefix)
    else: tmp_directory = 'temp/{}/certificates'.format(graph.prefix)

    # get a list of all the input filenames for this motif size
    input_filenames = sorted(glob.glob('{}/motif-size-{:03d}-*.txt'.format(tmp_directory, k)))

    # create a dictionary of certificates
    certificates = {}

    # set initial counter variables
    total_nsubgraphs, total_time = 0, 0

    # iterate over all the input files
    for input_filename in input_filenames:
        # open the output file
        with open(input_filename, 'r') as fd:
            for certificate_line in fd:
                segments = certificate_line.split()

                # update the mode that currently exists
                if segments[0] == 'Enumerated':
                    nsubgraphs, vertex, time = int(segments[1]), int(segments[5]), float(segments[7])

                    # update the counter variables that verify correctness
                    total_nsubgraphs += nsubgraphs
                    total_time += time
                    vertices.remove(vertex)

                    certificate_mode = False
                else:
                    certificate, nsubgraphs = segments[0].strip(':'), int(segments[1])

                    # update the certificate information
                    if not certificate in certificates: certificates[certificate] = nsubgraphs
                    else: certificates[certificate] += nsubgraphs

                    certificate_mode = True

        assert (not certificate_mode)

    # all vertices are found
    assert (not len(vertices))
    # the number of output subgraphs is correctly tabulated
    assert (sum(certificates.values()) == total_nsubgraphs)

    # create the output directory if it does not exist
    if community_based: output_directory = 'subgraphs-community-based/{}'.format(graph.prefix)
    else: output_directory = 'subgraphs/{}'.format(graph.prefix)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok = True)

    output_filename = '{}/motif-size-{:03d}-certificates.txt'.format(output_directory, k)
    with open(output_filename, 'w') as fd:
        for certificate, nsubgraphs in certificates.items():
            fd.write('{}: {}\n'.format(certificate, nsubgraphs))
            print ('{}: {}'.format(certificate, nsubgraphs))
        fd.write('Enumerated {} subgraphs in {:0.2f} seconds.'.format(total_nsubgraphs, total_time))
        print ('Enumerated {} subgraphs in {:0.2f} seconds.'.format(total_nsubgraphs, total_time))
