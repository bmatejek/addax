import os



import numpy as np



cimport cython
cimport numpy as np
import ctypes
from libcpp cimport bool



from addax.utilities.dataIO import ReadGraph, ReadPrefix



cdef extern from 'cpp-enumerate.h':
    void CppSetCommunityBased(bool community_based)
    void CppEnumerateSubgraphsSequentially(const char *prefix, short k)
    void CppEnumerateSubgraphsFromNode(const char *prefix, short k, long u)
    void CppEnumerateSubgraphsFromNodes(const char *prefix, short k, long *nodes, long nnodes)



def EnumerateSubgraphsSequentially(input_filename, k, community_based = False):
    """
    Enumerate all subgraphs in the graph specified by input_filename

    @param input_filename: location for the graph to enumerate
    @parak k: the motif subgraph size to find
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    """
    # create the temp directory if it does not exist
    if community_based: tmp_directory = 'temp/{}-community-based'.format(ReadPrefix(input_filename))
    else: tmp_directory = 'temp/{}'.format(ReadPrefix(input_filename))
    if not os.path.exists(tmp_directory):
        os.makedirs(tmp_directory, exist_ok = True)

    if community_based: timing_directory = 'temp/{}-community-based/timing'.format(ReadPrefix(input_filename))
    else: timing_directory = 'temp/{}/timing'.format(ReadPrefix(input_filename))
    if not os.path.exists(timing_directory):
        os.makedirs(timing_directory, exist_ok = True)

    # set the community basd flag
    CppSetCommunityBased(community_based)

    # enumerate the subgraph, cast the string into a character array
    CppEnumerateSubgraphsSequentially(input_filename.encode('utf-8'), k)



def EnumerateSubgraphsFromNode(input_filename, k, u, community_based = False):
    """
    Enumerate all subgraphs in the graph starting at the nodes array

    @param input_filename: location for the graph to enumerate
    @parak k: the motif subgraph size to find
    @param u: the vertex from which to enumerate
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    """
    # create the temp directory if it does not exist
    if community_based: tmp_directory = 'temp/{}-community-based'.format(ReadPrefix(input_filename))
    else: tmp_directory = 'temp/{}'.format(ReadPrefix(input_filename))
    if not os.path.exists(tmp_directory):
        os.makedirs(tmp_directory, exist_ok = True)

    if community_based: timing_directory = 'temp/{}-community-based/timing'.format(ReadPrefix(input_filename))
    else: timing_directory = 'temp/{}/timing'.format(ReadPrefix(input_filename))
    if not os.path.exists(timing_directory):
        os.makedirs(timing_directory, exist_ok = True)

    # set the community basd flag
    CppSetCommunityBased(community_based)

    # enumerate the subgraph, cast the string into a character array
    CppEnumerateSubgraphsFromNode(input_filename.encode('utf-8'), k, u)



def EnumerateSubgraphsFromNodes(input_filename, k, nodes, community_based = False):
    """
    Enumerate all subgraphs in the graph starting at the nodes array

    @param input_filename: location for the graph to enumerate
    @parak k: the motif subgraph size to find
    @param nodes: an array of nodes to enumerate starting at
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    """
    # create the temp directory if it does not exist
    if community_based: tmp_directory = 'temp/{}-community-based'.format(ReadPrefix(input_filename))
    else: tmp_directory = 'temp/{}'.format(ReadPrefix(input_filename))
    if not os.path.exists(tmp_directory):
        os.makedirs(tmp_directory, exist_ok = True)

    if community_based: timing_directory = 'temp/{}-community-based/timing'.format(ReadPrefix(input_filename))
    else: timing_directory = 'temp/{}/timing'.format(ReadPrefix(input_filename))
    if not os.path.exists(timing_directory):
        os.makedirs(timing_directory, exist_ok = True)

    # set the community basd flag
    CppSetCommunityBased(community_based)

    # convert the array of nodes into a c array
    nnodes = len(nodes)
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_nodes = np.ascontiguousarray(nodes, dtype=ctypes.c_int64)

    # enumerate the subgraph, cast the string into a character array
    CppEnumerateSubgraphsFromNodes(input_filename.encode('utf-8'), k, &(cpp_nodes[0]), nnodes)

    # free memory
    del cpp_nodes



def CombineEnumeratedSubgraphs(input_filename, k, community_based = False):
    """
    Enumerate all subgraphs in the graph specified by input_filename

    @param input_filename: location for the graph to enumerate
    @parak k: the motif subgraph size to find
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    """
    # read the graph
    graph = ReadGraph(input_filename, vertices_only = True)

    # create a global dictionary for each certificate
    certificates = {}

    # iterate over all vertices
    for vertex in graph.vertices.keys():
        # open the filename with the enumerated subgraphs
        if community_based: subgraph_filename = 'temp/{}-community-based/motif-size-{:03d}-node-{:016d}-certificates.txt'.format(graph.prefix, k, vertex)
        else: subgraph_filename = 'temp/{}/motif-size-{:03d}-node-{:016d}-certificates.txt'.format(graph.prefix, k, vertex)

        with open(subgraph_filename, 'r') as fd:
            for certificate_info in fd.readlines():
                # get the certificate and the number of subgraphs
                certificate = certificate_info.split(':')[0].strip()
                nsubgraphs = int(certificate_info.split(':')[1].strip())

                # update the number of subgraphs for this certificate
                if not certificate in certificates:
                    certificates[certificate] = nsubgraphs
                else:
                    certificates[certificate] += nsubgraphs

    # create the output directory if it does not exist
    if community_based: output_directory = 'subgraphs/{}-community-based'.format(graph.prefix)
    else: output_directory = 'subgraphs/{}'.format(graph.prefix)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok = True)

    # write the combined subgraphs to an output file
    output_filename = '{}/motif-size-{:03d}.txt'.format(output_directory, k)
    with open(output_filename, 'w') as fd:
        # iterate over all certificates
        for certificate, nsubgraphs in sorted(certificates.items(), key = lambda x : x[1], reverse = True):
            fd.write('{}: {}\n'.format(certificate, nsubgraphs))
