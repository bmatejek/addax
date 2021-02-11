import glob



import numpy as np



from addax.kavosh.enumerate import CreateDirectoryStructure
from addax.utilities.dataIO import ReadGraph



def PrintRunningTimeStatistics(input_filename, motif_size, vertex_colored = False, edge_colored = False, community_based = False):
    """
    Print statistics for the running times for each vertex in the graph

    @param input_filename: location for the graph that was enumerated
    @parak motif_size: the motif subgraph size to find
    @param vertex_colored: a boolean flag to allow for vertex colors
    @param edge_colored: a boolean flag to allow for edge colors
    @param community_based: a boolean flag to only enumerate subgraphs in the same community
    """
    # read the input graph
    graph = ReadGraph(input_filename, vertices_only = True)

    # get the temp directory
    temp_directory = CreateDirectoryStructure(input_filename, vertex_colored, edge_colored, community_based, False)

    # get a list of all the input filenames for this motif size
    filenames = sorted(glob.glob('{}/certificates/motif-size-{:03d}-*.txt'.format(temp_directory, motif_size)))

    # create a dictionary of running times
    running_times = {}

    # read all of the filenames corresponding to enumeration data
    for filename in filenames:
        # open the file
        with open(filename, 'r') as fd:
            # only worry about the lines with Enumerated summaries
            for line in fd.readlines():
                if not line.startswith('Enumerated'): continue

                # get the relevant stats from this line
                segments = line.split()
                _, vertex_index, time = int(segments[1]), int(segments[5]), float(segments[7])

                # add to the mapping of running times
                running_times[vertex_index] = time

    # make sure that every vertex is accounted for
    assert (sorted(list(running_times.keys())) == sorted(list(graph.vertices.keys())))

    # print statistics on the running time
    running_times = list(running_times.values())

    print ('Running Time Statistics for {} Motif Size {}'.format(graph.prefix, motif_size))
    print ('  Mean Time: {:0.2f} seconds'.format(np.mean(running_times)))
    print ('  Median Time: {:0.2f} seconds'.format(np.median(running_times)))
    print ('  Maximum Time: {:0.2f} seconds'.format(np.max(running_times)))
    print ('  Std Dev: {:0.2f} seconds'.format(np.std(running_times)))
    print ('  Total Time: {:0.2f} seconds'.format(sum(running_times)))
