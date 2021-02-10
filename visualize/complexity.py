import glob



import matplotlib
import matplotlib.pyplot as plt



import numpy as np


plt.style.use('seaborn')



def GetRunningTimes(graph, motif_sizes):
    """
    Get the running times for each nodes for this graph

    @param graph: the input graph to analyze
    @param motif_size: an array of motif sizes to consider
    """
    # get the input directory (hard-coded for this one type)
    input_directory = 'temp/{}/global/vertex-agnostic/edge-agnostic/certificates'.format(graph.prefix)

    # create an empty dictionary of runnings times
    running_times = {}

    # each vertex gets its own dictionary
    for vertex_index in graph.vertices.keys():
        running_times[vertex_index] = {}

    # get all the filenames for this configuration
    for motif_size in motif_sizes:
        filenames = glob.glob('{}/motif-size-{:03d}-*'.format(input_directory, motif_size))

        # iterate over all filenames
        for filename in filenames:
            # read the file
            with open(filename, 'r') as fd:
                for line in fd.readlines():
                    if not line.startswith('Enumerated'): continue

                    # break the line into segments to get the relevant information
                    segments = line.split()
                    nsubgraphs, vertex_index, time = int(segments[1]), int(segments[5]), float(segments[7])

                    running_times[vertex_index]['Motif Size {} Subgraphs'.format(motif_size)] = nsubgraphs
                    running_times[vertex_index]['Motif Size {} Time'.format(motif_size)] = time

        assert (sorted(list(running_times.keys())) == sorted(list(graph.vertices.keys())))

    return running_times



def GetSecondOrderStatistics(graph, running_times):
    """
    Get the second order statistics for the nodes for this graph.

    @param graph: the input graph to analyze
    @param running_times: a dictionary of running times for each vertex
    """
    # get the number of neighbors with a higher value than this vertex
    for vertex_index, vertex in graph.vertices.items():
        running_times[vertex_index]['Immediate Neighborhood'] = vertex.neighbors

    # get the number of neighbors of neighbors
    for vertex_index, vertex in graph.vertices.items():
        # create a new set with the existing neighbors
        running_times[vertex_index]['Extended Neighborhood'] = set(running_times[vertex_index]['Immediate Neighborhood'])
        # add the neighbors of the neighbors
        for neighbor_index in running_times[vertex_index]['Immediate Neighborhood']:
            # use update to union the two sets in place
            running_times[vertex_index]['Extended Neighborhood'].update(running_times[neighbor_index]['Immediate Neighborhood'])

    # condense the neighborhoods down to a single number
    for vertex_index in running_times.keys():
        # remove the vertices that cannot belong to this subgraph tree
        running_times[vertex_index]['Immediate Neighborhood'] = [neighbor_index for neighbor_index in running_times[vertex_index]['Immediate Neighborhood'] if vertex_index < neighbor_index ]
        running_times[vertex_index]['Immediate Neighborhood'] = len(running_times[vertex_index]['Immediate Neighborhood'])
        # remove the vertices that cannot belong to this subgraph tree
        running_times[vertex_index]['Extended Neighborhood'] = [neighbor_index for neighbor_index in running_times[vertex_index]['Extended Neighborhood'] if vertex_index < neighbor_index ]
        running_times[vertex_index]['Extended Neighborhood'] = len(running_times[vertex_index]['Extended Neighborhood'])



def PlotComplexities(graph, motif_size, running_times, attribute_one, attribute_two):
    """
    Plot the statistics about the computational complexity

    @param graph: the input graph to analyze
    @param motif_size: the size of motifs to consider
    @param running_times: a dictionary of running times for each vertex
    @param attribute: the attribute to compare against
    """
    # create a figure for this statistic
    fig = plt.figure(figsize=(10, 5))

    # create the x and y variables
    x = []
    y = []

    # iterate over all variables
    for vertex, attributes in running_times.items():
        x.append(attributes[attribute_one])
        y.append(attributes[attribute_two])


    # create an axis for the figure
    ax = fig.add_subplot(111)

    # plot the attribute
    ax.scatter(x, y, label=attribute_one)

    # set the titles
    ax.set_title('{} as a Function of {} Motifs Size: {}'.format(attribute_two, attribute_one, motif_size), fontsize=20)
    ax.set_xlabel(attribute_one, fontsize=18)
    ax.set_ylabel(attribute_two, fontsize=18)

    fig.tight_layout()

    plt.show()

    plt.clf()



def ComputationalComplexityPerNode(graph, motif_sizes):
    """
    Calculate the complexity of enumeration per node.

    @param graph: the input graph to analyze
    @param motif_sizes: an array of motif sizes to consider
    """
    running_times = GetRunningTimes(graph, motif_sizes)

    # create second order statistsics for the vertices
    GetSecondOrderStatistics(graph, running_times)

    for motif_size in motif_sizes:
        PlotComplexities(graph, motif_size, running_times, 'Motif Size {} Subgraphs'.format(motif_size), 'Motif Size {} Time'.format(motif_size))
        PlotComplexities(graph, motif_size, running_times, 'Immediate Neighborhood', 'Motif Size {} Time'.format(motif_size))
        PlotComplexities(graph, motif_size, running_times, 'Extended Neighborhood', 'Motif Size {} Time'.format(motif_size))
        PlotComplexities(graph, motif_size, running_times, 'Immediate Neighborhood', 'Motif Size {} Subgraphs'.format(motif_size))
        PlotComplexities(graph, motif_size, running_times, 'Extended Neighborhood', 'Motif Size {} Subgraphs'.format(motif_size))
        PlotComplexities(graph, motif_size, running_times, 'Immediate Neighborhood', 'Extended Neighborhood')

    for iv1 in range(len(motif_sizes)):
        for iv2 in range(iv1 + 1, len(motif_sizes)):
            PlotComplexities(graph, motif_size, running_times, 'Motif Size {} Time'.format(motif_sizes[iv1]), 'Motif Size {} Time'.format(motif_sizes[iv2]))
            PlotComplexities(graph, motif_size, running_times, 'Motif Size {} Subgraphs'.format(motif_sizes[iv1]), 'Motif Size {} Subgraphs'.format(motif_sizes[iv2]))
