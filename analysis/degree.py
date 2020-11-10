import math



import matplotlib
import matplotlib.pyplot as plt



def PlotDegrees(degrees):
    """
    Plot the degree distribution.

    @param degrees: a list of degrees from the graph
    """
    # get the number of vertices in the grpah
    nvertices = len(degrees)
    # get the maximum degree
    max_degree = max(degrees)

    # construct a degree distribution
    degree_distribution = [0 for _ in range(max_degree + 1)]
    for degree in degrees:
        degree_distribution[degree] += 1

    # get the probability of each degree
    cdf = [1.0 - sum(degree_distribution[iv:]) / nvertices for iv in range(1, max_degree + 1)]

    plt.title('CDF of Degree Distribution')
    plt.xlabel('Vertex Degree')
    plt.ylabel('Probability')

    plt.plot(cdf)
    plt.show()




def ModelDegreeDistribution(graph):
    """
    Model the distribution of edge degrees for each vertex

    @param graph: the graph for which to analyze the edge distribution
    """
    # get a list of all degrees
    incoming_degrees = []
    outgoing_degrees = []
    degrees = []

    # iterate over all vertices in the graph
    for index, vertex in graph.vertices.items():
        incoming_degrees.append(len(vertex.incoming_neighbors))
        outgoing_degrees.append(len(vertex.outgoing_neighbors))
        degrees.append(len(vertex.neighbors))

    PlotDegrees(degrees)
    PlotDegrees(incoming_degrees)
    PlotDegrees(outgoing_degrees)
