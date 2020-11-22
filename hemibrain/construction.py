import os
import csv
import time



from addax.data_structures.graph import Graph
from addax.utilities.dataIO import WriteGraph



def ConstructGraphFromHemiBrainCSV():
    """
    Construct a graph object for the hemi brain CSV files
    """
    # start statistics
    start_time = time.time()

    neurons = {}
    edges = {}
    synapses_per_region = {}

    # open the neuron csv file
    neuron_filename = 'CSVs/HemiBrain/traced-neurons.csv'
    with open(neuron_filename, newline='') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')

        # skip header
        next(csv_reader, None)
        for row in csv_reader:
            neuron_id = int(row[0])

            assert (not neuron_id in neurons)

            # create a region mapping for each neuron
            neurons[neuron_id] = {}

    # read the synapses by region of interest
    synapse_filename = 'CSVs/HemiBrain/traced-roi-connections.csv'
    with open(synapse_filename, newline='') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')

        # skip header
        next(csv_reader, None)
        for row in csv_reader:
            pre_neuron_id = int(row[0])
            post_neuron_id = int(row[1])
            roi = str(row[2])
            weight = int(row[3])

            assert (pre_neuron_id in neurons)
            assert (post_neuron_id in neurons)

            # increase the incidence for this region for the pre synaptic neuron
            if not roi in neurons[pre_neuron_id]:
                neurons[pre_neuron_id][roi] = weight
            else:
                neurons[pre_neuron_id][roi] += weight

            # increase the incidence for this region for the post synaptic neuron
            if not roi in neurons[post_neuron_id]:
                neurons[post_neuron_id][roi] = weight
            else:
                neurons[post_neuron_id][roi] += weight

            # update the connective strength between these two neurons
            if not (pre_neuron_id, post_neuron_id) in edges:
                edges[(pre_neuron_id, post_neuron_id)] = weight
            else:
                edges[(pre_neuron_id, post_neuron_id)] += weight

            # update region statistics
            if not roi in synapses_per_region:
                synapses_per_region[roi] = weight
            else:
                synapses_per_region[roi] += weight

    # determine the region (community/cluster) for each neuron
    for neuron in neurons.keys():
        max_region = None
        max_synaptic_weight = -1

        for region in neurons[neuron].keys():
            if neurons[neuron][region] > max_synaptic_weight:
                max_synaptic_weight = neurons[neuron][region]
                max_region = region

        neurons[neuron] = max_region

    regions_to_communities = {}
    for index, (region, _) in enumerate(sorted(synapses_per_region.items(), key = lambda x : x[1], reverse = True)):
        regions_to_communities[region] = index

    # construct a directed graph object
    graph = Graph('hemi-brain', directed = True, colored = False)

    # add vertices with the communities
    for enumeration_index, neuron in enumerate(sorted(neurons.keys())):
        # vertices start with  the default coloring of -1, enumeration_index of -1
        graph.AddVertex(neuron, enumeration_index, community = regions_to_communities[neurons[neuron]], color = -1)

    # add edges with the synaptic weights
    for (pre_neuron_id, post_neuron_id) in edges.keys():
        graph.AddEdge(pre_neuron_id, post_neuron_id, weight = edges[(pre_neuron_id, post_neuron_id)])

    # write the hemibrain graph to disk
    if not os.path.exists('graphs'):
        os.makedirs('graphs', exist_ok = True)
    output_filename = 'graphs/hemi-brain.graph.bz2'
    WriteGraph(graph, output_filename)

    # write the condensed CSV neuron and edge files
    with open('CSVs/HemiBrain/condensed-neurons.csv', 'w') as fd:
        fd.write('Neuron ID,Community\n')
        for neuron in graph.vertices.values():
            fd.write('{},{}\n'.format(neuron.index, neuron.community))

    with open('CSVs/HemiBrain/condensed-edges.csv', 'w') as fd:
        fd.write('Pre Synaptic Neuron ID,Post Synaptic Neuron ID,Weight\n')
        for edge in graph.edges:
            fd.write('{},{},{}\n'.format(edge.source_index, edge.destination_index, edge.weight))

    # print statistics
    print ('Wrote {} in {:0.2f} seconds.'.format(output_filename, time.time() - start_time))
