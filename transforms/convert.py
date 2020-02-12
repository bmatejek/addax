import csv


def HemibrainCSV2Graph():
    """
    Converts the hemibrain csv data provided by Janelia into a custom built
    graph format for better visualization and analysis.
    """

    neurons = set()
    edges = {}

    with open('csvs/hemibrain/traced-total-connections.csv') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        for row in csv_reader:
            # read the neurons and synaptic weight
            # the edge are directed from neuron one to neuron two
            neuron_one = row[0]
            neuron_two = row[1]

            synaptic_weight = row[2]

            # add the neurons to the set
            neurons.add(neuron_one)
            neurons.add(neuron_two)

            # add in the edge weights
            edges[(neuron_one, neuron_two)] = synaptic_weight

    print ('No. Neurons: {}'.format(len(neurons)))
    print ('No. Synapses: {}'.format(len(edges)))
