import time
import xlrd



from addax.data_structures.graph import Graph
from addax.data_structures.enumeration import CalculateAscendingEnumerationIndex
from addax.utilities.dataIO import ReadGraph, WriteGraph



def ConstructGraphsFromTimelapsedXLSX():
    """
    Construct a set of graphs from the timelapsed xlsx file.
    """
    input_filename = 'CSVs/Timelapsed/connectomes.xlsx'
    workbook = xlrd.open_workbook(input_filename)

    # iterate over all sheets
    nsheets = workbook.nsheets
    for iv in range(nsheets):
        # start statistics
        start_time = time.time()

        sheet = workbook.sheet_by_index(iv)

        # get the size of the adjacency matrix
        nrows = sheet.nrows
        ncols = sheet.ncols

        # create an empty graph
        graph = Graph('C-elegans-timelapsed-{:02}'.format(iv + 1), directed = True, vertex_colored = True, edge_colored = False)

        # craete a mapping from vertex names to indices
        vertex_index_to_name = {}
        for vertex_index, row_index in enumerate(range(2, nrows)):
            # all neuron names are in column B
            neuron_name = sheet.cell_value(row_index, 1)
            vertex_index_to_name[vertex_index] = neuron_name

        # add vertices with communities at -1
        for enumeration_index, (vertex_index, neuron_type) in enumerate(vertex_index_to_name.items()):
            # each vertex receives a unique color since the graphs are stereotyped
            graph.AddVertex(vertex_index, enumeration_index, community = -1, color = vertex_index)

        # add edges based on the wiring diagram
        for pre_synaptic_vertex_index, row_index in enumerate(range(2, nrows)):
            for post_synaptic_vertex_index, col_index in enumerate(range(2, ncols)):
                synaptic_strength = sheet.cell_value(row_index, col_index)
                # skip over edges that do not exist
                if not synaptic_strength: continue

                # skip over self loops
                if pre_synaptic_vertex_index == post_synaptic_vertex_index: continue 

                graph.AddEdge(pre_synaptic_vertex_index, post_synaptic_vertex_index, weight=synaptic_strength)

        graph.SetVertexTypeMapping(vertex_index_to_name)

        # get the output filename
        output_filename = 'graphs/{}.graph.bz2'.format(graph.prefix)
        WriteGraph(graph, output_filename)

        # write the condensed CSV neuron and edge files
        with open('CSVs/Timelapsed/C-elegans-timelapsed-{:02d}-condensed-neurons.csv'.format(iv + 1), 'w') as fd:
            fd.write('Neuron ID,Color Name\n')
            for neuron in graph.vertices.values():
                fd.write('{},{}\n'.format(neuron.index, graph.vertex_type_mapping[neuron.color]))

        with open('CSVs/Timelapsed/C-elegans-timelapsed-{:02d}-condensed-edges.csv'.format(iv + 1), 'w') as fd:
            fd.write('Pre Synaptic Neuron ID,Post Synaptic Neuron ID,Weight\n')
            for edge in graph.edges.values():
                fd.write('{},{},{}\n'.format(edge.source_index, edge.destination_index, edge.weight))


        # print statistics
        print ('No. Neurons: {}'.format(len(graph.vertices.keys())))
        print ('No. Edges: {}'.format(len(graph.edges)))

        print ('Wrote {} in {:0.2f} seconds.'.format(output_filename, time.time() - start_time))

        CalculateAscendingEnumerationIndex(graph, 'minimum')
