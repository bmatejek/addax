import os
import networkx as nx



import matplotlib
import matplotlib.pyplot as plt



plt.style.use('seaborn')



from addax.analysis.certificates import ReadCertificates, ReadSummaryStatistics
from addax.kavosh.classify import ParseCertificate
from addax.utilities.dataIO import ReadGraph



def FindFrequentMotifs():
    """
    Find the top ten most frequent motifs for all of the datasets and draw them
    """
    # the datasets
    input_filenames = [
        'graphs/hemi-brain-minimum.graph.bz2',
        'graphs/C-elegans-timelapsed-01-minimum.graph.bz2',
        'graphs/C-elegans-timelapsed-02-minimum.graph.bz2',
        'graphs/C-elegans-timelapsed-03-minimum.graph.bz2',
        'graphs/C-elegans-timelapsed-04-minimum.graph.bz2',
        'graphs/C-elegans-timelapsed-05-minimum.graph.bz2',
        'graphs/C-elegans-timelapsed-06-minimum.graph.bz2',
        'graphs/C-elegans-timelapsed-07-minimum.graph.bz2',
        'graphs/C-elegans-timelapsed-08-minimum.graph.bz2',
        'graphs/C-elegans-sex-male-minimum.graph.bz2',
        'graphs/C-elegans-sex-hermaphrodite-minimum.graph.bz2',
    ]

    # the number of motifs to consider per dataset
    nmotifs = 5

    # go through each dataset and get the top ten most frequent motifs
    for input_filename in input_filenames:
        graph = ReadGraph(input_filename, vertices_only = True)

        # only consider k values of 4 and 5
        for k in [4, 5]:
            # no colors or community based
            nsubgraphs, _ = ReadSummaryStatistics(input_filename, k, False, False, False)
            certificates, _, _ = ReadCertificates(input_filename, k, False, False, False, 5)

            # create a new figure with two rows of 5
            fig, ax = plt.subplots(1, 5)

            ratio = 3
            fig.set_figheight(1.25 * ratio)
            fig.set_figwidth(5 * ratio)

            for iv, (certificate, subgraphs) in enumerate(sorted(certificates.items(), key = lambda x: x[1], reverse = True)):
                # vertex and edges are not colored, graph is directed
                nx_graph = ParseCertificate(k, certificate, False, False, True)

                # get the position for networkx
                pos = nx.circular_layout(nx_graph)

                row = iv // 5
                col = iv % 5

                # draw the network to matplotlib
                nx.draw_networkx_nodes(nx_graph, pos, ax=ax[col], node_size = 800, linewidths=2, node_shape = 'o', node_color='#e06666')
                nx.draw_networkx_edges(nx_graph, pos, ax=ax[col], edge_color = 'black', arrowsize=50, width=3)#, min_source_margin=-1000, min_target_margin=-1000)

                ax[col].set_title('{:0.2f}%'.format(100 * subgraphs / nsubgraphs), fontsize=42)
                ax[col].axis('off')

                # draw a boundary around the nodes
                ax[col].collections[0].set_edgecolor('#666666')


            # create the output directory if it doesn't exist
            output_directory = 'figures/{}'.format(graph.prefix, k)
            if not os.path.exists(output_directory):
                os.makedirs(output_directory, exist_ok = True)
            output_filename = '{}//motifs-{:03d}.png'.format(output_directory, k)

            plt.tight_layout()

            plt.savefig(output_filename)
            plt.clf()
