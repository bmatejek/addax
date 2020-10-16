import struct



from addax.utilities.dataIO import ReadGraph



def ReadCommunities(input_filename):
    communities = {}

    with open(input_filename, 'rb') as fd:
        nvertices, = struct.unpack('q', fd.read(8))
        for _ in range(nvertices):
            vertex_id, community, = struct.unpack('qq', fd.read(16))

            communities[vertex_id] = community

    return communities



def CompareCommunities(ground_truth_community, proposed_community):
    vertices = sorted(ground_truth_community.keys())
    nvertices = len(vertices)

    TP = 0
    FN = 0
    FP = 0
    TN = 0

    for iv1 in range(nvertices):
        vertex_one = vertices[iv1]

        # get the ground truth and proposal for vertex one
        ground_truth_one = ground_truth_community[vertex_one]
        proposed_one = proposed_community[vertex_one]


        for iv2 in range(nvertices):
            vertex_two = vertices[iv2]

            # get the ground truth and proposal for vertex two
            ground_truth_two = ground_truth_community[vertex_two]
            proposed_two = proposed_community[vertex_two]

            ground_truth = (ground_truth_one == ground_truth_two)
            proposed = (proposed_one == proposed_two)

            if ground_truth and proposed: TP += 1
            elif ground_truth and not proposed: FN += 1
            elif not ground_truth and proposed: FP += 1
            else: TN += 1

    print ('True Positives: {}'.format(TP))
    print ('False Negatives: {}'.format(FN))
    print ('False Positives: {}'.format(FP))
    print ('True Negatives: {}'.format(TN))

    print ('Accuracy: {:0.2f}%'.format((TP + TN) / (TP + FN + FP + TN)))



def EvaluateLouvainCommunities():
    # read the graph and louvain communities
    graph = ReadGraph('graphs/hemi-brain.graph.bz2')
    ground_truth_communities = graph.Communities()

    louvain_communities = ReadCommunities('communities/hemi-brain-louvain.dict')

    print (len(set(louvain_communities.values())))
    return

    CompareCommunities(ground_truth_communities, louvain_communities)
