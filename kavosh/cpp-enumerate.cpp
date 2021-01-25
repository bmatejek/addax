#include <pwd.h>
#include <stdio.h>
#include <ctime>
#include <map>
#include <string>
#include <nauty.h>
#include "cpp-nauty.h"
#include "cpp-graph.h"



static long enumerated_subgraphs = 0;              // the number of enumerated subgraphs identified
static NyGraph *nauty_graph;                       // graph object for canonical labeling
static std::map<std::string, long> certificates;   // map of certificates

static FILE *certificate_fp = NULL;             // file descriptor to write all certificates
static FILE *subgraph_fp = NULL;    // file descriptor to write all subgraphs



// global parameter flags
static bool VERTEX_COLORED = false;
static bool EDGE_COLORED = false;
static bool COMMUNITY_BASED = false;
static bool WRITE_SUBGRAPHS = false;



std::vector<long> Validate(Graph *G,
                           std::unordered_set<long> &parents,
                           long u,
                           std::unordered_set<long> &visited)
{
    /*
    Find the valid vertices for the next recursive level. Excludes elements with indices smaller
    than the root vertex and those already visited at a previous recursive level

    @param G: graph
    @param parents: selected vertices of last layer
    @param u: root vertex
    visited: a list of vertices already visited

    Returns a sorted list of valid vertices for the next level of analysis
    */

    // create a set of valid vertices at this level
    std::unordered_set<long> valid_vertices = std::unordered_set<long>();

    // iterate over all the parents of this layer
    for (std::unordered_set<long>::iterator it1 = parents.begin(); it1 != parents.end(); ++it1) {
        long v = *it1;

        // iterate over all the neighbors of this parent
        std::unordered_set<long> neighbors = G->vertices[v]->neighbors;
        for (std::unordered_set<long>::iterator it2 = neighbors.begin(); it2 != neighbors.end(); ++it2) {
            long w = *it2;

            // only consider neighbors that are in the same community if community based
            if (COMMUNITY_BASED && G->vertices[v]->community != G->vertices[w]->community) continue;

            // if the root vertex is less than the neighbor and the neighbor has not been visited
            // we use <= rather than < since u is always in visited as the S[0] entry
            // By using <=, we can enumerate all subgraphs with duplication by setting the enumeration indices to be non unique
            if (G->vertices[u]->enumeration_index <= G->vertices[w]->enumeration_index && visited.find(w) == visited.end()) {
                valid_vertices.insert(w);
                visited.insert(w);
            }
        }
    }

    // convert the set into a vector
    std::vector<long> valid_vertices_vector = std::vector<long>();
    for (std::unordered_set<long>::iterator it = valid_vertices.begin(); it != valid_vertices.end(); ++it) {
        valid_vertices_vector.push_back(*it);
    }

    return valid_vertices_vector;
}



void Combination(std::vector<long> &vertices,
                 std::vector<std::vector<long> > &combinations,
                 long nvertices,
                 long r,
                 long combo_index,
                 long combination[],
                 long vertex_index)
{
    /*
    Generate a new combination

    @param nvertices: size of input
    @param combinations: a vector to return the combinations
    @param r: size of combination
    @param index: current index within the new combination
    @param data: temporary array to store current combination
    @param vertex_index: index of current element in the input vector
    */

    // current combination is ready
    if (combo_index == r) {
        std::vector<long> new_combination = std::vector<long>();

        for (long iv = 0; iv < r; ++iv) {
            new_combination.push_back(combination[iv]);
        }

        combinations.push_back(new_combination);

        // return to previous recursion level
        return;
    }

    // there are no more elements to add to the combination
    if (vertex_index >= nvertices) return;

    // vertex_index is included, put next combo_index at next location
    combination[combo_index] = vertices[vertex_index];
    Combination(vertices, combinations, nvertices, r, combo_index + 1, combination, vertex_index + 1);

    // current is excluded, replace it with the next vertex
    // we do not increment the combo_index, but the vertex_index
    Combination(vertices, combinations, nvertices, r, combo_index, combination, vertex_index + 1);

}



std::vector<std::vector<long> > Combinations(std::vector<long> &vertices, long r)
{
    /*
    Create a vector of combinations of length r from vertices

    @param vertices: a vector of vertex indices from which to generate combinations
    @param r: the size of the combinations to generate
    */

    std::vector<std::vector<long> > combinations = std::vector<std::vector<long> >();

    // get the number of elements in the vector
    long nvertices = vertices.size();
    // a temporary array to store combinations as they form
    long combination[r];

    // begin population all combinations
    Combination(vertices, combinations, nvertices, r, 0, combination, 0);

    return combinations;
}



void EnumerateVertex(Graph *G,
                     long u,
                     std::map<long, std::unordered_set<long> > &S,
                     short rem,
                     short i,
                     std::unordered_set<long> &visited)
{
    /*
    Enumerate all subgraphs of size rem that contain vertices in S[0] ... S[i - 1]

    @param G: graph
    @param u: root vertex
    @param S: selection (S = {S_0, S_i, ... S_{k - 1}}) is an array of the set of all S_i
    @param rem: number of remaining vertices to be selected
    @param i: current depth of the tree
    @param visited: a list of vertices already visited

    Returns a generator that continually gives the next subgraph that contains S[0] ... S[i - 1] of
    the appropriate size (k)
    */

    // if there are no remaining vertices to add, subgraph size limit reached
    if (!rem) {
        // the set of vertices in S[0] to S[i - 1] contain the subgraph
        // note that the sets in S[i] ... S[k] do not belong to the subgraph but a previous iteration

        // create a mapping from vertices to linear indices
        std::map<long, long> index_to_vertex = std::map<long, long>();
        // initialize a colorind mapping regardless of if vertex coloring exists,
        // will not be populated for uncolored graphs
        // maps vertex color -> list of subgraph indices with that color
        std::map<long, std::vector<long> > coloring = std::map<long, std::vector<long> >();
        std::map<long, long> index_to_coloring = std::map<long, long>();
        short index = 0;

        for (short level = 0; level <= i - 1; ++level) {
            for (std::unordered_set<long>::iterator it = S[level].begin(); it != S[level].end(); ++it, ++index) {
                // map this vertex to an index between 0 and k - 1
                index_to_vertex[index] = *it;


                if (VERTEX_COLORED) {
                    // get the color for this vertex
                    long color = G->vertices[*it]->color;

                    // if this color is not yet seen, create a new vector for these colors
                    if (coloring.find(color) == coloring.end()) {
                        coloring[color] = std::vector<long>();
                    }

                    // add this index to the coloring
                    coloring[color].push_back(index);

                    // create a mapping from the indices to the coloring
                    index_to_coloring[index] = color;
                }
            }
        }

        // the size of the motif
        long k = index;
        for (long out_index = 0; out_index < k; ++out_index) {
            long out_vertex = index_to_vertex[out_index];
            for (long in_index = 0; in_index < k; ++in_index) {
                long in_vertex = index_to_vertex[in_index];

                // there is an edge from out_vertex to in_vertex
                if (G->vertices[out_vertex]->outgoing_neighbors.find(in_vertex) != G->vertices[out_vertex]->outgoing_neighbors.end()) {
                    ADDELEMENT((GRAPHROW(nauty_graph->matrix, out_index, nauty_graph->no_setwords)), in_index);
                }
            }
        }

        /*
        If nauty_graph->options->defaultptn = FALSE the vertices have colors. The colors are determined by the
        arrays int *lab and int *ptn. If nauty_graph->options->getcanon = TRUE, nauty_graph->lab will list the vertices in g
        in the otder in which they need to be relabelled..
        */

        // set *ptn and *lab for graph coloring
        if (VERTEX_COLORED) {
            // keep a linear index for the permutation arrays (lab, ptn)
            long permutation_index = 0;

            // go through the colors in order
            for (std::map<long, std::vector<long> >::iterator it = coloring.begin(); it != coloring.end(); ++it) {
                for (unsigned long iv = 0; iv < it->second.size(); ++iv) {
                    long vertex_index = it->second[iv];

                    // set the labeling for this permutation index to this vertex between (0, k - 1)
                    nauty_graph->lab[permutation_index] = vertex_index;
                    // all values of ptn should be one except for the end of the coloring
                    // set all to one here, and after this loop set the previous index to 0
                    nauty_graph->ptn[permutation_index] = 1;

                    permutation_index += 1;
                }

                // the last element input should have a value of 0 since it ended the coloring
                nauty_graph->ptn[permutation_index - 1] = 0;
            }
        }

        // call the dense version of nauty
        densenauty(
                    nauty_graph->matrix,
                    nauty_graph->lab,
                    nauty_graph->ptn,
                    nauty_graph->orbits,
                    nauty_graph->options,
                    nauty_graph->stats,
                    nauty_graph->no_setwords,
                    nauty_graph->no_vertices,
                    nauty_graph->cmatrix
                );

        long nbytes = nauty_graph->no_vertices * nauty_graph->no_setwords * sizeof(setword);
        unsigned char certificate_bytes[nbytes];
        memcpy(&(certificate_bytes[0]), &(nauty_graph->cmatrix[0]), nbytes);

        // construct a string certificate - need to iteratively add chars to the string
        // to allow for null characters which are common in the cmatrix
        std::string certificate = std::string();
        for (long iv = 0; iv < nbytes; ++iv) {
            certificate.push_back(certificate_bytes[iv]);
        }

        // add the vertex coloring to the certificate
        if (VERTEX_COLORED) {
            // go through all vertices in the small subgraph
            for (long iv = 0; iv < k; ++iv) {
                // the value of int *lab after the call to nauty returns the vertices of
                // g in order in which they need to be relabelled to give the canonical graph
                // so lab[iv] gives the original index that maps to this location in the
                // canonical labeling, and index_to_coloring[labl[iv]] gives the color of that vertex
                long color = index_to_coloring[nauty_graph->lab[iv]];

                // convert the long into bytes
                short nbytes_per_long = 8;
                for (long ib = 0; ib < nbytes_per_long; ++ib) {
                    // first remove the bits in the previous bytes
                    // then remove the bits lower order than this
                    unsigned char byte = (color << 8 * ib) >> 56;
                    certificate.push_back(byte);
                }
            }
        }
        // add the edge coloring to the certificate
        if (EDGE_COLORED) {
            // go through all pairs of vertices in the small subgraph
            for (long iv1 = 0; iv1 < k; ++iv1) {
                long vertex_one = index_to_vertex[nauty_graph->lab[iv1]];

                for (long iv2 = 0; iv2 < k; ++iv2) {
                    long vertex_two = index_to_vertex[nauty_graph->lab[iv2]];

                    // get the edge color from vertex_one to vertex_two (default = -1)
                    long color = -1;
                    if (G->edges.find(std::pair<long, long>(vertex_one, vertex_two)) != G->edges.end()) {
                        color = G->edges[std::pair<long, long>(vertex_one, vertex_two)]->color;
                    }
                    
                    // convert the long into bytes
                    short nbytes_per_long = 8;
                    for (long ib = 0; ib < nbytes_per_long; ++ib) {
                        // first remove the bits in the previous bytes
                        // then remove the bits lower order than this
                        unsigned char byte = (color << 8 * ib) >> 56;
                        certificate.push_back(byte);
                    }
                }
            }
        }

        // add this enumerated subgraph to the grouping of certificates
        if (certificates.find(certificate) == certificates.end()) {
            certificates[certificate] = 1;
        }
        else {
            certificates[certificate] += 1;
        }

        // write the subgraph and labeling to disk if required
        if (WRITE_SUBGRAPHS) {
            // write the certificate for this subgraph
            const char *certificate_chars = certificate.c_str();

            for (unsigned long iv = 0; iv < certificate.length(); ++iv) {
                fprintf(subgraph_fp, "%02x", (unsigned char) certificate_chars[iv]);
            }
            // create separation for certificate to vertices
            fprintf(subgraph_fp, ": ");

            // go through all vertices in the small subgraph
            for (long iv = 0; iv < k; ++iv) {
                // as above, the value of int *lab after the call to nauty returns the vertices of
                // g in order in which they need to be relablled to give the canonical graph
                // so lab[iv] gives the original index that maps to this location in the canonical
                // labeling, and index_to_vertex[lab[iv]] gives the original vertex value
                long vertex = index_to_vertex[nauty_graph->lab[iv]];

                fprintf(subgraph_fp, "%ld ", vertex);
            }

            // create a new line between this and the next subgraph
            fprintf(subgraph_fp, "\n");
        }

        // clear the graph
        EMPTYGRAPH(nauty_graph->matrix, nauty_graph->no_setwords, nauty_graph->no_vertices);

        // update the total enumerated subgraphs
        enumerated_subgraphs += 1;

        // final recursion limit reached for this subgraph
        return;
    }

    // get the valid vertices
    std::vector<long> valid_vertices = Validate(G, S[i - 1], u, visited);

    // the max number of vertices for this layer is the minimum of the number of children or the remaining
    short n_i = std::min(valid_vertices.size(), (unsigned long)rem);

    for (short k_i = 1; k_i <= n_i; ++k_i) {
        // get all combinations of size k_i for the valid vertices
        std::vector<std::vector<long> > combinations = Combinations(valid_vertices, k_i);

        for (unsigned long iv1 = 0; iv1 < combinations.size(); ++iv1) {
            std::vector<long> combination = combinations[iv1];

            // update the set S value at this level with this combination
            S[i] = std::unordered_set<long>();
            for (unsigned long iv2 = 0; iv2 < combination.size(); ++iv2) {
                S[i].insert(combination[iv2]);
            }

            EnumerateVertex(G, u, S, rem - k_i, i + 1, visited);
        }
    }

    // remove all the valid vertices from the list of visited
    // finished all subgraphs for this level and proceed back to level above (closer to root)
    for (unsigned long it = 0; it < valid_vertices.size(); ++it) {
        visited.erase(valid_vertices[it]);
    }

}



void EnumerateSubgraphsFromNode(Graph *G, short k, long u)
{
    /*
    Enumerate all subgraphs of a given motif size rooted at a given vertex

    @param G: graph
    @param k: motif size
    @param u: root vertex index

    Returns a generator that continually gives the next subgraph in the graph rooted as this vertex
    */
    // start statistics
    unsigned int start_time = clock();

    enumerated_subgraphs = 0;

    // create an empty NyGraph (Nauty Graph) data structure
    if (VERTEX_COLORED) nauty_graph = new NyGraph(k, true);
    else nauty_graph = new NyGraph(k, false);


    // create an empty certificates dictionary
    certificates = std::map<std::string, long>();

    // make sure this vertex appears in the graph
    assert (G->vertices.find(u) != G->vertices.end());

    // keep track globally (through parameter passing) the vertices visited at higher enumeration steps
    std::unordered_set<long> visited = std::unordered_set<long>();

    // add the root vertex to the visited set
    visited.insert(u);

    // create the selection (first layer only ever has the root vertex)
    std::map<long, std::unordered_set<long> > S = std::map<long, std::unordered_set<long> >();
    S[0] = std::unordered_set<long>();
    S[0].insert(u);

    // enumerate all subgraphs of size k - 1 that contain the root u
    EnumerateVertex(G, u, S, k - 1, 1, visited);

    // don't include any I/O time in the total time
    float total_time = (float) (clock() - start_time) / CLOCKS_PER_SEC;

    for (std::map<std::string, long>::iterator it = certificates.begin(); it != certificates.end(); ++it) {
        const char *certificate = it->first.c_str();

        for (unsigned long iv = 0; iv < it->first.length(); ++iv) {
            fprintf(certificate_fp, "%02x", (unsigned char) certificate[iv]);
        }
        fprintf(certificate_fp, ": %ld\n", it->second);
    }

    // clear the certificates
    certificates.clear();

    // free memory
    delete nauty_graph;

    // print statistics
    fprintf(certificate_fp, "Enumerated %ld subgraphs for node %ld in %0.2f seconds.\n", enumerated_subgraphs, u, total_time);
}



void CppSetVertexColored(bool input_vertex_colored) {
    // set the vertex colored flag
    VERTEX_COLORED = input_vertex_colored;
}



void CppSetEdgeColored(bool input_edge_colored) {
    // set the edge colored flag
    EDGE_COLORED = input_edge_colored;
}



void CppSetCommunityBased(bool input_community_based) {
    // set the community based flag
    COMMUNITY_BASED = input_community_based;
}



void CppSetWriteSubgraphs(bool input_write_subgraphs) {
    // set the write subgraphs flag
    WRITE_SUBGRAPHS = input_write_subgraphs;
}



void CppEnumerateSubgraphsSequentially(const char *input_filename, const char *temp_directory, short k)
{
    // read the input file
    Graph *G = ReadBZ2Graph(input_filename);
    if (!G) exit(-1);

    // home directory
    passwd *pw = getpwuid(getuid());

    // create a new file for writing the certificates
    char output_filename[4096];
    snprintf(output_filename, 4096, "%s/motifs/%s/certificates/motif-size-%03d-certificates.txt", pw->pw_dir, temp_directory, k);

    // open the file
    certificate_fp = fopen(output_filename, "w");
    if (!certificate_fp) { fprintf(stderr, "Failed to open %s\n", output_filename); exit(-1); }

    // create a new file for writing subgraphs if needed
    if (WRITE_SUBGRAPHS) {
        char subgraph_filename[4096];
        snprintf(subgraph_filename, 4096, "%s/motifs/%s/subgraphs/motif-size-%03d-subgraphs.txt", pw->pw_dir, temp_directory, k);

        // open the file
        subgraph_fp = fopen(subgraph_filename, "w");
        if (!subgraph_fp) { fprintf(stderr, "Failed to open %s\n", subgraph_filename); exit(-1); }
    }

    // iterate over all vertices in the graph
    for (std::map<long, Vertex *>::iterator it = G->vertices.begin(); it != G->vertices.end(); ++it) {
        long u = it->first;
        EnumerateSubgraphsFromNode(G, k, u);
    }

    // close the files
    fclose(certificate_fp);
    if (WRITE_SUBGRAPHS) fclose(subgraph_fp);

    // free memory
    delete G;
}



void CppEnumerateSubgraphsFromNodes(const char *input_filename, const char *temp_directory, short k, long *nodes, long nnodes, long output_suffix)
{
    // read the input file
    Graph *G = ReadBZ2Graph(input_filename);
    if (!G) exit(-1);

    // home directory
    passwd *pw = getpwuid(getuid());

    // create a new file for writing the certificates
    char output_filename[4096];
    snprintf(output_filename, 4096, "%s/motifs/%s/certificates/motif-size-%03d-output-%08ld-certificates.txt", pw->pw_dir, temp_directory, k, output_suffix);

    // open the file
    certificate_fp = fopen(output_filename, "w");
    if (!certificate_fp) { fprintf(stderr, "Failed to open %s\n", output_filename); exit(-1); }

    // create a new file for writing subgraphs if needed
    if (WRITE_SUBGRAPHS) {
        char subgraph_filename[4096];
        snprintf(subgraph_filename, 4096, "%s/motifs/%s/subgraphs/motif-size-%03d-output-%08ld-subgraphs.txt", pw->pw_dir, temp_directory, k, output_suffix);

        // open the file
        subgraph_fp = fopen(subgraph_filename, "w");
        if (!subgraph_fp) { fprintf(stderr, "Failed to open %s\n", subgraph_filename); exit(-1); }
    }

    for (long iv = 0; iv < nnodes; ++iv) {
        EnumerateSubgraphsFromNode(G, k, nodes[iv]);
    }

    // close the files
    fclose(certificate_fp);
    if (WRITE_SUBGRAPHS) fclose(subgraph_fp);

    // free memory
    delete G;
}
