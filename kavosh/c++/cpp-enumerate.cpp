#include <stdio.h>
#include <ctime>
#include <map>
#include <string>
#include <nauty.h>
#include "cpp-nauty.h"
#include "cpp-graph.h"



static long enumerated_subgraphs = 0;
static NyGraph *nauty_graph;
static std::map<std::string, long> certificates;



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

            // if the root vertex is less than the neighbor and the neighbor has not been visited
            if (G->vertices[u]->enumeration_index < G->vertices[w]->enumeration_index && visited.find(w) == visited.end()) {
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
        short index = 0;

        for (short level = 0; level <= i - 1; ++level) {
            for (std::unordered_set<long>::iterator it = S[level].begin(); it != S[level].end(); ++it, ++index) {
                // map this vertex to an index between 0 and k - 1
                index_to_vertex[index] = *it;
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
        array int *ptn. If nauty_graph->options->getcanon = TRUE, nauty_graph->lab will list the vertices in g
        in the otder in which they need to be relabled.
        */

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

        enumerated_subgraphs += 1;

        // clear the graph
        EMPTYGRAPH(nauty_graph->matrix, nauty_graph->no_setwords, nauty_graph->no_vertices);

        long nbytes = nauty_graph->no_vertices * nauty_graph->no_setwords * sizeof(setword);
        unsigned char certificate_bytes[nbytes];
        memcpy(&(certificate_bytes[0]), &(nauty_graph->cmatrix[0]), nbytes);

        // construct a string certificate - need to iteratively add chars to the string
        // to allow for null characters which are common in the cmatrix
        std::string certificate = std::string();
        for (long iv = 0; iv < nbytes; ++iv) {
            certificate.push_back(certificate_bytes[iv]);
        }

        // add this enumerated subgraph to the grouping of certificates
        if (certificates.find(certificate) == certificates.end()) {
            certificates[certificate] = 1;
        }
        else {
            certificates[certificate] += 1;
        }

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

        for (long iv1 = 0; iv1 < combinations.size(); ++iv1) {
            std::vector<long> combination = combinations[iv1];

            // update the set S value at this level with this combination
            S[i] = std::unordered_set<long>();
            for (long iv2 = 0; iv2 < combination.size(); ++iv2) {
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

    enumerated_subgraphs = 0;

    // create an empty NyGraph (Nauty Graph) data structure
    nauty_graph = new NyGraph(k);

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

    char output_filename[128];
    snprintf(output_filename, 128, "/home/brian/motifs/temp/hemi-brain-c++/motif-size-%03d-node-%016ld-certificates.txt", k, u);

    // open file
    FILE *fp = fopen(output_filename, "w");
    if (!fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

    for (std::map<std::string, long>::iterator it = certificates.begin(); it != certificates.end(); ++it) {
        const char *certificate = it->first.c_str();

        for (long iv = 0; iv < it->first.length(); ++iv) {
            fprintf(fp, "%02x", (unsigned char) certificate[iv]);
        }
        fprintf(fp, ": %ld\n", it->second);
    }

    // close file
    fclose(fp);

    // clear the certificates
    certificates.clear();

    // free memory
    delete nauty_graph;
}



void EnumerateSubgraphsSequentially(Graph *G, short k)
{
    /*
    Enumerate all subgraphs of a given size in the graph

    @param G: graph
    @param k: motif size

    Returns a generator that continually gives the next subgraph in the graph
    */

    // iterate over all vertices in the graph
    for (std::map<long, Vertex *>::iterator it = G->vertices.begin(); it != G->vertices.end(); ++it) {
        unsigned int start_time = clock();

        long u = it->first;
        EnumerateSubgraphsFromNode(G, k, u);

        printf("Enumerated %ld subgraphs for node %ld in %0.2f seconds.\n", enumerated_subgraphs, u, (float)(clock() - start_time) / CLOCKS_PER_SEC);
    }
}



int main(void)
{
    unsigned int start_time = clock();
    //Graph *graph = ReadGraph("/home/brian/Dropbox/motifs-temp/graphs/hemi-brain.graph");
    Graph *graph = ReadBZ2Graph("/home/brian/motifs/graphs/hemi-brain.graph.bz2");

    if (!graph) exit(-1);
    printf("Read graph in %0.2lf seconds.\n", (float)(clock() - start_time) / CLOCKS_PER_SEC);

    int k = 3;
    EnumerateSubgraphsSequentially(graph, k);

    // free memory
    delete graph;

    return 1;
}
