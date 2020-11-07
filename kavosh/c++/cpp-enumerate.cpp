#include <stdio.h>
#include <ctime>
#include <unordered_map>
#include "cpp-graph.h"


long enumerated_subgraphs = 0;



std::vector<long> Validate(Graph *G,
                           std::unordered_set<long> &parents,
                           long u, std::unordered_set<long> &visited)
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
            if (u < w && visited.find(w) == visited.end()) {
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
                     std::unordered_map<long, std::unordered_set<long> > &S,
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
        enumerated_subgraphs += 1;

        // return a sorted tuple for this graph
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

    // make sure this vertex appears in the graph
    assert (G->vertices.find(u) != G->vertices.end());

    // keep track globally (through parameter passing) the vertices visited at higher enumeration steps
    std::unordered_set<long> visited = std::unordered_set<long>();

    // add the root vertex to the visited set
    visited.insert(u);

    // create the selection (first layer only ever has the root vertex)
    std::unordered_map<long, std::unordered_set<long> > S = std::unordered_map<long, std::unordered_set<long> >();
    S[0] = std::unordered_set<long>();
    S[0].insert(u);

    // enumerate all subgraphs of size k - 1 that contain the root u
    EnumerateVertex(G, u, S, k - 1, 1, visited);
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
    Graph *graph = ReadGraph("/home/bmatejek/Dropbox/motifs-temp/graphs/hemi-brain.graph");
    if (!graph) exit(-1);
    printf("Read graph in %0.2lf seconds.\n", (float)(clock() - start_time) / CLOCKS_PER_SEC);

    int k = 4;
    EnumerateSubgraphsSequentially(graph, k);

    return 1;
}
