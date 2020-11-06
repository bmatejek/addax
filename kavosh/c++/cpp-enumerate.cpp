#include <stdio.h>
#include "cpp-graph.h"



int main(void)
{
    Graph graph = Graph("hemi-brain", true);

    graph.AddVertex(3);
    graph.AddVertex(4);

    graph.AddEdge(3, 4);

    return 1;
}
