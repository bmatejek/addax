#ifndef __CPP_ENUMERATE_H__
#define __CPP_ENUMERATE_H__


// set global flags
void CppSetCommunityBased(bool input_community_based);
void CppSetWriteSubgraphs(bool input_write_subgraphs);

// enumeration functions
void CppEnumerateSubgraphsSequentially(const char *input_filename, short k);
void CppEnumerateSubgraphsFromNodes(const char *input_filename, short k, long *nodes, long nnodes, long output_suffix);



#endif
