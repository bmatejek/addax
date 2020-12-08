#ifndef __CPP_ENUMERATE_H__
#define __CPP_ENUMERATE_H__


void CppSetCommunityBased(bool input_community_based);
void CppEnumerateSubgraphsSequentially(const char *input_filename, short k);
void CppEnumerateSubgraphsFromNodes(const char *input_filename, short k, long *nodes, long nnodes, long output_suffix);



#endif
