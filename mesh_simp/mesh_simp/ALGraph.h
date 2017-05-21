/************************************************************************/  
/* 
*
* The Adjacent List format of a graph and its most related functions
* 
* Author: Chao Wang
* Date: 12.18.2013
*/  
/************************************************************************/

#ifndef ALGRAPH_H
#define ALGRAPH_H

#include <iostream>
#include <vector>
#include "global.h"
#include "tools.h"

using namespace std;


// Edge Node
typedef struct node 
{  
	VertexType index;	// index
	struct node *next;  // pointer to another EdgeNode
}EdgeNode;

// Graph in Adjacent List Format
class ALGraph
{
public:
	ALGraph(int numNode);
	~ALGraph();

	// Clear all variables
	void clearAll();

	// Initialization
	void init();

	// Delete a node
	void delete_node(VertexID i);

	// Check if a node is already visited
	inline bool isVisited(VertexID i){return visited->at(i);}

	// Clear the node number as 0 (Must be used before a new search)
	inline void clear_node_number(){countNum = 0;}

	// Add an edge (v1,v2) into the graph
	void addEdge(VertexType v1, VertexType v2);

	// Get number of nodes counted after a search
	inline int get_node_number_in_search(){return countNum;} 

	// Create an ALGraph from input faces
	void createALGraph(vector<int> *vecFaceIndex, int *faceArray);

	// Depth-First Search starting from one node
	// Note: (1) The number of nodes for this DFS search is saved in "countNum"
	//		 (2) The node indexes for this search are saved in "vecVtxTraverse" in the DFS order
	void DFS(int nodeID);

	// Traverse all the nodes of the graph using DFS
	// Note: must clear the "countNum" using get_node_number_in_search() before a new search
	void DFSTraverse();

	// Depth-First Search starting from one node (recursive version) (non-recursive version)
	void DFS_NonRecursive(int nodeID);

	// Traverse all the nodes of the graph using non-recursive DFS
	void DFSTraverse_NonRecursive();

	// Get a connected component starting from one node and return the number of this component
	int GetConnectedComponent(int nodeID);

	// Get all node numbers of connected component starting from nodes in "vecStartNodeID" and save 
	// the results into the "vecContCompNum"
	void GetAllConnectedComponentNumber(vector<int> *vecStartNodeID, vector<int> *vecContCompNum);

public:
	vector<EdgeNode*> *adjlist;			// adjacent list
	vector<bool> *visited;				// visited flag
	int numVtx;							// number of VERTICES
	int countNum;
	vector<VertexType> *vecVtxTraverse;	// VERTICES obtained in a search process
};



#endif


