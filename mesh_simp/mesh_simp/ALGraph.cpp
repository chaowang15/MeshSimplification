#include "ALGraph.h"
#include <stack>

ALGraph::ALGraph(int numNode)
{
	numVtx = numNode;
	init();
}

ALGraph::~ALGraph()
{
	clearAll();
}

void ALGraph::init()
{
	adjlist = new vector<EdgeNode*>;
	adjlist->resize(numVtx, NULL);
	visited = new vector<bool>;
	visited->resize(numVtx, false);
	countNum = 0;
	vecVtxTraverse = new vector<VertexType>;
}

void ALGraph::clearAll()
{
	if (!adjlist->empty())
	{
		for (int i = 0; i != adjlist->size(); ++i)
		{
			delete_node(i);
		}
		adjlist->clear();
		ClearVector(*adjlist);
	}
	visited->clear();
	ClearVector(*visited);
	vecVtxTraverse->clear();
	ClearVector(*vecVtxTraverse);
}

void ALGraph::delete_node(VertexID i)
{
	if (adjlist->at(i) == NULL)
	{
		return;
	}
	EdgeNode *edgeNow = adjlist->at(i);
	EdgeNode *edgeTemp;
	while (edgeNow->next != NULL)
	{
		edgeTemp = edgeNow->next;
		delete edgeNow;
		edgeNow = edgeTemp;
	}
	delete edgeNow;
}

void ALGraph::addEdge(VertexType v1, VertexType v2)
{
	EdgeNode *edgeNodeNow = adjlist->at(v1);
	if (edgeNodeNow != NULL)
	{
		bool flag = true;
		while (true)
		{
			if (edgeNodeNow->index == v2)
			{
				flag = false;
				break;
			}
			
			if (edgeNodeNow->next == NULL)
				break;

			edgeNodeNow = edgeNodeNow->next;
		}
		if (flag)
		{
			EdgeNode *edge1 = new EdgeNode;
			edge1->index = v2;
			edge1->next = NULL;
			edgeNodeNow->next = edge1;
		}
	}
	else
	{
		EdgeNode *edge1 = new EdgeNode;
		edge1->index = v1;
		EdgeNode *edge2 = new EdgeNode;
		edge2->index = v2;
		edge2->next = NULL;
		edge1->next = edge2;
		adjlist->at(v1) = edge1;
	}
}

void ALGraph::createALGraph(vector<int> *vecFaceIndex, int *faceArray)
{
	for (int i = 0; i != vecFaceIndex->size(); ++i)
	{
		int faceID = vecFaceIndex->at(i);
		int *f = &(faceArray[3*faceID]);
		for (int j = 0; j != 3; ++j)
		{
			if (adjlist->at(f[j]) == NULL)
			{
				EdgeNode *edge0 = new EdgeNode;
				edge0->index = f[j];
				EdgeNode *edge1 = new EdgeNode;
				edge1->index = f[(j+1)%3];
				EdgeNode *edge2 = new EdgeNode;
				edge2->index = f[(j+2)%3];
				edge2->next = NULL;
				edge1->next = edge2;
				edge0->next = edge1;
				adjlist->at(f[j]) = edge0;
			}
			else
			{
				for (int k = 1; k <= 2; ++k)
				{
					// Add f[j+1]
					EdgeNode *edgeNow = adjlist->at(f[j]);
					int nodeID = f[(j+k)%3];
					int flag = true;
					do
					{// search all adjacent nodes of f[j]
						edgeNow = edgeNow->next;
						if (edgeNow->index == nodeID)
						{// If the node f[j+1] has already been saved in the adjacent list of f[j], exit
							flag = false;
							break;
						}
					}while (edgeNow->next != NULL);
					if (flag)
					{
						EdgeNode *edge1 = new EdgeNode;
						edge1->index = nodeID;
						edge1->next = NULL;
						edgeNow->next = edge1;
					}
				}	
			}
		}
	}
}

void ALGraph::DFS(int nodeID)
{
	EdgeNode *edgeNow = adjlist->at(nodeID);
	if (edgeNow != NULL)
	{
		if (!visited->at(edgeNow->index))
		{
			visited->at(edgeNow->index) = true;
			countNum++;
			vecVtxTraverse->push_back(edgeNow->index);
			//cout<<edgeNow->index<<endl;
			while (edgeNow->next != NULL)
			{
				edgeNow = edgeNow->next;
				if (!visited->at(edgeNow->index))
				{
					DFS(edgeNow->index);
				}
			}
		}
	}
}

void ALGraph::DFSTraverse()
{
	if (!visited->empty())
	{
		visited->clear();
	}
	visited->resize(numVtx, false);
	for (int i = 0; i != numVtx; ++i)
	{
		EdgeNode *edgeNow = adjlist->at(i);
		if (edgeNow != NULL)
		{
			if (!visited->at(edgeNow->index))
			{
				DFS(edgeNow->index);
				cout<<"Number of nodes for this connected component:"<<countNum<<endl;
				clear_node_number();
			}
		}
	}
}

void ALGraph::DFS_NonRecursive(int nodeID)
{
	EdgeNode *edgeNow = adjlist->at(nodeID);
	if (edgeNow != NULL)
	{
		if (!visited->at(edgeNow->index))
		{
			visited->at(edgeNow->index) = true;
			countNum++;
			vecVtxTraverse->push_back(edgeNow->index);

			stack<VertexID> stack_dfs;
			stack_dfs.push(nodeID);
			while (!stack_dfs.empty())
			{
				VertexID v_top = stack_dfs.top();
				edgeNow = adjlist->at(v_top);
				bool flag = true;
				while (edgeNow->next != NULL)
				{
					edgeNow = edgeNow->next;
					if (!visited->at(edgeNow->index))
					{
						countNum++;
						visited->at(edgeNow->index) = true;
						stack_dfs.push(edgeNow->index);
						vecVtxTraverse->push_back(edgeNow->index);
						flag = false;
						break;
					}
				}
				if (flag)
				{
					stack_dfs.pop();
				}
			}
			stack<VertexID> empty;
			std::swap( stack_dfs, empty );
		}
	}
}

void ALGraph::DFSTraverse_NonRecursive()
{
	if (!visited->empty())
	{
		visited->clear();
	}
	visited->resize(numVtx, false);
	for (int i = 0; i != numVtx; ++i)
	{
		EdgeNode *edgeNow = adjlist->at(i);
		if (edgeNow != NULL)
		{
			if (!visited->at(edgeNow->index))
			{
				DFS_NonRecursive(edgeNow->index);
				cout<<"Number of nodes for this connected component:"<<countNum<<endl;
				countNum = 0;
			}
		}
	}
}

int ALGraph::GetConnectedComponent(int nodeID)
{
	if (!visited->empty())
	{
		visited->clear();
	}
	if (!vecVtxTraverse->empty())
	{
		vecVtxTraverse->clear();
	}
	visited->resize(numVtx, false);
	EdgeNode *edgeNow = adjlist->at(nodeID);
	countNum = 0;
	if (edgeNow != NULL)
	{
		if (!visited->at(edgeNow->index))
		{
			DFS(edgeNow->index);
		}
	}
	return countNum;
}

void ALGraph::GetAllConnectedComponentNumber(vector<int> *vecStartNodeID, vector<int> *vecContCompNum)
{
	if (!visited->empty())
	{
		visited->clear();
	}
	if (!vecVtxTraverse->empty())
	{
		vecVtxTraverse->clear();
	}
	visited->resize(numVtx, false);
	for (int i = 0; i != numVtx; ++i)
	{
		EdgeNode *edgeNow = adjlist->at(i);
		if (edgeNow != NULL)
		{
			if (!visited->at(edgeNow->index))
			{
				countNum = 0;
				DFS(edgeNow->index);
				vecStartNodeID->push_back(edgeNow->index);
				vecContCompNum->push_back(countNum);
			}
		}
	}
}