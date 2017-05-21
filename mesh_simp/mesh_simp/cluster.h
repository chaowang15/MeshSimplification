/************************************************************************/  
/* 
* Cluster.h
*	Clusters used in QEM
* 
* Author: Chao Wang
* Date: 12.18.2013
*/  
/************************************************************************/

#ifndef CLUSTER_H
#define CLUSTER_H

#include <set>
#include <ctime>
#include "Eigen/Dense"
#include "MxGeoPrims.h"
using namespace std;

typedef unsigned int uint;

// Point Cluster
class PCluster
{
public:
	set<MxVertexID> * point_clusters;	// all clusters (each may contain different VERTICES)
	long * cluster_belong;			// the cluster index each point belongs to
	int numClusters;				// number of clusters

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
	PCluster(int n)
	{
		numClusters = n;
		point_clusters = new set<MxVertexID>[n];	// initially n clusters for n VERTICES
		cluster_belong = new long[n];			// each point belongs to the cluster itself
		for (int i=0; i< n; i++)
		{
			point_clusters[i].clear();
		}

		//at first, each domain only contain and belong to itself
		for (int i=0; i< n; i++)
		{
			point_clusters[i].insert(i);
			cluster_belong[i] = i;
		}
		
	}

	~PCluster()
	{
		for (uint i = 0; i != point_clusters->size(); ++i)
		{
			point_clusters[i].clear();
		}		
		delete [] point_clusters;
		delete [] cluster_belong;
	}

	// Merge all VERTICES in cluster v2 into cluster v1
	void mergeCluster(int v1, int v2)
	{		
		//update containing and belonging information
		for (set<MxVertexID>::iterator it = point_clusters[v2].begin(); it !=point_clusters[v2].end(); it++)
		{
			point_clusters[v1].insert(*it);
			cluster_belong[*it] = v1;
		}
		point_clusters[v2].clear();	
	}
};

#endif