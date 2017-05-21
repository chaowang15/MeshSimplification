#include "VertexQuadrics.h"
#include "tools.h"

VertexQuadrics::VertexQuadrics(VertexID id, uint dimVal) : target(dimVal)
{
	vertexId = id;
	edges = new vector<EdgeID>;
	vtxQuadrics = new MxQuadric(dimVal);
	initialQuadrics = new MxQuadric(dimVal);
	cluster = new set<VertexID>;
	cluster->insert(id);
	neighbors = new set<VertexID>;
	faces = new vector<FaceID>;
	clusterColor[0] = clusterColor[1] = clusterColor[2] = 0;
	clusterColor[3] = 1.0;
	cluster_color_flag = false;
	vertexQuality = 0.0;

	//isotropicQuadrics = new MxQuadric(dimVal);
}

VertexQuadrics::~VertexQuadrics()
{
	edges->clear();
	ClearVector(*edges);
	delete vtxQuadrics;
	delete initialQuadrics;
	cluster->clear();
	ClearSet(*cluster);
	neighbors->clear();
	ClearSet(*neighbors);
	faces->clear();
	ClearVector(*faces);

	//delete isotropicQuadrics;
}

void VertexQuadrics::clearAll()
{
	clear_all_edges();
	//ClearVector(*edges);
	
	vertexId = 0;
	vtxQuadrics->clear(0);
	cluster->clear();

	clear_all_faces();
	//ClearVector(*faces);
	cluster_color_flag = false;
}

bool VertexQuadrics::is_edge_included(EdgeID id)
{
	if (edges->empty())
	{
		return false;
	}
	else
	{
		for (EdgeID i = 0; i != edgeNumber; ++i)
		{
			if (edges->at(i) == id)
			{
				return true;
			}
		}
		return false;
	}
}

void VertexQuadrics::merge_cluster(VertexQuadrics *vtxQuad)
{
	for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
		iter != vtxQuad->cluster->end(); ++iter)
	{
		cluster->insert(*iter);
	}
	vtxQuad->cluster->clear();
}


void VertexQuadrics::addEdge(EdgeID id)
{
	edges->push_back(id);
	edgeNumber = edges->size();
}


void VertexQuadrics::addEdge(EdgeInfo *e)
{
	edges->push_back(e->no);
	edgeNumber = edges->size();
}

void VertexQuadrics::addFace(FaceInfo *info_f)
{
	faces->push_back(info_f->no);
	faceNumber = faces->size();
}

void VertexQuadrics::addFace(FaceID id)
{
	faces->push_back(id);
	faceNumber = faces->size();
}


void VertexQuadrics::delete_edge(EdgeID id)
{
	if (is_edge_included(id))
	{
		uint i = 0;
		for ( i = 0 ; i != edgeNumber; ++i)
		{
			if (edges->at(i) == id)
			{
				break;
			}
		}
		edges->erase(edges->begin()+i);
		edgeNumber = edges->size();
	}
}

void VertexQuadrics::delete_face(FaceID id)
{
	uint i = 0;
	for ( i = 0 ; i != faces->size(); ++i)
	{
		if (faces->at(i) == id)
		{
			break;
		}
	}
	faces->erase(faces->begin()+i);
	faceNumber = faces->size();
}
