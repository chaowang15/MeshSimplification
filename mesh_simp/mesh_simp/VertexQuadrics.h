
/************************************************************************/
/* 
*	VertexQuadrics.h
*		Define vertex quadrics class for each vertex
*
*	Author: Chao Wang
*	Date: 2.25.2013
*/
/************************************************************************/

#ifndef VERTEXQUADRICS_H
#define VERTEXQUADRICS_H
#include <vector>
#include <set>
#include "MxHeap.h"
#include "global.h"
#include "MxQMetric.h"


using namespace std;

//typedef Vector4d VectorType;
//typedef Matrix4d MatrixType;
typedef MxVector VectorType;

class FaceInfo
{
public:
	VertexID v[3];			// three vertex indexes
	FaceID no;				// face No.
	bool flag_exist;		// true if it exists, and false if it is deleted

public:
	FaceInfo(){v[0]=v[1]=v[2]=no=0;flag_exist=true;}
	FaceInfo(VertexID a, VertexID b, VertexID c, FaceID nno=0): no(nno), flag_exist(true)
	{v[0] = a; v[1] = b; v[2] = c;}
	FaceInfo(int t[3], FaceID nno=0): no(nno), flag_exist(true)
	{v[0] = t[0]; v[1] = t[1]; v[2] = t[2];}
	// Update the vertex id a to b
	inline void update_vertex(VertexID a, VertexID b){
		if (flag_exist)
			for (uint i = 0; i != 3; ++i)
				if (a == v[i])
					v[i] = b;
	}
	// Delete this face by setting its flag as false
	inline void delete_face(){flag_exist=false;}
	inline bool is_exist(){return flag_exist;}
};


// Edge information for each edge
class EdgeInfo : public MxHeapable
{
public:
	VertexID v1, v2;		// two endpoint indexes of this edge
	double quadError;		// quadrics error for this edge
	EdgeID no;				// edge No.
							
public:
	// set initial target vertex coordinate as 0, where D is the dimension for the coordinates (may include extra data)
	EdgeInfo(VertexID vv1=0, VertexID vv2=0, EdgeID noVal=0) : v1(vv1), v2(vv2),quadError(0),no(noVal){}
	EdgeInfo(EdgeInfo *e) : v1(e->v1), v2(e->v2), quadError(e->quadError),no(e->no){}
	~EdgeInfo(){}
	
	// Compare an input edge with this edge, and return true if they are the same
	bool compareEdge(EdgeInfo *e){
		assert( e->v1 < e->v2 );
		return (this->v1 == e->v1 && this->v2 == e->v2);
	}
	// Swap the two endpoints if v1 > v2
	void swap_edge_endpoints(){
		if (v1 > v2){
			VertexID temp = v2; v2 = v1; v1 = temp;
		}
	}
};

struct HalfEdge
{
	VertexID v1, v2;		// two endpoints
	bool inside_face_flag;	// flag to denote if this half edge is already contained in a face
	bool visited_flag;		// flag for visited or not
};

//class p_switches : public MxHeapable
//{
//public:
//	VertexID from_id, to_id;
//	VertexID pixel_id;
//
//public:
//	p_switches(VertexID pixelid, VertexID fromid, VertexID toid):
//		pixel_id(pixelid), from_id(fromid), to_id(toid){ }
//};

class p_switches
{
public:
	VertexID from_id, to_id;
	VertexID vertex_id;

public:
	p_switches(VertexID vertexid, VertexID fromid, VertexID toid):
		vertex_id(vertexid), from_id(fromid), to_id(toid){ }
};


// Vertex quadrics information for each vertex
class VertexQuadrics
{
public:
	VertexQuadrics(VertexID id, uint dimVal);
	~VertexQuadrics();

	// Set cluster color
	inline void set_cluster_color(float f[3])
	{clusterColor[0] = f[0];clusterColor[1] = f[1];clusterColor[2] = f[2];}
	inline void set_cluster_color_flag(bool flag){cluster_color_flag = flag;}
	inline bool is_cluster_color_set(){return cluster_color_flag;}

	/////////////////////////////////////////////
	// Add an edge for this vertex
	void addEdge(EdgeID id);
	// Overload function: add an edge
	void addEdge(EdgeInfo *e);
	// Delete edge through its edge id
	void delete_edge(EdgeID id);
	// Get Edge number
	int getEdgeNumber() const {return edgeNumber;}
	// Clear edge list for this vertex
	inline void clear_all_edges(){edges->clear();edgeNumber=0;}
	// Clear edge list for this vertex
	inline void clear_all_faces(){faces->clear();faceNumber=0;}

	/////////////////////////////////////////////
	// Add a face by adding its faceID
	void addFace(FaceInfo *info_f);
	// Overload function: add a face by adding its faceID
	void addFace(FaceID id);
	// Delete face by deleting its faceID
	void delete_face(FaceID id);

	// Add a neighbor
	void add_Neighbor(VertexID id){neighbors->insert(id);}
	// Clear this vertex quadrics element.
	// Note: this function does not change the initial quadrics.
	void clearAll();
	// Get cluster size
	inline uint get_cluster_size(){return cluster->size();}
	// Copy vertex quadrics
	// void copy_one_quadrics(Quad &Q){*vtxQuadrics = Q;}
	
	// Determine if an edge is already included
	bool is_edge_included(EdgeID id);
	// Merge new cluster into current cluster
	void merge_cluster(VertexQuadrics *vtxQuad);
	// Delete one element from cluster
	void delete_vertex_from_cluster(VertexID id){cluster->erase(id);}
	// Add one element to the cluster
	void add_vertex_to_cluster(VertexID id){cluster->insert(id);}


public:
	vector<EdgeID> *edges;	// the edge indexes for this vertex
	vector<FaceID> *faces;	// the face indexes for this vertex
	//Quad* vtxQuadrics;		// current vertex quadrics
	//Quad* initialQuadrics;	// initial vertex quadrics (combination of all quadrics of adjacent faces of this vertex)
	MxQuadric* vtxQuadrics;
	MxQuadric* initialQuadrics;
	//MxQuadric* isotropicQuadrics;	// isotropic quadrics, i.e. alpha*I
	VertexID vertexId;		// vertex index in the model
	uint edgeNumber;		// edge number for this vertex
	uint faceNumber;		// face number for this vertex
	VectorType target;		// target vertex coordinate after contraction
	set<VertexID> *cluster;	// cluster for saving all vertex indexes of this cluster
	set<VertexID> *neighbors;		// neighbors in the ORIGINAL mesh
	bool cluster_color_flag;// cluster color no
	float clusterColor[4];	// cluster color
	float vertexQuality;	// final vertex quality value
};

inline bool compareEdgeInfo(EdgeInfo *e1, EdgeInfo *e2)
{
	return (e1->quadError < e2->quadError);
}


#endif