/************************************************************************/
/* 
*	acvt_mesh.h
*		Define acvt class based on a mesh model
*
*/
/************************************************************************/


#ifndef ACVTMESH_H
#define ACVTMESH_H
#include "PLYReader.h"
#include "VertexQuadrics.h"
#include "tools.h"
#include <set>
#include <algorithm>
#include "gfx.h"
#include "MxHeap.h"
#include "CrestLine.h"
#include "global.h"
#include <map>
#include <ctime>
#include "DualQEM.h"
#include <stack>

typedef map<pair<VertexID,VertexID>, uint> EdgePairType;

class ACVT
{
public:
	int dimVal;								// dimension value for each quadrics matrix
	//PCluster *pcluster;					// clusters
	uint valid_verts;						// valid number of VERTICES
	uint valid_faces;						// valid number of faces
	vector<VertexQuadrics*> *vtxQuadricsList;// quadrics list for all VERTICES
	PLYReader *mesh;						// mesh
	uint vertex_number;						// original vertex number
	uint face_number;						// original face number
	vector<EdgeInfo*> *edgeList;			// edge list
	vector<FaceInfo*> *faceList;			// updated face list
	vector<FaceInfo*> *faceList_new;		// result face list (obtained based on the result clusters)
	MxHeap heap;							// A maximum heap for operating all quadrics errors
	vector<float*> *normals_original;		// normals for VERTICES
	vector<float*> *normals_current;		// current normals for vertices
	vector<double> *vtxQuadricsWeights;		// vertex quadrics weights

	EdgePairType edgePairs;			// original edge pairs (first pair) and its number of frequency (the second paramter) in faces
	EdgePairType edgePairs_new;		// new edge pairs (first pair) and its number of frequency (the second paramter) in faces
	EdgePairType borderEdgePairs;	// border edges (first pair) and their faceIds (the second paramter)
	
	vector<VertexID> *cluster_belong;			// the cluster no. each vertex belongs to
	vector<set<VertexID>> *neighbor_clusters;// neighbor clusters

	vector<vector<HalfEdge*>*> *halfEdges;	// half edges for computing output trianglar faces

	//uint count_halfedge;
	//EdgePairType edge_visited;

	// Crest Line Variables
	CrestLine *ridge;						// ridges
	CrestLine *ravine;						// ravines

	// Constant parameters
	double etaVal;							// eta
	double gammaVal;						// gamma
	double vtxWeightThreshold;				// vertex weight threshold

	// Flags
	bool use_crest_line_flag;				// flag for using crest lines or not
	bool use_vertex_quality_flag;			// flag for using vertex quality weight or not
	
public:
	ACVT();
	~ACVT();

public:

	/*********************************************************************************************************/
	/*
	* Commmon Functions (Small functions)
	*/
	// Set eta and gamma values for the crest line weight and vertex quality weight, respectively
	inline void set_eta(double val){etaVal = val;}
	inline void set_gamma(double val){gammaVal = val;}
	inline void set_vertex_weight_threshold(double val){vtxWeightThreshold = val;}

	// Get dimensional value for each quadrics matrix
	inline uint dim() const { return dimVal; }

	// Get the quadric number
	uint quadrics_count() const {return vtxQuadricsList->size();}

	// Get one quadrics
	//Quad& Quadrics(uint i){return *(vtxQuadricsList->at(i)->vtxQuadrics);}
	MxQuadric& Quadrics(uint i){return *(vtxQuadricsList->at(i)->vtxQuadrics);}

	// Get one target vertex
	VectorType& TargetVertex(uint i){return vtxQuadricsList->at(i)->target;}

	// Get Isotropic Quadrics
	//MxQuadric& IsotropicQuadrics(uint i){return *(vtxQuadricsList->at(i)->isotropicQuadrics);}

	// Get one initial quadrics
	//Quad& InitialQuadrics(uint i){return *(vtxQuadricsList->at(i)->initialQuadrics);}
	MxQuadric& InitialQuadrics(uint i){return *(vtxQuadricsList->at(i)->initialQuadrics);}

	// Erase one vertex by initializing all its information
	void erase_vertex(VertexID i){vtxQuadricsList->at(i)->clearAll();}

	// Get original normals and current normals
	inline float* getOriginalNormal(VertexID i){return normals_original->at(i);}
	inline float* getCurrentNormal(VertexID i){return normals_current->at(i);}

	// Set flags
	inline void setCrestLineFlag(bool flag){use_crest_line_flag = flag;}
	inline void setVertexQualityFlag(bool flag){use_vertex_quality_flag = flag;}

	inline VertexID getClusterBelong(VertexID i){return cluster_belong->at(i);}

	// Get current vertex number and face number
	uint get_current_vertex_number();
	//uint get_current_face_number();

	/*********************************************************************************************************/
	/*
	* Large-Step Functions (Always connected to the buttons in the UI directly)
	*/
	// Read mesh model
	void readMesh(char *infilename);

	// Initialize the qslim method
	void init_qslim();

	// Compute initial quadrics for all VERTICES
	void collect_quadrics();	

	// Compute quadrics for one face with input face index id
	// Note: here "indexInFace" is the vertex index (0, 1, or 2) in this face for computing the face quadrics 
	void compute_face_quadrics(FaceID id, VertexID indexInFace, MxQuadric &Q);

	// Simplify
	void simplify();

	// Save the result VERTICES and faces
	// flag: 1 denotes saving mesh using QEM connectivity
	//		 2 denotes saving mesh using cluster connectivity
	void save_result_model(uint flag);

	// Optimization Step in one iteration
	// swap all qualified cluster border vertices in all clusters to reduce the total energy
	bool swap_once();

	// QEM optimization
	//void optimization();

	// Assign Cluster colors
	void assign_cluster_colors();

	// Add isotropic constraint, i.e. a diagonal weight for each quadrics to reduce the 
	// "inversion" circumstances between the border points of neighbor clusters in the result clusters
	void add_isotropic_constraint();

	// Update border constraint
	void update_border_constraint();
	
	// Compute the vertex weight
	void compute_vertex_weight();

	// Compute the final vertex qualities for each vertex
	void compute_final_vertex_quality();

	/*********************************************************************************************************/
	/*
	* Quadrics-related Functions
	*/
	// Assign additional constraint weight to each vertex quadrics
	void assign_weight_to_quadrics();

	// Constrain boundary by adding a large weight to the boundary edge VERTICES
	void constrain_boundary();	

	// Collect all original edge pairs and put the data in the "edgePairs" element
	// Meanwhile, save all original faces into "faceList"
	void compute_original_edge_pairs();

	// Compute the vertex extra data (usually vertex quality value) as a weight for quadrics
	void compute_vertex_quality_weight();

	// Compute the crest-line weights for each vertex
	// Note: the crest-line weight for each vertex is in inverse proportion to the nearest 
	// distance from the vertex to the crest lines
	void constain_crest_line_features();	
	
	// Apply contraction on a edge by updating the new quadrics Q and the new contraction
	// vertex coordinates of the first endpoint vertex of this edge
	// Note: for edge info=(v1,v2), I just contract the vertex v2 to v1 by deleting v2 and keeping v1,
	// then update the quadrics of v1 by adding v2's quadrics to it.
	void apply_contraction(EdgeInfo *info);
	
	// Given an input vertex index, get this vertex with extra data 
	void getVertexWithExtraData(uint id, VectorType& vec);
	
	// Get optimal edge: edge with smallest error
	//EdgeInfo * getOptimalEdge();

	// Update quadrics of all edges
	// void update_quadrics_errors(VertexID vid);
	
	// Compute total energy (summation of all quadric errors)
	double compute_total_energy();

	// Compute target vertex
	void compute_target(VertexID i);
	
	/*********************************************************************************************************/
	/*
	* Edge and Face Functions
	*/
	// Update edge information, including
	//   -- update edge endpoints of related VERTICES
	//   -- merge clusters
	//   -- update the quadrics of related edges
	void update_edge(EdgeInfo *info);

	// Update faces related to an input edge (Used to get connectivity inherited from original connectivity in QEM)
	void update_face(EdgeInfo *info);

	// Update the heap
	void update_edge_in_heap(EdgeInfo *e);

	// Determine if the edge (ev1, ev2) is included in the adjacent list of vertex vid
	bool is_edge_included(VertexID vid, VertexID ev1, VertexID ev2);

	// Compute initial edge information, including:
	// -- Compute the contraction quadrics errors of all edges 
	// -- Insert all edges into a maximum heap for QEM algorithm, with each edge's key the contraction quadrics error
	// -- Get all the neighbor faces of each vertex and save them
	void collect_edges();

	// Compute edge information, including the vertex after contraction for each edge
	void compute_edge_info(EdgeInfo *info);	

	// Delete an edge from the adjacent list of a vertex
	void delete_edge(VertexID vid, EdgeInfo *e);

	// Compute original vertex normals and initialize the current normals
	void initialize_vertex_normals();

	// Compute current normals for use of current mesh
	void compute_current_normals();

	// Compute normals for the QEM connectivity (not cluster-based connectivity)
	void compute_QEM_normals();

	/*
	* Compute new faces based on the connectivity of clusters using half-edge method.
	* This function contains some "correction" modifications against non-manifold clusters
	* to some extent. Because of this reason, however, this function is very slow.
	*
	* Note: this "half-edge" method may obtain all current normals inversely. If so,
	* you can use function "reverse_current_normals()" to reverse all normals.
	*/
	void compute_new_faces();

	/*
	* Compute new faces based on the connectivity of clusters using half-edge method
	* under the condition of all manifold clusters.
	* 
	* Note: 
	* (1) The difference between this function and the above function "compute_new_faces()" 
	* is that this function runs on all manifold clusters, while the latter function runs on 
	* arbitrary clusters.
	* (2) This method may obtain all current normals inversely. If so, you can use function
	* "reverse_current_normals()" to reverse all normals.
	*/
	void compute_new_faces_nonrecursive();

	/* 
	* Compute the difference energy between initial and result by removing vertex vid from cluster i to cluster j.
	* If the result energy is smaller than the initial energy, then return a positive delta value; 
	* otherwise return a negative value.
	*/
	double compute_delta_energy(VertexID vid, VertexID i, VertexID j);

	// Reverse current normals
	void reverse_current_normals();

	// Create new edgelist and heap for a new QEM process
	void compute_new_edgelist_heap();

	// Compute new edge pairs and save them into variable "edgePairs"
	// Note: You must run function "compute_all_neighbor_clusters()" before running this function.
	void compute_new_edge_pairs();

	/*********************************************************************************************************/
	/*
	* Half Edge Functions (Used in a recursion method to create triangles from edges)
	*/
	// Run the "half-edge" method to get triangles recursively starting from one half-edge.
	// This function will save the newly created faces in the "faceList_new"
	void run_half_edge_method(HalfEdge *halfedge_start, vector<vector<HalfEdge*>*>* halfEdges_all);

	// Non-recursive version of the "half-edge" method to get triangles starting from one half-edge.
	// This function will save the newly created faces in the "faceList_new"
	void run_half_edge_method_nonrecursive(HalfEdge *halfedge_start);

	// Get a Half Edge starting from v1 and ends in v2
	HalfEdge* get_half_edge(VertexID v1, VertexID v2, vector<vector<HalfEdge*>*>* halfedges_all);

	// Find the half edge of an input half edge
	HalfEdge* find_opposite_halfedge(HalfEdge *halfedge_input, vector<vector<HalfEdge*>*>* halfedges_all);

	// Check if an input half edge can form a new face if combined with the other two half edges in an input
	// half edge list. If so, save the other two half edges return true; otherwise, save "NULL" in the other two
	// half edges and return false
	bool check_halfedge_form_new_face(HalfEdge* halfedge_1, vector<HalfEdge*> *halfedge_list, VertexID& halfedge2_index, VertexID& halfedge3_index);

	/*********************************************************************************************************/
	/*
	* Cluster Operation Functions
	*/
	// Split clusters into separate connected components, including the following steps:
	// -- (1) Check the connectivity for each cluster to find such clusters whose vertices do not compose a connected component 
	// -- (2) Split such clusters into different sub-clusters, with each connected component a sub-cluster
	// -- (3) Check if one cluster is totally inside another "mother" cluster (all the neighbors of border vertices of this cluster belongs to a same cluster). 
	//        If so, merge it into its "mother" cluster
	void split_clusters();

	// Split one cluster into two sub-clusters using the first PCA principal axis
	void split_cluster_into_two(VertexID cluster_id, vector<vector<VertexID>> &subclusters);

	// Check non manifold clusters
	// Note: must compute all neighbor clusters before running this function
	void check_nonmanifold_clusters(set<VertexID>& cluster_ids);

	// Split non manifold clusters into different sub-clusters. If there are still non manifold clusters to be split, 
	// return true; otherwise return false.
	bool split_nonmanifold_clusters();

	// Check if a cluster is totally inside another cluster
	bool check_cluster_inside_another(vector<VertexID> &vecTemp);

	// Merge a cluster to a new cluster
	void merge_cluster(VertexID from_cluster_id, VertexID to_cluster_id);

	// Merge all those small clusters totally inside some other clusters
	void merge_all_clusters_inside_another();

	// Get the empty cluster id list
	void get_empty_cluster_list(vector<VertexID> &cluster_list);

	// Create new clusters and put them in the empty list.
	// Input: 
	//	empty_cluster_list: empty cluster list (can be obtained by the "get_empty_cluster_list" function);
	//	subclusters: new clusters to be created
	void create_new_clusters(vector<VertexID> &empty_cluster_list, vector<vector<VertexID>> &newclusters);

	// Save cluster informations and output this file
	// Note: The output file "*.cluster" is composed of N lines, where N is the number
	// of vertices in the original model. In the ith line, the first value is the number of 
	// of vertices for the ith cluster, and the second to the last values of the line denote the vertex ids
	// which the ith cluster contains. If there is no vertex for the ith cluster, then there is only 
	// a 0 in this ith line.
	void save_cluster_in_file(char *filename);

	// Open Cluster file "*.cluster"
	// Note: You must open the mesh first before opening the "*.cluster" file
	void open_cluster_file(char *filename, vector<vector<VertexID>*> *cluster_all);

	// Set all clusters based on an input "*.cluster" file
	void set_clusters_from_input(vector<vector<VertexID>*> *cluster_all);

	// Compute clusters each vertex belongs to
	void compute_cluster_belong();

	// Get total number of vertices in clusters
	uint get_vertex_number_in_clusters();

	// Update the belonging cluster of a vertex vid to cid
	inline void update_cluster_belong(VertexID vid, VertexID cid){cluster_belong->at(vid) = cid;}

	// Determine if a vertex is a cluster border vertex
	bool is_cluster_border(VertexID vid);

	// Compute all neighbor clusters for each cluster and save the data in "neighbor_clusters"
	void compute_all_neighbor_clusters();

	// Compute an input cluster's neighbor clusters
	void compute_neighbor_cluster(VertexID i,set<VertexID>& neighbors);
};


#endif // !ACVTMESH_H
