#include "acvt_mesh.h"
#include <fstream>


ACVT::ACVT()
{
	mesh = NULL;
	dimVal = 3;		// initial dimension for quadrics is 3
	valid_verts = valid_faces = 0;
	vertex_number = face_number = 0;
	vtxQuadricsList = new vector<VertexQuadrics*>;
	edgeList = new vector<EdgeInfo*>;
	faceList = new vector<FaceInfo*>;
	faceList_new = new vector<FaceInfo*>;
	cluster_belong = new vector<uint>;
	neighbor_clusters = new vector<set<VertexID>>;
	normals_original = new vector<float*>;
	normals_current = new vector<float*>;
	vtxQuadricsWeights = new vector<double>;
	use_vertex_quality_flag = use_crest_line_flag = false;
	halfEdges = new vector<vector<HalfEdge*>*>;
	ridge = NULL;
	ravine = NULL;
	etaVal = 3;
	gammaVal = 5;
	vtxWeightThreshold = 0.1;
}

ACVT::~ACVT()
{
	if (!vtxQuadricsList->empty())
	{
		for (vector<VertexQuadrics*>::iterator iter = vtxQuadricsList->begin();
			iter != vtxQuadricsList->end(); ++iter)
		{
			delete *iter;
		}
	}
	vtxQuadricsList->clear();
	ClearVector(*vtxQuadricsList);

	if (!edgeList->empty())
	{
		for (vector<EdgeInfo*>::iterator iter = edgeList->begin();
			iter != edgeList->end(); ++iter)
		{
			delete *iter;
		}
	}
	edgeList->clear();
	ClearVector(*edgeList);

	if (!faceList->empty())
	{
		for (vector<FaceInfo*>::iterator iter = faceList->begin();
			iter != faceList->end(); ++iter)
		{
			delete *iter;
		}
	}
	faceList->clear();
	ClearVector(*faceList);

	if (!faceList_new->empty())
	{
		for (vector<FaceInfo*>::iterator iter = faceList_new->begin();
			iter != faceList_new->end(); ++iter)
		{
			delete *iter;
		}
	}
	faceList_new->clear();
	ClearVector(*faceList_new);

	if (normals_original != NULL)
	{
		for (uint i = 0; i != normals_original->size(); ++i)
		{
			delete [](normals_original->at(i));
		}
	}
	normals_original->clear();
	ClearVector(*normals_original);

	if (normals_current != NULL)
	{
		for (uint i = 0; i != normals_current->size(); ++i)
		{
			delete [](normals_current->at(i));
		}
	}
	normals_current->clear();
	ClearVector(*normals_current);

	cluster_belong->clear();
	ClearVector(*cluster_belong);
	vtxQuadricsWeights->clear();
	ClearVector(*vtxQuadricsWeights);

	// Clear relative data structure
	if (!neighbor_clusters->empty())
	{
		for (uint i = 0; i != vertex_number; ++i)
		{
			neighbor_clusters->at(i).clear();
		}
	}
	neighbor_clusters->clear();
	ClearVector(*neighbor_clusters);

	edgePairs.clear();
	edgePairs_new.clear();
	borderEdgePairs.clear();

	// Clear half edges
	for (uint i = 0; i != halfEdges->size(); ++i)
	{
		for (uint j = 0; j != halfEdges->at(i)->size(); ++j)
		{
			delete halfEdges->at(i)->at(j);
		}
		halfEdges->at(i)->clear();
		ClearVector(*(halfEdges->at(i)));
	}
	halfEdges->clear();
	ClearVector(*halfEdges);

#ifdef CRESTLINE
	delete ridge;
	delete ravine;
#endif
}

void ACVT::readMesh(char* infilename)
{
	mesh = new PLYReader(infilename, 1.0);
	valid_verts = mesh->getNumVertices();
	valid_faces = mesh->getNumTriangles();
	vertex_number = valid_verts;
	face_number = valid_faces;
	vtxQuadricsWeights->resize(vertex_number, 1.0);

	// Extra dimension for each vertex
	//if (mesh->extraTypeNumber != 0)
	//{
	//	dimVal += mesh->extraTypeNumber;
	//}

#ifdef CRESTLINE
	char *ridgename = new char[strlen(mesh->modelname)+10];
	strncpy(ridgename, mesh->modelname, strlen(mesh->modelname)-3);
	ridgename[strlen(mesh->modelname)-3] = '\0';
	strcat(ridgename, "ridges");
	//ridgename[strlen(mesh->modelname)+5] = '\0';
	char *ravinename = new char[strlen(mesh->modelname)+8];
	strncpy(ravinename, mesh->modelname, strlen(mesh->modelname)-3);
	ravinename[strlen(mesh->modelname)-3]='\0';
	strcat(ravinename, "ravines");

	ridge = new CrestLine;
	if (ridge->readCrestLineFile(ridgename))
	{
		ridge->update_crestlines();
	}
	ravine = new CrestLine;
	if (ravine->readCrestLineFile(ravinename))
	{
		ravine->update_crestlines();
	}
	delete []ridgename;
	delete []ravinename;
#endif

	// Note: The step to compute normals must be put behind the read-crest-line-file step,
	// since it contains computing normals for crest lines.
	initialize_vertex_normals();

	// Collect all edge pairs and their frequencies
	// Note: this step must be put before the following constrain steps
	compute_original_edge_pairs();
}

void ACVT::init_qslim()
{
	for (uint i = 0; i != valid_verts; ++i)
	{
		//float *v = mesh->getVertex(i);
		VertexQuadrics *vtxQuad = new VertexQuadrics(i, dim());
		getVertexWithExtraData(i, vtxQuad->target);
		vtxQuadricsList->push_back(vtxQuad);
	}

	// Firstly, Compute all vertex quadrics
	collect_quadrics();

	// Add isotropic constraint at last (must be the last step after computing initial quadrics step and 
	// adding weights to all quadrics step)
	add_isotropic_constraint();

	// Step 4: assign additional weight for each quadrics
	// Note: this step must be put behind the constrain steps which computes weights for quadrics
	assign_weight_to_quadrics();

	// Next, constrain the boundary by adding a quadrics to each pair of edge endpoints
	//constrain_boundary();


	// Final step: create all edges and insert the edge constrain error into the maximum heap
	collect_edges();

	// Make a copy of initial quadrics
	for(uint i=0; i != valid_verts; i++)
	{
		InitialQuadrics(i) = Quadrics(i);
	}
}

void ACVT::collect_quadrics()
{
	for (uint i = 0; i != mesh->getNumTriangles(); ++i)
	{
		MxQuadric Q(dim());
		int *t = mesh->getTriangle(i);
		for (uint j = 0; j != 3; ++j)
		{
			// The following line incorrectly computes quadrics using the same vertex 
			// for all the three vertices in the same face (also the same as Yiqi's code).
			//compute_face_quadrics(i, 0, Q);
			
			// This is the correct method to compute quadrics for all the three vertices
			// in the same face
			compute_face_quadrics(i, j, Q);
			Quadrics(t[j]) += Q;
		}
	}
}

void ACVT::compute_face_quadrics(FaceID id, VertexID indexInFace, MxQuadric &Q)
{
	VectorType v0(dim());
	VectorType v1(dim());
	VectorType v2(dim());

	int *t = mesh->getTriangle(id);
	getVertexWithExtraData(t[indexInFace], v0);
	getVertexWithExtraData(t[(indexInFace+1)%3], v1);
	getVertexWithExtraData(t[(indexInFace+2)%3], v2);

	// The forth parameter r is a weight constraint to the face
	double r = 1.0;
	MxQuadric Q_p(v0, v1, v2, r);
	Q = Q_p;

	// Add isotropic constraint for each face quadrics (this step is integrated into
	// a new function as one step in the "init_qslim" function)
	//MxQuadric Q_iso(v0, isotropic_coefficient);
	//Q += Q_iso;

	//IsotropicQuadrics(t[indexInFace]) += Q_iso;
}

/************************************************************************/
/* 
* Edge Functions
*/
void ACVT::compute_original_edge_pairs()
{
	// Collect all initial edges of input mesh
	// Use std::set to get unique edge
	for (uint i = 0; i != valid_faces; ++i)
	{
		int *t = mesh->getTriangle(i);

		// Save faces
		FaceInfo *info_f = new FaceInfo(t,i);
		faceList->push_back(info_f);

		// Get each unique edge
		for (uint j = 0; j != 3; ++j)
		{
			int a = t[j], b = t[(j+1)%3];
			assert(a != b);
			if (a < b)
			{
				edgePairs[pair<uint,uint>(a,b)]++;
			}
			else
			{
				edgePairs[pair<uint,uint>(b,a)]++;
			}
		}
	}

	// Determine all the face ids of all border edges
	for (uint i = 0; i != valid_faces; ++i)
	{
		int *t = mesh->getTriangle(i);
		// Get each unique edge
		for (uint j = 0; j != 3; ++j)
		{
			int a = t[j], b = t[(j+1)%3];
			sort2(a,b);
			if (edgePairs[pair<uint,uint>(a,b)] == 1)
			{
				borderEdgePairs[pair<uint,uint>(a,b)] = i;
			}
		}
	}
}

void ACVT::collect_edges()
{
	// Create all edges and compute contraction quadrics information
	uint count = 0;
	for (map<pair<uint,uint>, uint>::iterator iter = edgePairs.begin(); 
		iter != edgePairs.end(); ++iter)
	{
		// Create this edge
		EdgeInfo *e = new EdgeInfo((*iter).first.first, (*iter).first.second, count);
		count++;
		
		// Compute the edge information including the quadrics error for this edge
		compute_edge_info(e);
		
		// Put the edge key inside the maximum heap as the negative quadrics error
		e->heap_key(e->quadError);
		
		// Save the edge
		edgeList->push_back(e);
		
		// Push the edge into the heap
		heap.insert(e);

		// Save this edge into the edgelist of related vertex
		vtxQuadricsList->at(e->v1)->addEdge(e);
		vtxQuadricsList->at(e->v2)->addEdge(e);
		vtxQuadricsList->at(e->v1)->add_Neighbor(e->v2);
		vtxQuadricsList->at(e->v2)->add_Neighbor(e->v1);
	}

	// Set adjacent faces for each vertex
	for (uint i = 0; i != valid_faces; ++i)
	{
		FaceInfo *info_f = faceList->at(i);
		vtxQuadricsList->at(info_f->v[0])->addFace(info_f);
		vtxQuadricsList->at(info_f->v[1])->addFace(info_f);
		vtxQuadricsList->at(info_f->v[2])->addFace(info_f);
	}
}

void ACVT::compute_edge_info(EdgeInfo *info)
{
	VertexID i = info->v1, j = info->v2;
	//const Quad Qi=Quadrics(i), Qj=Quadrics(j);
	// For the new contract vertex of this edge, Q = Qi + Qj
	MxQuadric Q = Quadrics(i);
	Q += Quadrics(j);

	double err;
	VectorType vtxCoord(dimVal);
	if( Q.optimize(vtxCoord) )
	{
		// If the Quadrics matrix is invertible, then just compute the target vertex
		// coordinate vtxCoord =A^(-1)*b and compute the error as the Q^(-1) * v
		err = Q(vtxCoord);
	}
	else
	{
		// If the Quadrics matrix is not invertible, then get the smaller error on these
		// two vertices
		VectorType v0(dim()), v1(dim());
		//getVertexWithExtraData(i, v0);
		//getVertexWithExtraData(j, v1);
		v0 = TargetVertex(i);
		v1 = TargetVertex(j);
		double err_0 = Q(v0), err_1 = Q(v1);
		err = (err_0 < err_1) ? err_0 : err_1;
		//cout<<"Special Edge ("<<i<<","<<j<<"):"<<err<<endl;
	}
	// The punishment of negative error, make it a great positive error 
	//if (err < 0) 
	//	err += 1e200;
	
	//err += 1e-10 * Q.area(); // Make sure that area works even the error is zero
	//     if( weight_by_area )
	// 	err / Q.area();

	// Add the error
	// Note: The QEM needs to get the smallest error each time, but since we use a maximum heap
	// here for this pop out operation, then we have to put the negative error into the heap 
	// and pop out the largest one (corresponding to the smallest error)
	info->quadError = -err;
}

/************************************************************************/


void ACVT::getVertexWithExtraData(uint id, VectorType& vec)
{
	float *v = mesh->getVertex(id);
	for (uint i = 0; i != dimVal; ++i)
	{
		vec[i] = v[i];
	}
	/*
	* Note: Using vertex quality information during QEM is very likely to 
	* obtain bad result, since this project has already modified QEM by adding isotropic
	* constraint in the initialization of computing quadrics of each vertex. Therefore,
	* here we comment the following lines and put all the extra data dimension values as 0.
	*/
	
	// Extra Dimensional with all 0 values
	//if (mesh->extraTypeNumber != 0)
	//{
	//	for (uint i = 0; i != mesh->extraTypeNumber; ++i)
	//	{
	//		vec[i+3] = 0;
	//	}
	//}

	// Extra Dimensional with exact values from input mesh
	//if (mesh->extraTypeNumber != 0)
	//{
	//	for (uint i = 0; i != mesh->extraTypeNumber; ++i)
	//	{
	//		vec[i+3] = mesh->getVtxExtraData(id, i);
	//	}
	//}
}

void ACVT::simplify()
{
	EdgeInfo *info = (EdgeInfo *)heap.extract();
	uint v1 = info->v1, v2 = info->v2;
	//cout<<"Edge ("<<v1<<","<<v2<<"):"<<info->heap_key()<<endl;

	if (!vtxQuadricsList->at(v1)->cluster->empty() && !vtxQuadricsList->at(v2)->cluster->empty())
	{
		// If the cluster belonging to each of the two VERTICES of this edge is not empty,
		// then apply contraction by merging v2 into v1
		apply_contraction(info);

		// Update all edges related to this edge
		update_edge(info);

		// Update faces related to this edge
		update_face(info);

		// Erase the vertex v2
		erase_vertex(v2);

		valid_verts--;
	}
}

void ACVT::apply_contraction(EdgeInfo *info)
{
	MxQuadric Q = Quadrics(info->v1), Q2 = Quadrics(info->v2);
	Q += Q2;
	Quadrics(info->v1) = Q;

	// Update the new target vertex
	VectorType v(dim());
	bool flag = Q.optimize(v);
	if (flag)
	{
		TargetVertex(info->v1) = v;
	}
	else
	{
		getVertexWithExtraData(info->v1, TargetVertex(info->v1));
	}
}

void ACVT::update_edge(EdgeInfo *info)
{
	VertexID v1 = info->v1, v2 = info->v2;
	VertexQuadrics *v1VtxQuad = vtxQuadricsList->at(v1);
	VertexQuadrics *v2VtxQuad = vtxQuadricsList->at(v2);
	
	// Firstly, remove this edge from heap 
	if (info->is_in_heap())
	{
		heap.remove(info);
	}

	// Delete this edge from the edge list of the two endpoints
	v1VtxQuad->delete_edge(info->no);
	v2VtxQuad->delete_edge(info->no);	
	//delete_edge(v1, info);
	//delete_edge(v2, info);

	// Update all quadrics of all adjacent edges of vertex v1
	for (uint i = 0; i != v1VtxQuad->edgeNumber; ++i)
	{
		EdgeID id = v1VtxQuad->edges->at(i);
		EdgeInfo *e = edgeList->at(id);
		update_edge_in_heap(e);
	}

	// Update the endpoints of all adjacent edges of vertex v2, since
	// all the indexes of these edges are changed 
	for (uint i = 0; i != v2VtxQuad->edgeNumber; ++i)
	{
		EdgeID id = v2VtxQuad->edges->at(i);
		EdgeInfo *e = edgeList->at(id);

		// New edge that e will be converted to
		VertexID a = (e->v1 == v2) ? e->v2 : e->v1;
		VertexID b = v1;
		sort2(a,b);

		// the endpoint different from v2
		VertexID c = (e->v1 == v2) ? e->v2 : e->v1;		
	
		if (a != b)
		{
			// Update the new edge e in the adjacent list of vertex v1
			if (is_edge_included(v1, a, b))
			{
				// If this edge already exists in the adjacent list of v1
				if (e->is_in_heap())
				{
					heap.remove(e);
				}
				// Delete this edge from the adjacent list of the different endpoint from v2
				delete_edge(c, e);
			}
			else
			{
				// Update the endpoints of the new edge
				e->v1 = a;
				e->v2 = b;
				update_edge_in_heap(e);
				if (!v1VtxQuad->is_edge_included(e->no))
				{
					v1VtxQuad->addEdge(e);
				}
			}
		}
	}

	// Merge cluster from v2 to v1
	v1VtxQuad->merge_cluster(v2VtxQuad);
}

void ACVT::update_face(EdgeInfo *info)
{
	VertexID v1 = info->v1, v2 = info->v2;
	VertexQuadrics *v1VtxQuad = vtxQuadricsList->at(v1);
	VertexQuadrics *v2VtxQuad = vtxQuadricsList->at(v2);
	for (uint i = 0; i != v2VtxQuad->faces->size(); ++i)
	{
		FaceID id = v2VtxQuad->faces->at(i);
		FaceInfo *info_f = faceList->at(id);
		bool flag = true;
		// For each face from the adjacent list of v2
		for (uint j = 0; j != v1VtxQuad->faces->size(); ++j)
		{
			if (id == v1VtxQuad->faces->at(j))
			{
				// if this face is also adjacent to v1, then delete it
				// from adjacent list of v1
				v1VtxQuad->delete_face(id);

				// Update the adjacent list of the 3rd vertex of this face
				// the third vertex
				VertexID v3 = info_f->v[0] + info_f->v[1] + info_f->v[2] - v1 - v2;
				VertexQuadrics *v3VtxQuad = vtxQuadricsList->at(v3);
				v3VtxQuad->delete_face(id);
				
				// Delete this face from the face list by setting its flag as false
				info_f->delete_face();

				flag = false;
				break;
			}
		}
		if (flag)
		{
			// if this face is not adjacent to v1, then update its face VERTICES
			for (uint i = 0; i != 3; ++i)
			{
				if (info_f->v[i] == v2)
				{
					info_f->v[i] = v1;
					break;
				}
			}
			// Add the face to v1
			v1VtxQuad->addFace(id);
		}
	}
}

void ACVT::delete_edge(VertexID vid, EdgeInfo *e)
{
	VertexID v1 = e->v1, v2 = e->v2;
	VertexQuadrics *vtxQuad = vtxQuadricsList->at(vid);
	if (vtxQuad->edgeNumber != 0)
	{
		uint i = 0;
		for (i = 0; i != vtxQuad->edgeNumber; ++i)
		{
			EdgeInfo *info = edgeList->at(vtxQuad->edges->at(i));
			if ((info->v1 == v1) && (info->v2 == v2))
			{
				break;
			}
		}
		vtxQuad->edges->erase(vtxQuad->edges->begin()+i);
		vtxQuad->edgeNumber = vtxQuad->edges->size();
	}
}


bool ACVT::is_edge_included(VertexID vid, VertexID ev1, VertexID ev2)
{
	VertexQuadrics *vtxQuad = vtxQuadricsList->at(vid);
	if (vtxQuad->edgeNumber == 0)
	{
		return false;
	}
	else
	{
		for (uint i = 0; i != vtxQuad->edgeNumber; ++i)
		{
			EdgeInfo *e = edgeList->at(vtxQuad->edges->at(i));
			if ((e->v1 == ev1) && (e->v2 == ev2))
			{
				return true;
			}
		}
		return false;
	}
}


void ACVT::update_edge_in_heap(EdgeInfo *e)
{
	if (e->is_in_heap())
	{
		MxQuadric Q = Quadrics(e->v1);
		Q += Quadrics(e->v2);
		VectorType optVtx(dim());
		// Update the quadric error of this new edge
		double err = 0;
		if( Q.optimize(optVtx) )
		{// Q is invertible 
			err = Q(optVtx);
		}
		else
		{
			// If Q is not invertible, use one endpoint that results in
			// smaller energy as the target vertex
			VectorType v0(dim()), v1(dim());
			VertexID i = e->v1, j = e->v2;
			v0 = TargetVertex(i);
			v1 = TargetVertex(j);
			double err_0 = Q(v0), err_1 = Q(v1);
			err = (err_0 < err_1) ? err_0 : err_1;
		}
		e->quadError = -err;
		e->heap_key(e->quadError);
		heap.update(e);
	}
}

double ACVT::compute_total_energy()
{
	double energy = 0;
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vertQuad = vtxQuadricsList->at(i);
		if (!vertQuad->cluster->empty())
		{
			MxQuadric Q = Quadrics(i);
			energy += Q(TargetVertex(i));
			//cout<<energy<<endl;
		}
	}
	return energy;
}


// NOTE: this function still has bugs!!!!!!!!!
void ACVT::constrain_boundary()
{
	for (map<pair<uint,uint>, uint>::iterator iter = borderEdgePairs.begin();
		iter != borderEdgePairs.end(); ++iter)
	{
		// This edge only contains in on face, this means it is a border edge
		VertexID v0id = (*iter).first.first;
		VertexID v1id = (*iter).first.second;
		VectorType v0(dimVal), v1(dimVal);
		getVertexWithExtraData(v0id, v0);
		getVertexWithExtraData(v1id, v1);
		MxQuadric Q(v0, v1);

	//	// Border edge vector e
	//	VectorType e = v1 - v0;
	//	// Get face normal of this border edge
	//	VertexID fid = borderEdgePairs[pair<uint,uint>(v0id, v1id)];
	//	int *v = mesh->getTriangle(fid);
	//	float *v_0 = mesh->getVertex(v[0]);
	//	float *v_1 = mesh->getVertex(v[1]);
	//	float *v_2 = mesh->getVertex(v[2]);
	//	float e0[3], e1[3], n[3];
	//	for (uint j = 0; j != 3; ++j)
	//	{
	//		e0[j] = v_1[j] - v_0[j];
	//		e1[j] = v_2[j] - v_0[j];
	//	}

	//	/*************************************************************************
	//	* Note: here we have three methods to compute the quadrics for boundaries
	//	*/
	//	//////////////////////////////////////////////////////////
	//	// Method 1: Yiqi's method (a little weird)
	//	//e0[2] = e1[2] = 0;
	//	//crossProduct(e0, e1, n);
	//	//Quad Q(e, n, v0);	

	//	//////////////////////////////////////////////////////////
	//	// Method 2: adding a vertical plane
	//	crossProduct(e0, e1, n);
	//	float n2[3];
	//	e0[0] = e[0];
	//	e0[1] = e[1];
	//	e0[2] = e[2];
	//	crossProduct(e0, n, n2);	
	//	MxQuadric Q(n2, v0);

	//	//////////////////////////////////////////////////////////
	//	// Method 3: adding a vertical plane and an existing plane which contains the face
	//	// of the boundary edge (method in the paper Quadric-based simplification in any dimension)
	//	//Quad Q(v0, v1);

	//	Q *= penalty_factor;

	//	//cout<<"Border Edge ("<<v0id<<","<<v1id<<") in face "<<fid<<endl;
	//	//cout<<Q.A<<endl;
	//	//cout<<Q.b<<endl;
	//	//cout<<Q.c<<endl<<endl;


		Q *= penalty_factor;

		Quadrics(v0id) += Q;
		Quadrics(v1id) += Q;
	}
}

void ACVT::update_border_constraint()
{
	compute_cluster_belong();

	for (map<pair<uint,uint>, uint>::iterator iter = borderEdgePairs.begin();
		iter != borderEdgePairs.end(); ++iter)
	{
		// This edge only contains in on face, this means it is a border edge
		VertexID v0id = (*iter).first.first;
		VertexID v1id = (*iter).first.second;
		VectorType v0(dim()), v1(dim());
		getVertexWithExtraData(v0id, v0);
		getVertexWithExtraData(v1id, v1);

		/*
		* Note: here we have two methods to compute the quadrics for boundary edges
		*/
		// Method 1: adding a vertical plane and an existing plane which contains the face
		// of the boundary edge (method in the paper Quadric-based simplification in any dimension)
		MxQuadric Q(v0, v1);

		// Method 2: adding a vertical plane
		//VectorType e = v1 - v0;
		//VertexID fid = borderEdgePairs[pair<uint,uint>(v0id, v1id)];
		//int *v = mesh->getTriangle(fid);
		//float *v_0 = mesh->getVertex(v[0]);
		//float *v_1 = mesh->getVertex(v[1]);
		//float *v_2 = mesh->getVertex(v[2]);
		//float e0[3], e1[3], n[3];
		//for (uint j = 0; j != 3; ++j)
		//{
		//	e0[j] = v_1[j] - v_0[j];
		//	e1[j] = v_2[j] - v_0[j];
		//}
		//crossProduct(e0, e1, n);
		//float n2[3];
		//e0[0] = e[0];
		//e0[1] = e[1];
		//e0[2] = e[2];
		//crossProduct(e0, n, n2);	
		//Quad Q(n2, v0);

		Q *= penalty_factor;	
		// Note: add the border quadrics to the clusters of the two endpoints instead of 
		// the two endpoints themselves
		VertexID v0_clusterID = cluster_belong->at(v0id);
		VertexID v1_clusterID = cluster_belong->at(v1id);
		Quadrics(v0_clusterID) += Q;
		Quadrics(v1_clusterID) += Q;
		// Need to update the initial quadrics of border edge endpoints for future use
		InitialQuadrics(v0id) += Q;
		InitialQuadrics(v1id) += Q;
	}

	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			MxQuadric Q = Quadrics(i);
			Q.optimize(TargetVertex(i));
		}
	}
}

void ACVT::assign_weight_to_quadrics()
{
	for (uint i = 0; i != vertex_number; ++i)
	{
		//VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		Quadrics(i) *= vtxQuadricsWeights->at(i);
	}
}

// Note: this function still has bugs!!!!!!!!!!!!!!!!!!
void ACVT::add_isotropic_constraint()
{
	for (uint i = 0; i != mesh->getNumTriangles(); ++i)
	{
		int *t = mesh->getTriangle(i);
		for (uint j = 0; j != 3; ++j)
		{
			VectorType v(dim());
			getVertexWithExtraData(t[j], v);
			
			MxQuadric Q(v, isotropic_coefficient);
			Quadrics(t[j]) += Q;
		}
	}
}

void ACVT::compute_vertex_quality_weight()
{
	vector<double> weights;
	for (uint i = 0; i != vertex_number; ++i)
	{
		weights.push_back(mesh->getVtxExtraData(i));
	}
	double maxVal = *max_element(weights.begin(), weights.end());
	const double threVal = (1.0-vtxWeightThreshold)/2;
	for (uint i = 0; i != vertex_number; ++i)
	{
		double w = weights[i] / maxVal;
		if (w - vtxWeightThreshold < 1e-5)
		{
			weights[i] = pow(w + threVal, gammaVal);
		}
		else
		{
			weights[i] = gammaVal * 5;
		}
		vtxQuadricsWeights->at(i) *= weights[i];
	}
	weights.clear();
	ClearVector(weights);
}

void ACVT::compute_vertex_weight()
{
	vtxQuadricsWeights->clear();
	vtxQuadricsWeights->resize(vertex_number, 1.0);

	// Compute the crest-line-weight 
#ifdef CRESTLINE
	if (use_crest_line_flag)
	{
		constain_crest_line_features();
	}
#endif
	// Compute vertex quality weight
	if (use_vertex_quality_flag)
	{
		compute_vertex_quality_weight();
	}
}


void ACVT::constain_crest_line_features()
{
	struct kdtree *kd = kd_create(3);
	struct kdres *presults;
	float *v;
	double vtx[3];
	double dis = 0.0;
	int *data;
	// Add ridge vertices into the kd-tree
	for (uint i = 0; i != ridge->vertexNumber; ++i)
	{
		if (ridge->is_feature(i))
		{
			float *v = ridge->getVertex(i);
			vtx[0] = (double)v[0];
			vtx[1] = (double)v[1];
			vtx[2] = (double)v[2];
			data = new int;
			*data = i;
			kd_insert3(kd, vtx[0],vtx[1],vtx[2],data);
		}
	}
	// Add ravine vertices into the kd-tree
	for (uint i = 0; i != ravine->vertexNumber; ++i)
	{
		if (ravine->is_feature(i))
		{
			float *v = ravine->getVertex(i);
			vtx[0] = (double)v[0];
			vtx[1] = (double)v[1];
			vtx[2] = (double)v[2];
			data = new int;
			// Note: the id of all ravine vertices are different from ridge vertices
			*data = i+ridge->vertexNumber;
			kd_insert3(kd, vtx[0],vtx[1],vtx[2],data);
		}
	}

	// Get the nearest VERTICES of each vertex from the original model
	double pos[3];
	vector<double> weights;
	for (uint i = 0; i != vertex_number; ++i)
	{
		v = mesh->getVertex(i);
		vtx[0] = (double)v[0];
		vtx[1] = (double)v[1];
		vtx[2] = (double)v[2];
		// Find the nearest vertex of vtx
		presults = kd_nearest (kd, vtx);
		int *ptr_nearestVtx = (int*)kd_res_item( presults, pos );
		// Nearest distance from an arbitrary vertex to crest lines
		dis = sqrt((vtx[0]-pos[0])*(vtx[0]-pos[0])
			+ (vtx[1]-pos[1])*(vtx[1]-pos[1])
			+ (vtx[2]-pos[2])*(vtx[2]-pos[2]));
		weights.push_back(dis);
		kd_res_free( presults );
	}

	//delete kd-tree;
	kd_free( kd );


	// Compute final weight for the quadrics of each vertex
	double maxWeight = *max_element(weights.begin(), weights.end());
	for (uint i = 0; i != weights.size(); ++i)
	{
		weights[i] = pow(weights[i]/maxWeight + epsilon_const, etaVal);
	}

	for (uint i = 0; i != vertex_number; ++i)
	{
		// Crest line weight is in inverse proportion to the nearest distance from
		// each vertex to the crest lines
		vtxQuadricsWeights->at(i) /= weights[i];
	}
	weights.clear();
	ClearVector(weights);
}

void ACVT::save_result_model(uint flag)
{
	if (flag != 1 && flag != 2)
	{
		cout<<"You must save mesh using QEM connectivity or cluster connectivity!"<<endl;
		return;
	}
	mesh->clearNewVertexFaceVector();
	mesh->newFaces = new vector<int*>;
	mesh->newVertices = new vector<float*>;

	vector<int*> *face_result = mesh->newFaces;
	vector<float*> *vertices_result = mesh->newVertices;
	vector<vector<float>> *vtx_extra_result;
	if (mesh->extraTypeNumber != 0)
	{
		mesh->vtxNewExtraData = new vector<vector<float>>;
		vtx_extra_result = mesh->vtxNewExtraData;
		vtx_extra_result->resize(1);
	}
	
	// Save result VERTICES
	vector<int> vertex_result_id;
	vertex_result_id.resize(vertex_number,-1);
	int count = 0;
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			// Save coordinates and vertex qualities
			VectorType v = vtxQuad->target;
			float *vtx = new float[3];
			vtx[0] = v[0];
			vtx[1] = v[1];
			vtx[2] = v[2];
			vertices_result->push_back(vtx);
			if (mesh->extraTypeNumber != 0)
			{
				vtx_extra_result->at(0).push_back(vtxQuad->vertexQuality);
			}			
			vertex_result_id[i] = count;
			count++;
		}
	}

	if (flag == 1)
	{	// QEM connectivity
		for (uint i = 0; i != faceList->size(); ++i)
		{
			FaceInfo *info_f = faceList->at(i);
			if (info_f->is_exist())
			{
				int *t = new int[3];
				for (uint j = 0; j != 3; ++j)
				{
					t[j] = vertex_result_id[info_f->v[j]];	
				}
				face_result->push_back(t);
			}
		}
	}
	else
	{	// Cluster connectivity
		for (uint i = 0; i != faceList_new->size(); ++i)
		{
			FaceInfo *info_f = faceList_new->at(i);
			if (info_f->is_exist())
			{
				int *t = new int[3];
				for (uint j = 0; j != 3; ++j)
				{
					t[j] = vertex_result_id[info_f->v[j]];	
				}
				face_result->push_back(t);
			}
		}
	}
}

void ACVT::initialize_vertex_normals()
{
	////////////////////////////////////////////////////////
	// Compute Original Normals
	for (uint i = 0; i != vertex_number; ++i)
	{
		float *v = new float[3];
		for (uint j = 0; j != 3; ++j)
		{
			v[j] = 0;
		}
		normals_original->push_back(v);
	}

	for (uint i = 0; i != face_number; ++i)
	{
		int *v = mesh->getTriangle(i);
		float *v0 = mesh->getVertex(v[0]);
		float *v1 = mesh->getVertex(v[1]);
		float *v2 = mesh->getVertex(v[2]);
		float e0[3], e1[3], e2[3];
		for (uint j = 0; j != 3; ++j)
		{
			e0[j] = v1[j] - v0[j];
			e1[j] = v2[j] - v0[j];
		}
		crossProduct(e0, e1, e2);

		float *n0 = normals_original->at(v[0]);
		float *n1 = normals_original->at(v[1]);
		float *n2 = normals_original->at(v[2]);
		for (uint j = 0; j != 3; ++j)
		{
			n0[j] += e2[j];
			n1[j] += e2[j];
			n2[j] += e2[j];			
		}
	}
	for (uint i = 0; i != vertex_number; ++i)
	{
		float *v = normals_original->at(i);
		normalize(v);
	}


	////////////////////////////////////////////////////////
	// Initialize Current Normals
	for (uint i = 0; i != vertex_number; ++i)
	{
		float *v = new float[3];
		for (uint j = 0; j != 3; ++j)
		{
			v[j] = 1;
		}
		normals_current->push_back(v);
	}
}

void ACVT::compute_current_normals()
{
	for (uint i = 0; i != vertex_number; ++i)
	{
		float *v = normals_current->at(i);
		for (uint j = 0; j != 3; ++j)
		{
			v[j] = 0;
		}
	}
	if (!faceList_new->empty())
	{
		for (uint i = 0; i != faceList_new->size(); ++i)
		{
			FaceInfo* info = faceList_new->at(i);
			if (info->is_exist())
			{
				VectorType v0 = TargetVertex(info->v[0]);
				VectorType v1 = TargetVertex(info->v[1]);
				VectorType v2 = TargetVertex(info->v[2]);			
				float e0[3], e1[3], e2[3];
				for (uint j = 0; j != 3; ++j)
				{
					e0[j] = v1[j] - v0[j];
					e1[j] = v2[j] - v0[j];
				}
				crossProduct(e0, e1, e2);

				float *n0 = normals_current->at(info->v[0]);
				float *n1 = normals_current->at(info->v[1]);
				float *n2 = normals_current->at(info->v[2]);
				for (uint j = 0; j != 3; ++j)
				{
					n0[j] += e2[j];
					n1[j] += e2[j];
					n2[j] += e2[j];			
				}
			}
		}
		for (uint i = 0; i != vertex_number; ++i)
		{
			float *v = normals_current->at(i);
			normalize(v);
		}
	}
}

void ACVT::compute_QEM_normals()
{
	for (uint i = 0; i != vertex_number; ++i)
	{
		float *v = normals_current->at(i);
		for (uint j = 0; j != 3; ++j)
		{
			v[j] = 0;
		}
	}

	for(uint i = 0 ; i != faceList->size() ; i++)
	{
		FaceInfo *info_f = faceList->at(i);
		if (info_f->is_exist())
		{
			uint *t = info_f->v;
			VectorType v0 = TargetVertex(t[0]);
			VectorType v1 = TargetVertex(t[1]);
			VectorType v2 = TargetVertex(t[2]);
			float e0[3], e1[3], e2[3];
			for (uint j = 0; j != 3; ++j)
			{
				e0[j] = v1[j] - v0[j];
				e1[j] = v2[j] - v0[j];
			}
			crossProduct(e0, e1, e2);
			float *n0 = normals_current->at(t[0]);
			float *n1 = normals_current->at(t[1]);
			float *n2 = normals_current->at(t[2]);
			for (uint j = 0; j != 3; ++j)
			{
				n0[j] += e2[j];
				n1[j] += e2[j];
				n2[j] += e2[j];			
			}
		}
	}
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			float *v = normals_current->at(i);
			normalize(v);
		}
	}
}

void ACVT::compute_cluster_belong()
{
	cluster_belong->clear();
	ClearVector(*cluster_belong);
	cluster_belong = new vector<uint>;
	cluster_belong->resize(vertex_number, 0);
	//uint count = 0;
	//uint count_vtx = 0;
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			for (set<VertexID>::iterator iter = vtxQuad->cluster->begin();
				iter != vtxQuad->cluster->end(); ++iter)
			{
				//if (i == 0)
				//{
				//	count++;
				//	cout<<"Case 1: Vertex "<<i<<endl;
				//}
				cluster_belong->at(*iter) = i;
				//count_vtx++;
			}
		}
	}
	//cout<<"Vertex Number in all clusters:"<<count_vtx<<endl;
	//uint count1 = 0;
	//for (uint i = 0; i != vertex_number; ++i)
	//{
	//	if (cluster_belong->at(i) == 0)
	//	{
	//		cout<<"Case 2: Vertex "<<i<<endl;
	//		count1++;
	//	}
	//}
	//if (count != count1)
	//{
	//	cout<<"Cluster 0:"<<count<<endl;
	//	cout<<"Number of vertices belonging to cluster 0:"<<count1<<endl;
	//	cout<<"Warning: Something wrong in splitting the clusters."<<endl;
	//}
}

uint ACVT::get_vertex_number_in_clusters()
{
	uint count = 0;
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			for (set<VertexID>::iterator iter = vtxQuad->cluster->begin();
				iter != vtxQuad->cluster->end(); ++iter)
			{
				count++;
			}
		}
	}
	return count;
}

bool ACVT::is_cluster_border(VertexID vid)
{
	VertexQuadrics *vtxQuad = vtxQuadricsList->at(vid);
	VertexID clusterID = cluster_belong->at(vid);
	for (set<VertexID>::iterator iter = vtxQuad->neighbors->begin();
		iter != vtxQuad->neighbors->end(); ++iter)
	{
		if (cluster_belong->at(*iter) != clusterID)
		{
			// If there is one neighbor belonging to a different cluster,
			// then this vertex is a cluster border vertex.
			return true;
		}
	}
	return false;
}

double ACVT::compute_delta_energy(VertexID vid, VertexID i, VertexID j)
{
	VertexID clusterID_v0 = i, clusterID_v1 = j;
	MxQuadric Q0 = Quadrics(clusterID_v0);
	MxQuadric Q1 = Quadrics(clusterID_v1);
	MxQuadric Q_delta = InitialQuadrics(vid);
	double err_initial = Q0(TargetVertex(clusterID_v0));
	err_initial += Q1(TargetVertex(clusterID_v1));

	// Get result quadrics error after moving one vertex to another
	Q0 -= Q_delta;
	Q1 += Q_delta;
	VectorType newVtx0(dim());
	if (!Q0.optimize(newVtx0))
	{
		// if Q is invertible, use the mass center as target vertex
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(clusterID_v0);
		for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
			iter != vtxQuad->cluster->end(); ++iter)
		{
			if (*iter != vid)
			{
				VectorType v_temp(dim());
				getVertexWithExtraData(*iter, v_temp);
				newVtx0 += v_temp;
			}
		}
		newVtx0 /= (vtxQuad->cluster->size() - 1);
		cout<<"Warning: cluster "<<clusterID_v0<<"'s quadrics is invertible!"<<endl;
	}
	VectorType newVtx1(dim());
	if (!Q1.optimize(newVtx1))
	{
		// if Q is invertible, use the mass center as target vertex
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(clusterID_v1);
		for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
			iter != vtxQuad->cluster->end(); ++iter)
		{
			VectorType v_temp(dim());
			getVertexWithExtraData(*iter, v_temp);
			newVtx1 += v_temp;
		}
		VectorType vtx_vid(dim());
		getVertexWithExtraData(vid, vtx_vid);
		newVtx1 += vtx_vid;
		newVtx1 /= (vtxQuad->cluster->size() + 1);
		cout<<"Warning: cluster "<<clusterID_v1<<"'s quadrics is invertible!"<<endl;
	}
	// Get the result quadrics error
	double err_result = Q0(newVtx0);
	err_result += Q1(newVtx1);

	// Save the related information of the circumstance which results in the largest 
	// decrease of quadrics error after moving this vertex v0 to its neighbor cluster
	return (err_initial - err_result);
}

bool ACVT::swap_once()
{
	/***********************************************************************************/
	/* 
	* Here we will find all the cluster border vertices which could result in an energy 
	* decrease if we swap them into their neighbor clusters.
	* Note: here we firstly find all the qualified vertices and then SWAP THEM ALL AT ONCE.
	*/
	
	/********************************************/
	/*
	* Firstly, get the vertices to be swapped
	*/
	/* Get all border points to be swapped */
	vector<VertexID> *cluster_belong_new = new vector<VertexID>;
	for (uint i = 0; i != vertex_number; ++i)
	{
		cluster_belong_new->push_back(cluster_belong->at(i));
	}
	uint count = 0;
	uint count_border = 0;
	vector<VertexID> vertex_to_be_swapped;
	set<VertexID> clusters_to_be_updated;
	//double t = get_cpu_time();
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			if (vtxQuad->cluster->size() > 1)
			{
				// For each vertex in cluster i, determine if it is a cluster border vertex
				uint clusterID_v0 = i;
				for (set<VertexID>::iterator iter = vtxQuad->cluster->begin();
					iter != vtxQuad->cluster->end(); ++iter)
				{
					// For a vertex v0
					VertexID v0_id = *iter;
					VertexQuadrics *v0Quad = vtxQuadricsList->at(v0_id);
					VertexID cluster_opt = clusterID_v0;
					double energy_saved_most = 0;
					bool border_flag = false;
					for (set<VertexID>::iterator iter1 = v0Quad->neighbors->begin();
						iter1 != v0Quad->neighbors->end(); ++iter1)
					{
						// If one of v0's neighbor belongs to a different cluster,
						// then it is a border vertex
						VertexID v1_id = *iter1;
						uint clusterID_v1 = cluster_belong->at(v1_id);
						if (clusterID_v1 != clusterID_v0)
						{
							border_flag = true;
							break;
						}
					}
					if (border_flag)
					{
						count_border++;
						/*
						* Swap this cluster border vertex v0. We have two methods:
						* 1) Swap v0 to a neighbor cluster of v0's cluster
						* 2) Swap v0 to a neighbor cluster of v0
						* Here I use the first one. Please distinguish the difference between the two.
						* Both methods are OK, but the result will be a little different.
						*/ 
						for (set<VertexID>::iterator iter1 = neighbor_clusters->at(clusterID_v0).begin();
							iter1 != neighbor_clusters->at(clusterID_v0).end(); ++iter1)
						{
							uint clusterID_v1 = *iter1;
							double energy_delta = compute_delta_energy(v0_id, clusterID_v0, clusterID_v1);
							if (energy_delta - energy_saved_most > 1e-8)
							{
								energy_saved_most = energy_delta;
								cluster_opt = clusterID_v1;
							}
						}
					}
					// After checking all the neighbors of v0, find the best neighbor cluster
					// and assign it as the new cluster of v0
					if (cluster_opt != clusterID_v0)
					{
						count++;
						vertex_to_be_swapped.push_back(v0_id);
						clusters_to_be_updated.insert(clusterID_v0);
						clusters_to_be_updated.insert(cluster_opt);
					}
					cluster_belong_new->at(v0_id) = cluster_opt;
				}
			}
		}
	}
	cout<<"\t"<<count<<" vertices will be swapped"<<endl;
	
	/********************************************/
	/*
	* Secondly, swap these vertices
	*/
	for (uint i = 0; i != vertex_to_be_swapped.size(); ++i)
	{
		// If a vertex's new cluster id is different from old, then we should move it to the new cluster
		VertexID v_id = vertex_to_be_swapped[i];
		VertexID clusterID_v0 = cluster_belong->at(v_id);
		VertexID clusterID_v1 = cluster_belong_new->at(v_id);

		// Remove the vertex from its old cluster
		VertexQuadrics *vtxQuad_old = vtxQuadricsList->at(clusterID_v0);
		vtxQuad_old->delete_vertex_from_cluster(v_id);
		Quadrics(clusterID_v0) -= InitialQuadrics(v_id);

		// Add this vertex into its new cluster
		VertexQuadrics *vtxQuad_new = vtxQuadricsList->at(clusterID_v1);
		vtxQuad_new->add_vertex_to_cluster(v_id);		
		Quadrics(clusterID_v1) += InitialQuadrics(v_id);
	}
	// Update the target vertex coordinates for each cluster (i.e. vertex in result mesh)
	for (set<VertexID>::iterator iter = clusters_to_be_updated.begin(); 
		iter != clusters_to_be_updated.end(); ++iter)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(*iter);
		if (!vtxQuad->cluster->empty())
		{
			compute_target(*iter);
		}
	}

	/***********************************************************************************/
	/*
	* Next, there may exist some cluster border points that if we swap them back to their
	* original clusters, the total energy will decrease. Then, we have to check and find these 
	* vertices and swap them back ONE BY ONE. Here we swap these points back in natural number order
	* from smallest index to largest index.
	*/
	uint count1 = 0;
	vector<VertexID> vtx_not_swapped;
	for (uint i = 0; i != vertex_to_be_swapped.size(); ++i)
	{
		// If a vertex's new cluster id is different from old, then we should move it to the new cluster
		VertexID v_id = vertex_to_be_swapped[i];
		VertexID clusterID_v0 = cluster_belong_new->at(v_id);
		VertexID clusterID_v1 = cluster_belong->at(v_id);

		assert(clusterID_v0 != clusterID_v1);

		// If the total energy will decrease after swapping this vertex back,
		// then swap this vertex back by updating the connectivity and cluster-belong information
		// Remove the vertex from its old cluster
		double energy_delta = compute_delta_energy(v_id, clusterID_v0, clusterID_v1);
		if (energy_delta > 0)
		{
			// Firstly, update the two related clusters
			VertexQuadrics *vtxQuad_old = vtxQuadricsList->at(clusterID_v0);
			vtxQuad_old->delete_vertex_from_cluster(v_id);
			Quadrics(clusterID_v0) -= InitialQuadrics(v_id);
			compute_target(clusterID_v0);

			VertexQuadrics *vtxQuad_new = vtxQuadricsList->at(clusterID_v1);
			vtxQuad_new->add_vertex_to_cluster(v_id);
			Quadrics(clusterID_v1) += InitialQuadrics(v_id);
			compute_target(clusterID_v1);

			count1++;
			cluster_belong_new->at(v_id) = cluster_belong->at(v_id);
			vtx_not_swapped.push_back(i);
		}
	}
	//
	// Delete those vertices swapped back and only keep those vertices to be swapped 
	//
	for (vector<VertexID>::reverse_iterator iter = vtx_not_swapped.rbegin(); 
		iter != vtx_not_swapped.rend(); ++iter)
	{
		vertex_to_be_swapped.erase(vertex_to_be_swapped.begin() + (*iter));
	}
	/*
	* Check if there are such clusters whose vertices are all swapped outside.
	* If so, reinstate these clusters to their original states, since we
	* want to preserve the total vertex number during this swapping step.
	*/
	for (set<VertexID>::iterator iter = clusters_to_be_updated.begin();
		iter != clusters_to_be_updated.end(); ++iter)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(*iter);
		if (vtxQuad->cluster->empty())
		{
			for (uint i = 0; i != vertex_to_be_swapped.size(); ++i)
			{
				VertexID v_id = vertex_to_be_swapped[i];
				VertexID clusterID_v0 = cluster_belong->at(v_id);
				VertexID clusterID_v1 = cluster_belong_new->at(v_id);
				// Get all the swapping-outside vertices and reinstate this cluster by
				// swapping these vertices back to their original clusters
				if (clusterID_v0 == *iter && clusterID_v1 != *iter)
				{ 
					VertexQuadrics *vtxQuad_old = vtxQuadricsList->at(clusterID_v1);
					vtxQuad_old->delete_vertex_from_cluster(v_id);
					Quadrics(clusterID_v1) -= InitialQuadrics(v_id);
					compute_target(clusterID_v1);

					VertexQuadrics *vtxQuad_new = vtxQuadricsList->at(clusterID_v0);
					vtxQuad_new->add_vertex_to_cluster(v_id);
					Quadrics(clusterID_v0) += InitialQuadrics(v_id);
					
					cluster_belong_new->at(v_id) = cluster_belong->at(v_id);
				}
			}
			compute_target(*iter);
		}
	}
	// 
	// Update cluster belong
	//
	for (uint i = 0; i != vertex_to_be_swapped.size(); ++i)
	{
		VertexID v_id = vertex_to_be_swapped[i];
		cluster_belong->at(v_id) = cluster_belong_new->at(v_id);
	}
	//
	// Free memories
	//
	cluster_belong_new->clear();
	ClearVector(*cluster_belong_new);
	vertex_to_be_swapped.clear();
	ClearVector(vertex_to_be_swapped);
	vtx_not_swapped.clear();
	ClearVector(vtx_not_swapped);
	cout<<"\t"<<count1<<" vertices will be swapped back"<<endl;
	//
	// If the swapping back number is smaller than the swapping front number, then 
	// this swap step is good.
	//
	if (count1 < count)
	{		
		// Update the neighbor clusters of related clusters
		set<VertexID> setTemp;
		for (set<VertexID>::iterator iter = clusters_to_be_updated.begin();
			iter != clusters_to_be_updated.end(); ++iter)
		{
			for (set<VertexID>::iterator iter0 = neighbor_clusters->at(*iter).begin(); 
				iter0 != neighbor_clusters->at(*iter).end(); ++iter0)
			{
				setTemp.insert(*iter0);
			}
		}
		for (set<VertexID>::iterator iter = setTemp.begin();
			iter != setTemp.end(); ++iter)
		{
			clusters_to_be_updated.insert(*iter);
		}
		for (set<VertexID>::iterator iter = clusters_to_be_updated.begin();
			iter != clusters_to_be_updated.end(); ++iter)
		{
			compute_neighbor_cluster(*iter, neighbor_clusters->at(*iter));
		}
		setTemp.clear();
		ClearSet(setTemp);

		clusters_to_be_updated.clear();
		ClearSet(clusters_to_be_updated);
		return true;
	}
	else
	{
		// No need to swap again.
		clusters_to_be_updated.clear();
		ClearSet(clusters_to_be_updated);
		return false;
	}
}


void ACVT::compute_all_neighbor_clusters()
{
	// Firstly, compute the cluster ids each vertex belongs to
	compute_cluster_belong();

	// Next, compute the neighbor clusters for each cluster
	// Clear relative data structure
	if (!neighbor_clusters->empty())
	{
		for (uint i = 0; i != vertex_number; ++i)
		{
			neighbor_clusters->at(i).clear();
		}
		neighbor_clusters->clear();
		ClearVector(*neighbor_clusters);
		neighbor_clusters = new vector<set<VertexID>>;
	}
	neighbor_clusters->resize(vertex_number);

	for (uint i = 0; i != vertex_number; ++i)
	{
		compute_neighbor_cluster(i, neighbor_clusters->at(i));
	}
}

void ACVT::compute_neighbor_cluster(VertexID i, set<VertexID>& neighbors)
{
	VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
	neighbors.clear();
	if (!vtxQuad->cluster->empty())
	{
		VertexID v0_clusterID = i;
		for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
			iter != vtxQuad->cluster->end(); ++iter)
		{
			VertexID v0_id = *iter;
			VertexQuadrics *v0_vtxQuad = vtxQuadricsList->at(v0_id);

			// For vertex v0, detect all its neighbors to check if there exists
			// such neighbors belonging to different clusters
			for (set<VertexID>::iterator iter0 = v0_vtxQuad->neighbors->begin(); 
				iter0 != v0_vtxQuad->neighbors->end(); ++iter0)
			{
				VertexID v1_id = *iter0;
				VertexID v1_clusterID = cluster_belong->at(v1_id);
				if (v1_clusterID != v0_clusterID)
				{
					neighbors.insert(v1_clusterID);
				}
			}
		}
	}
}

void ACVT::compute_new_edge_pairs()
{
	edgePairs_new.clear();
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			for (set<VertexID>::iterator iter = neighbor_clusters->at(i).begin();
				iter != neighbor_clusters->at(i).end(); ++iter)
			{
				VertexID a = i;
				VertexID b = *iter;
				sort2(a,b);
				edgePairs_new[pair<uint,uint>(a,b)]++;
			}
		}
	}
}

void ACVT::compute_new_faces()
{
	/*
	* Clear the old objects
	*/
	if (!faceList_new->empty())
	{
		for (vector<FaceInfo*>::iterator iter = faceList_new->begin();
			iter != faceList_new->end(); ++iter)
		{
			delete *iter;
		}
		faceList_new->clear();
		ClearVector(*faceList_new);
	}
	faceList_new = new vector<FaceInfo*>;
	if (!halfEdges->empty())
	{
		for (uint i = 0; i != halfEdges->size(); ++i)
		{
			for (uint j = 0; j != halfEdges->at(i)->size(); ++j)
			{
				delete halfEdges->at(i)->at(j);
			}
			halfEdges->at(i)->clear();
			ClearVector(*(halfEdges->at(i)));
		}
		halfEdges->clear();
		ClearVector(*halfEdges);
	}
	halfEdges = new vector<vector<HalfEdge*>*>;

	/* Firstly compute all cluster neighbors and new edges */
	compute_all_neighbor_clusters();
	compute_new_edge_pairs();	// save all new edges in "edgePairs_new" object

	/* 
	* Compute a "cluster" normal for each cluster vertex. Here the normal for each cluster vertex 
	* is the average normal of all cluster vertices from original mesh. This normal will be used
	* to check if a created face is reversed from original direction.
	*/
	for (uint i = 0; i != vertex_number; ++i)
	{
		float *v = normals_current->at(i);
		for (uint j = 0; j != 3; ++j)
		{
			v[j] = 0;
		}
	}
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			// Save the new "cluster normal" in "normal_current" Object
			float *nn = getCurrentNormal(i);
			// A cluster normal is the average normal of all cluster vertices
			for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
				iter != vtxQuad->cluster->end(); ++iter)
			{
				float *n = getOriginalNormal(*iter);
				for (uint j = 0; j != 3; ++j)
				{
					nn[j] += n[j];
				}
			}
			normalize(nn);
		}
	}

	/* 
	* Secondly, get all half edges
	*/
	for (uint i = 0; i != vertex_number; ++i)
	{
		vector<HalfEdge*>* vecTemp = new vector<HalfEdge*>;
		halfEdges->push_back(vecTemp);
	}
	HalfEdge *halfedge_now;
	for (EdgePairType::iterator iter = edgePairs_new.begin(); 
		iter != edgePairs_new.end(); ++iter)
	{
		// Each edge denotes two half edges with inverse directions
		/* forward half edge */
		halfedge_now = new HalfEdge;
		halfedge_now->v1 = (*iter).first.first;
		halfedge_now->v2 = (*iter).first.second;
		halfedge_now->inside_face_flag = false;
		halfedge_now->visited_flag = false;
		halfEdges->at(halfedge_now->v1)->push_back(halfedge_now);

		/* backward half edge */
		halfedge_now = new HalfEdge;
		halfedge_now->v1 = (*iter).first.second;
		halfedge_now->v2 = (*iter).first.first;
		halfedge_now->inside_face_flag = false;
		halfedge_now->visited_flag = false;
		halfEdges->at(halfedge_now->v1)->push_back(halfedge_now);
	}

	/* 
	* Next, Run half-edge algorithm on all half edges to determine trianglar faces.
	*/
	for (uint i = 0; i != halfEdges->size(); ++i)
	{
		if (halfEdges->at(i) != NULL)
		{
			for (uint j = 0; j != halfEdges->at(i)->size(); ++j)
			{
				HalfEdge *halfedge_start = halfEdges->at(i)->at(j);
				if (!halfedge_start->inside_face_flag)
				{
					// Firstly, run Recursion step on each unused half edge
					run_half_edge_method(halfedge_start, halfEdges);

					// Next, check if there exist unqualified faces to be split.
					bool flag_runrecursion = true;
					while (flag_runrecursion)
					{
						flag_runrecursion = false;
						
						// Get unused half edges that are not inside triangles
						vector<HalfEdge*> *halfedges_unused = new vector<HalfEdge*>;
						for (uint i = 0; i != halfEdges->size(); ++i)
						{
							if (halfEdges->at(i) != NULL)
							{
								for (uint j = 0; j != halfEdges->at(i)->size(); ++j)
								{
									HalfEdge *halfedge_now = halfEdges->at(i)->at(j);
									if (!halfedge_now->inside_face_flag)
									{
										halfedges_unused->push_back(halfedge_now);
										//HalfEdge* halfedge_opposite = find_opposite_halfedge(halfedge_now, halfEdges);
										//if (!halfedge_opposite->inside_face_flag)
										//{
										//	halfedges_unused->push_back(halfedge_now);
										//}
									}
								}
							}
						}
						//cout<<"Half edges that are not inside faces: "<<halfedges_unused->size()<<endl;
						
						// Get all triangles containing endpoints of unused half edges
						set<VertexID> endpoints, faces_adjacent;
						for (uint i = 0; i != halfedges_unused->size(); ++i)
						{
							HalfEdge *halfedge_now = halfedges_unused->at(i);
							endpoints.insert(halfedge_now->v1);
							endpoints.insert(halfedge_now->v2);
						}
						for (uint i = 0; i != faceList_new->size(); ++i)
						{
							FaceInfo* face_now = faceList_new->at(i);
							if (face_now->is_exist())
							{
								for (uint j = 0; j != 3; ++j)
								{
									set<VertexID>::iterator iter = endpoints.find(face_now->v[j]);
									if (iter != endpoints.end())
									{
										faces_adjacent.insert(i);
										break;
									}
								}			
							}
						}
						
						// For each adjacent face, test if this face can be split to three half edges, meanwhile 
						// each half edge can compose a new face if combined with another two
						// unused half edges in the unused half edge list obtained above.
						vector<VertexID> faces_to_be_deleted;
						for (set<VertexID>::iterator iter = faces_adjacent.begin(); iter != faces_adjacent.end(); ++iter)
						{
							FaceInfo* face_now = faceList_new->at(*iter);
							HalfEdge* halfedge_1 = get_half_edge(face_now->v[0],face_now->v[1], halfEdges);
							HalfEdge* halfedge_2 = get_half_edge(face_now->v[1],face_now->v[2], halfEdges);
							HalfEdge* halfedge_3 = get_half_edge(face_now->v[2],face_now->v[0], halfEdges);
							if (halfedge_1 != NULL && halfedge_2 != NULL && halfedge_3 != NULL)
							{
								VertexID halfedge_1_2_index = 0, halfedge_1_3_index = 0;
								bool flag1 = check_halfedge_form_new_face(halfedge_1, halfedges_unused, halfedge_1_2_index, halfedge_1_3_index);
								VertexID halfedge_2_2_index = 0, halfedge_2_3_index = 0;
								bool flag2 = check_halfedge_form_new_face(halfedge_2, halfedges_unused, halfedge_2_2_index, halfedge_2_3_index);
								VertexID halfedge_3_2_index = 0, halfedge_3_3_index = 0;
								bool flag3 = check_halfedge_form_new_face(halfedge_3, halfedges_unused, halfedge_3_2_index, halfedge_3_3_index);
								if (flag1 == true && flag2 == true && flag3 == true)
								{
									cout<<"Face to be deleted: ("<<face_now->v[0]<<","<<face_now->v[1]<<","<<face_now->v[2]<<")"<<endl;
									faces_to_be_deleted.push_back(*iter);

									HalfEdge *halfedge_now = halfedges_unused->at(halfedge_1_2_index);
									FaceInfo* newFace = new FaceInfo(halfedge_1->v1, halfedge_1->v2, halfedge_now->v2);
									faceList_new->push_back(newFace);
									halfedge_now->inside_face_flag = true;
									halfedge_now = halfedges_unused->at(halfedge_1_3_index);
									halfedge_now->inside_face_flag = true;
									cout<<"Face to be created: ("<<newFace->v[0]<<","<<newFace->v[1]<<","<<newFace->v[2]<<")"<<endl;

									halfedge_now = halfedges_unused->at(halfedge_2_2_index);
									newFace = new FaceInfo(halfedge_2->v1, halfedge_2->v2, halfedge_now->v2);
									faceList_new->push_back(newFace);
									halfedge_now->inside_face_flag = true;
									halfedge_now = halfedges_unused->at(halfedge_2_3_index);
									halfedge_now->inside_face_flag = true;
									cout<<"Face to be created: ("<<newFace->v[0]<<","<<newFace->v[1]<<","<<newFace->v[2]<<")"<<endl;

									halfedge_now = halfedges_unused->at(halfedge_3_2_index);
									newFace = new FaceInfo(halfedge_3->v1, halfedge_3->v2, halfedge_now->v2);
									faceList_new->push_back(newFace);
									halfedge_now->inside_face_flag = true;
									halfedge_now = halfedges_unused->at(halfedge_3_3_index);
									halfedge_now->inside_face_flag = true;
									cout<<"Face to be created: ("<<newFace->v[0]<<","<<newFace->v[1]<<","<<newFace->v[2]<<")"<<endl<<endl;

									flag_runrecursion = true;
								}
							}
						}
						// Delete triangular faces which need to be split into three half edges
						if (flag_runrecursion)
						{
							for (vector<VertexID>::reverse_iterator iter = faces_to_be_deleted.rbegin();
								iter != faces_to_be_deleted.rend(); ++iter)
							{
								FaceInfo* facenow = faceList_new->at(*iter);
								faceList_new->erase(faceList_new->begin()+(*iter));
								delete facenow;
							}

							// Run this recursion again on the same half edge to update this connected component 
							// by this half edge
							run_half_edge_method(halfedge_start, halfEdges);
						}
						halfedges_unused->clear();	
						ClearVector(*halfedges_unused);
						endpoints.clear();
						ClearSet(endpoints);
						faces_adjacent.clear();
						ClearSet(faces_adjacent);
						faces_to_be_deleted.clear();
						ClearVector(faces_to_be_deleted);
					}
				}
			}
		}
	}

	cout<<"Vertex Number: "<<get_current_vertex_number()<<endl;
	cout<<"Face Number: "<<faceList_new->size()<<endl;

	/*
	* Compute current normals and free memory
	*/
	compute_current_normals();
}

void ACVT::compute_new_faces_nonrecursive()
{
	/*
	* Clear the old objects
	*/
	if (!faceList_new->empty())
	{
		for (vector<FaceInfo*>::iterator iter = faceList_new->begin();
			iter != faceList_new->end(); ++iter)
		{
			delete *iter;
		}
		faceList_new->clear();
		ClearVector(*faceList_new);
	}
	faceList_new = new vector<FaceInfo*>;
	if (!halfEdges->empty())
	{
		for (uint i = 0; i != halfEdges->size(); ++i)
		{
			for (uint j = 0; j != halfEdges->at(i)->size(); ++j)
			{
				delete halfEdges->at(i)->at(j);
			}
			halfEdges->at(i)->clear();
			ClearVector(*(halfEdges->at(i)));
		}
		halfEdges->clear();
		ClearVector(*halfEdges);
	}
	halfEdges = new vector<vector<HalfEdge*>*>;

	/* 
	* Firstly compute all cluster neighbors and new edges 
	*/
	compute_all_neighbor_clusters();
	compute_new_edge_pairs();	// save all new edges in "edgePairs_new" object
	/* 
	* Compute a "cluster" normal for each cluster vertex. Here the normal for each cluster vertex 
	* is the average normal of all cluster vertices from original mesh. This normal will be used
	* to check if a created face is reversed from original direction.
	*/
	for (uint i = 0; i != vertex_number; ++i)
	{
		float *v = normals_current->at(i);
		for (uint j = 0; j != 3; ++j)
		{
			v[j] = 0;
		}
	}
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			// Save the new "cluster normal" in "normal_current" Object
			float *nn = getCurrentNormal(i);
			// A cluster normal is the average normal of all cluster vertices
			for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
				iter != vtxQuad->cluster->end(); ++iter)
			{
				float *n = getOriginalNormal(*iter);
				for (uint j = 0; j != 3; ++j)
				{
					nn[j] += n[j];
				}
			}
			normalize(nn);
		}
	}	
	/* 
	* Secondly, get all half edges
	*/
	for (uint i = 0; i != vertex_number; ++i)
	{
		vector<HalfEdge*>* vecTemp = new vector<HalfEdge*>;
		halfEdges->push_back(vecTemp);
	}
	HalfEdge *halfedge_now;
	for (EdgePairType::iterator iter = edgePairs_new.begin(); 
		iter != edgePairs_new.end(); ++iter)
	{
		/* Each edge denotes two half edges with inverse directions */
		/* forward half edge */
		halfedge_now = new HalfEdge;
		halfedge_now->v1 = (*iter).first.first;
		halfedge_now->v2 = (*iter).first.second;
		halfedge_now->inside_face_flag = false;
		halfedge_now->visited_flag = false;
		halfEdges->at(halfedge_now->v1)->push_back(halfedge_now);

		/* backward half edge */
		halfedge_now = new HalfEdge;
		halfedge_now->v1 = (*iter).first.second;
		halfedge_now->v2 = (*iter).first.first;
		halfedge_now->inside_face_flag = false;
		halfedge_now->visited_flag = false;
		halfEdges->at(halfedge_now->v1)->push_back(halfedge_now);
	}

	/*
	* Next, run half-edge algorithm on each half edge
	*/
	for (uint i = 0; i != halfEdges->size(); ++i)
	{
		if (halfEdges->at(i) != NULL)
		{
			for (uint j = 0; j != halfEdges->at(i)->size(); ++j)
			{
				HalfEdge *halfedge_start = halfEdges->at(i)->at(j);
				run_half_edge_method_nonrecursive(halfedge_start);
			}
		}
	}
	cout<<"Vertex Number: "<<get_current_vertex_number()<<endl;
	cout<<"Face Number: "<<faceList_new->size()<<endl;
	/*
	* Finally, compute current normals
	*/
	compute_current_normals();
}

void ACVT::run_half_edge_method(HalfEdge *halfedge_start, vector<vector<HalfEdge*>*>* halfEdges_all)
{
	// If all half-edges are used, return true
	if (halfedge_start->inside_face_flag == true || halfedge_start->visited_flag == true)
		return;

	// All 6 half-edges
	HalfEdge* halfedge_1 = halfedge_start;	// (v1->v2)
	HalfEdge* halfedge_2 = NULL;			// (v2->v3)
	HalfEdge* halfedge_3 = NULL;			// (v3->v1)
	HalfEdge* halfedge_1_backward = NULL;	// (v2->v1)
	HalfEdge* halfedge_2_backward = NULL;	// (v3->v2)
	HalfEdge* halfedge_3_backward = NULL;	// (v1->v3)

	VertexID v1 = halfedge_1->v1;	// We have v1->v2 initially
	VertexID v2 = halfedge_1->v2;
	VertexID v3 = 0;
	bool flag = false;
	for (uint i = 0; i != halfEdges_all->at(v2)->size(); ++i)
	{
		halfedge_2 = halfEdges_all->at(v2)->at(i);
		v3 = halfedge_2->v2;
		if (v3 == v1 || halfedge_2->inside_face_flag == true)
			continue;
		// Now we have v2->v3, and needs to find v3->v1
		for (uint j = 0; j != halfEdges_all->at(v3)->size(); ++j)
		{
			halfedge_3 = halfEdges_all->at(v3)->at(j);
			if (halfedge_3->v2 == v2 || halfedge_3->inside_face_flag == true)
				continue;

			// If these three edges compose a loop, i.e. a new triangle face
			if (halfedge_3->v2 == v1)
			{
				// Here we must firstly check if this face has already been created.
				// Note: there is only one chance that this new face A is already created:
				//	there is a same face just created in the last iteration.
				if (faceList_new->size() != 0)
				{
					// Check the face just created in the last iteration
					FaceInfo* info = faceList_new->at(faceList_new->size()-1);
					VertexID a = info->v[0], b = info->v[1], c = info->v[2];
					sort3(a,b,c);
					VertexID aa = v1, bb = v2, cc = v3;
					sort3(aa,bb,cc);
					if (aa == a && bb == b && cc == c)
						continue;
				}
				// If the face is different from the last face, then it is a real face. Quit
				flag = true;
				break;
			}
		}
		if (flag)
		{// Have found a new face, quit.
			break;
		}
	}

	// Find the 1st backward half edge 
	halfedge_1_backward = find_opposite_halfedge(halfedge_1, halfEdges_all);

	if (flag)
	{// If a new triangular face is found with the input half edge, then test if its
		// availability and save it.

		// Test if this face's normal is reversed 
		VectorType v1_target = TargetVertex(v1), v2_target = TargetVertex(v2), v3_target = TargetVertex(v3);
		float n[3];
		float e1[3], e2[3];
		for (uint i = 0; i != 3; ++i)
		{
			e1[i] = v2_target[i] - v1_target[i];
			e2[i] = v3_target[i] - v1_target[i];
		}
		crossProduct(e1,e2,n);
		float *n1 = getCurrentNormal(v1), *n2 = getCurrentNormal(v2), *n3 = getCurrentNormal(v3);
		if (dot3(n,n1) < thresVal_normal_direction 
			&& dot3(n,n2) < thresVal_normal_direction 
			&& dot3(n,n3) < thresVal_normal_direction)
		{
			return;
		}

		// If a new face is found, then create this face and run recursion
		// on the three backward half-edges
		FaceInfo* info = new FaceInfo(v1, v2, v3);
		faceList_new->push_back(info);
		halfedge_1->inside_face_flag = true;
		halfedge_2->inside_face_flag = true;
		halfedge_3->inside_face_flag = true;
		//cout<<v2<<" "<<v3<<endl;
		//cout<<v3<<" "<<v1<<endl;

		// run recursion process on the first reverse half edge
		run_half_edge_method(halfedge_1_backward, halfEdges_all);

		// Find the 2nd and 3rd backward half edges
		halfedge_2_backward = find_opposite_halfedge(halfedge_2, halfEdges_all);
		halfedge_3_backward = find_opposite_halfedge(halfedge_3, halfEdges_all);
		run_half_edge_method(halfedge_2_backward, halfEdges_all);
		run_half_edge_method(halfedge_3_backward, halfEdges_all);
	}
	else
	{
		// If a new triangular face is not found, then this input half-edge is very likely to be a 
		// border half edge outside the boundary. Then, only run recursion on the 1st backward half-edge.
		//halfedge_1->inside_face_flag = true;
		halfedge_1->visited_flag = true;
		run_half_edge_method(halfedge_1_backward, halfEdges_all);
	}
}

void ACVT::run_half_edge_method_nonrecursive(HalfEdge *halfedge_start)
{
	if (halfedge_start == NULL)
	{
		return;
	}
	if (halfedge_start->inside_face_flag == true || halfedge_start->visited_flag == true)
	{
		return;
	}

	// Use a queue to save half edges to be visited
	stack<pair<VertexID,VertexID>> stack_tobevisited_halfedges;
	stack_tobevisited_halfedges.push(make_pair(halfedge_start->v1, halfedge_start->v2));
	while (!stack_tobevisited_halfedges.empty())
	{
		pair<VertexID,VertexID> halfedge_temp = stack_tobevisited_halfedges.top();
		HalfEdge* halfedge_1 = get_half_edge(halfedge_temp.first, halfedge_temp.second, halfEdges);// half edge (v1->v2)
		/* 
		* For halfedge_1, check if there exist another two half edges which can compose a triangular face 
		*/
		HalfEdge* halfedge_2 = NULL;	// (v2->v3)
		HalfEdge* halfedge_3 = NULL;	// (v3->v1)
		VertexID v1 = halfedge_1->v1;	// We have half edge v1->v2 initially
		VertexID v2 = halfedge_1->v2;
		VertexID v3 = 0;
		bool flag = false;
		for (uint i = 0; i != halfEdges->at(v2)->size(); ++i)
		{
			halfedge_2 = halfEdges->at(v2)->at(i);
			v3 = halfedge_2->v2;
			if (v3 == v1 || halfedge_2->inside_face_flag == true)
				continue;
			// Now we have half edge v2->v3, and needs to find v3->v1
			for (uint j = 0; j != halfEdges->at(v3)->size(); ++j)
			{
				halfedge_3 = halfEdges->at(v3)->at(j);
				if (halfedge_3->v2 == v2 || halfedge_3->inside_face_flag == true)
					continue;
				// If these three edges compose a loop, i.e. a new triangle face
				if (halfedge_3->v2 == v1)
				{
					/*
					* Now we have half edge v3->v1 and a triangular face, but we must 
					* firstly check if this face has already been created and saved.
					* Note that if this face has already been created, it must be created
					* in the last iteration and saved in the bottom of "faceList_new".
					*/
					if (faceList_new->size() != 0)
					{
						// Check the face just created in the last iteration
						FaceInfo* info = faceList_new->at(faceList_new->size()-1);
						VertexID a = info->v[0], b = info->v[1], c = info->v[2];
						sort3(a,b,c);
						VertexID aa = v1, bb = v2, cc = v3;
						sort3(aa,bb,cc);
						if (aa == a && bb == b && cc == c)
							continue;
					}
					// If the face is new, set flag and quit.
					flag = true;
					break;
				}
			}
			if (flag)
			{// Have found a new face, quit.
				break;
			}
		}

		// Find the 1st backward half edge and push it into the queue
		halfedge_1->visited_flag = true;
		stack_tobevisited_halfedges.pop();
		HalfEdge* halfedge_1_backward = find_opposite_halfedge(halfedge_1, halfEdges);
		if (flag) 
		{
			/* 
			* If found a new face, firstly test if this face's normal is reversed. If not, create this face.
			*/
			//VectorType v1_target = TargetVertex(v1), v2_target = TargetVertex(v2), v3_target = TargetVertex(v3);
			//float n[3];
			//float e1[3], e2[3];
			//for (uint i = 0; i != 3; ++i)
			//{
			//	e1[i] = v2_target[i] - v1_target[i];
			//	e2[i] = v3_target[i] - v1_target[i];
			//}
			//crossProduct(e1,e2,n);
			//float *n1 = getCurrentNormal(v1), *n2 = getCurrentNormal(v2), *n3 = getCurrentNormal(v3);
			//if (dot3(n,n1) > thresVal_normal_direction || dot3(n,n2) > thresVal_normal_direction 
			//	|| dot3(n,n3) > thresVal_normal_direction)
			//{
			//	/* 
			//	* If one of the face's three vertex normals is not reversed, then the face is not reversed.
			//	* Then create this face.
			//	*/
			//	FaceInfo* info = new FaceInfo(v1, v2, v3);
			//	faceList_new->push_back(info);
			//	halfedge_1->inside_face_flag = true;
			//	halfedge_2->inside_face_flag = true;
			//	halfedge_2->visited_flag = true;
			//	halfedge_3->inside_face_flag = true;
			//	halfedge_3->visited_flag = true;
			//	// Find 3rd backward half edge into the stack first, the 2nd next, the 1st at last
			//	HalfEdge* halfedge_3_backward = find_opposite_halfedge(halfedge_3, halfEdges);
			//	if (halfedge_3_backward != NULL)
			//	{
			//		if(!halfedge_3_backward->visited_flag)
			//		{
			//			stack_tobevisited_halfedges.push(make_pair(halfedge_3_backward->v1, halfedge_3_backward->v2));
			//		}
			//	}				
			//	HalfEdge* halfedge_2_backward = find_opposite_halfedge(halfedge_2, halfEdges);
			//	if (halfedge_2_backward != NULL)
			//	{
			//		if(!halfedge_2_backward->visited_flag)
			//		{
			//			stack_tobevisited_halfedges.push(make_pair(halfedge_2_backward->v1, halfedge_2_backward->v2));
			//		}
			//	}
			//}

			FaceInfo* info = new FaceInfo(v1, v2, v3);
			faceList_new->push_back(info);
			halfedge_1->inside_face_flag = true;
			halfedge_2->inside_face_flag = true;
			halfedge_2->visited_flag = true;
			halfedge_3->inside_face_flag = true;
			halfedge_3->visited_flag = true;
			// Find 3rd backward half edge into the stack first, the 2nd next, the 1st at last
			HalfEdge* halfedge_3_backward = find_opposite_halfedge(halfedge_3, halfEdges);
			if (halfedge_3_backward != NULL)
			{
				if(!halfedge_3_backward->visited_flag)
				{
					stack_tobevisited_halfedges.push(make_pair(halfedge_3_backward->v1, halfedge_3_backward->v2));
				}
			}				
			HalfEdge* halfedge_2_backward = find_opposite_halfedge(halfedge_2, halfEdges);
			if (halfedge_2_backward != NULL)
			{
				if(!halfedge_2_backward->visited_flag)
				{
					stack_tobevisited_halfedges.push(make_pair(halfedge_2_backward->v1, halfedge_2_backward->v2));
				}
			}
		}
		if (halfedge_1_backward != NULL)
		{
			if (!halfedge_1_backward->visited_flag)
			{
				stack_tobevisited_halfedges.push(make_pair(halfedge_1_backward->v1, halfedge_1_backward->v2));
			}
		}
	}
}

bool ACVT::check_halfedge_form_new_face(HalfEdge* halfedge_1, vector<HalfEdge*> *halfedge_list,	VertexID& halfedge2_index, VertexID& halfedge3_index)
{
	VertexID v1 = halfedge_1->v1, v2 = halfedge_1->v2;
	VertexID v3;
	uint i,j;
	bool flag = false;
	for (i = 0; i != halfedge_list->size(); ++i)
	{
		if (halfedge_list->at(i)->v1 == v2)
		{
			v3 = halfedge_list->at(i)->v2;
			for (j = 0; j != halfedge_list->size(); ++j)
			{
				if (halfedge_list->at(j)->v1 == v3 && halfedge_list->at(j)->v2 == v1)
				{
					flag = true;
					break;		
				}
			}
			if (flag)
			{
				break;
			}
		}
	}
	if (flag)
	{
		halfedge2_index = i;
		halfedge3_index = j;
		return true;
	}
	else
	{
		halfedge2_index = 0;
		halfedge3_index = 0;
		return false;
	}
}

HalfEdge* ACVT::find_opposite_halfedge(HalfEdge *halfedge_input, vector<vector<HalfEdge*>*>* halfedges_all)
{
	VertexID v1 = halfedge_input->v1, v2 = halfedge_input->v2;
	for (uint j = 0; j != halfedges_all->at(v2)->size(); ++j)
	{
		if (halfedges_all->at(v2)->at(j)->v2 == v1)
		{	// backward half edge 3: v1->v3
			return halfedges_all->at(v2)->at(j);
		}
	}
	return NULL;
}

HalfEdge* ACVT::get_half_edge(VertexID v1, VertexID v2, vector<vector<HalfEdge*>*>* halfedges_all)
{
	for (uint j = 0; j != halfedges_all->at(v1)->size(); ++j)
	{
		if (halfedges_all->at(v1)->at(j)->v2 == v2)
		{	// backward half edge 3: v1->v3
			return halfedges_all->at(v1)->at(j);
		}
	}
	return NULL;
}

void ACVT::reverse_current_normals()
{
	for (uint i = 0; i != vertex_number; ++i)
	{
		float *v = normals_current->at(i);
		for (uint j = 0; j != 3; ++j)
		{
			v[j] = -v[j];
		}
	}
}

void ACVT::assign_cluster_colors()
{
	/////////////////////////////////////////////////////////////////////
	// Method 1: Set cluster color for each cluster using four-color-problem method
	// Note: this code may still contain bugs.
	//set<int> clusterColors;
	//for (uint i = 0; i != vertex_number; ++i)
	//{
	//	VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
	//	if (!vtxQuad->cluster->empty())
	//	{
	//		clusterColors.clear();

	//		// For cluster i, save the cluster colors of all its neighbors
	//		for (set<VertexID>::iterator iter = neighbor_clusters->at(i).begin(); 
	//			iter != neighbor_clusters->at(i).end(); ++iter)
	//		{
	//			VertexQuadrics *vtxQuad1 = vtxQuadricsList->at(*iter);
	//			clusterColors.insert(vtxQuad1->get_cluster_color_no());
	//		}

	//		// Set a color for cluster i
	//		for (uint j = 0; j != color_number; ++j)
	//		{
	//			bool flag = true;
	//			for (set<int>::iterator iter1 = clusterColors.begin(); 
	//				iter1 != clusterColors.end(); ++iter1)
	//			{
	//				if (j == *iter1)
	//				{
	//					flag = false;
	//					break;
	//				}
	//			}
	//			if (flag)
	//			{
	//				// If this color j is different from all of cluster i's neighbors,
	//				// then use j as the color of cluster i
	//				vtxQuad->set_cluster_color_no(j);
	//				break;
	//			}
	//		}
	//	}
	//}

	//clusterColors.clear();
	//ClearSet(clusterColors);


	/////////////////////////////////////////////////////////////////////
	// Method 2: Set random color for each cluster
	srand ((unsigned)time(0));
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			if (!vtxQuad->is_cluster_color_set())
			{
				float f[3];
				for (uint j = 0; j != 3; ++j)
				{
					f[j] = rand()/float(RAND_MAX);
				}
				vtxQuad->set_cluster_color(f);
				vtxQuad->set_cluster_color_flag(true);
			}
		}
	}

}

uint ACVT::get_current_vertex_number()
{
	uint count = 0;
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			count++;
		}
	}
	return count;
}

void ACVT::split_clusters()
{
	//////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Firstly get all sub-clusters by splitting each cluster into different connected components
	//
	ALGraph *algraph = new ALGraph(vertex_number);
	vector<vector<VertexID>> subclusters;	// save all connected components of sub-clusters
	vector<VertexID> parentclusters;	// the parent clusters of all subclusters
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			// Firstly get all vertices of this cluster
			vector<VertexID> *clusterNodes = new vector<VertexID>;
			for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
				iter != vtxQuad->cluster->end(); ++iter)
			{
				clusterNodes->push_back(*iter);
			}
			// Create a graph based on all original edges of these vertices
			for (uint j = 0; j != clusterNodes->size(); ++j)
			{
				for (uint k = j+1; k != clusterNodes->size(); ++k)
				{
					VertexID v1 = clusterNodes->at(j);
					VertexID v2 = clusterNodes->at(k);
					sort2(v1,v2);
					EdgePairType::iterator iter = edgePairs.find(pair<VertexID,VertexID>(v1,v2));
					if (iter != edgePairs.end())
					{
						// Edge (v1,v2) exists in the original connectivity
						algraph->addEdge(v1,v2);
						algraph->addEdge(v2,v1);
					}
				}
			}
			// Get all sub-clusters for this cluster by checking its connected component
			bool flag = true;
			for (uint j = 0; j != clusterNodes->size(); ++j)
			{
				// Firstly run DFS on each vertex
				//algraph->DFS(clusterNodes->at(j));
				algraph->DFS_NonRecursive(clusterNodes->at(j));
				
				if (j == 0)
				{
					if (algraph->get_node_number_in_search() == vtxQuad->cluster->size())
					{
						// This cluster is a connected component, then do not split it.
						flag = false;
						break;
					}	
				}		
				if (algraph->get_node_number_in_search() != 0)
				{
					// The vertex number of DFS tree is smaller than the total number in original cluster,
					// then this tree is a sub-cluster.

					vector<VertexID> vecContComp;
					for (uint k = 0; k != algraph->vecVtxTraverse->size(); ++k)
					{
						vecContComp.push_back(algraph->vecVtxTraverse->at(k));
					}
					subclusters.push_back(vecContComp);
				}
				else
				{
					if (!algraph->isVisited(clusterNodes->at(j)))
					{	
						// This means that the vertex now is a floating point in the cluster, 
						// then treat it as an isolated island.
						vector<VertexID> vecContComp;
						vecContComp.push_back(clusterNodes->at(j));
						subclusters.push_back(vecContComp);
					}
				}
				algraph->vecVtxTraverse->clear();
				algraph->clear_node_number();
			}
			if (flag)
			{// Save the cluster to be split
				parentclusters.push_back(i);
			}			

			algraph->clearAll();
			algraph->init();
			clusterNodes->clear();
			ClearVector(*clusterNodes);
		}
	}
	delete algraph;

	//////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Create new clusters to hold these sub-clusters
	//
	// Firstly delete the parent clusters which need to be split
	for (uint i = 0; i != parentclusters.size(); ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(parentclusters[i]);
		if (!vtxQuad->cluster->empty())
		{
			vtxQuad->clearAll();
			valid_verts--;
		}
	}
	cout<<"Clusters to be split: "<<parentclusters.size()<<endl;
	cout<<"Sub-clusters found: "<<subclusters.size()<<endl;

	// Next, create a new cluster for each sub-cluster
	if (parentclusters.size() != 0 && subclusters.size() != 0)
	{
		// get the empty cluster list
		vector<VertexID> cluster_list;
		get_empty_cluster_list(cluster_list);

		if (cluster_list.size() >= subclusters.size())
		{
			// Create new clusters based on the sub-cluster list and empty cluster list
			create_new_clusters(cluster_list, subclusters);
		}
		else
		{
			cout<<"Does not have empty cluster positions to create new clusters!"<<endl;
		}

		// Clear related variables
		cluster_list.clear();
		ClearVector(cluster_list);
	}

	for (uint j = 0; j != subclusters.size(); ++j)
	{
		subclusters[j].clear();
	}
	subclusters.clear();
	ClearVector(subclusters);
	parentclusters.clear();
	ClearVector(parentclusters);
}


void ACVT::merge_all_clusters_inside_another()
{
	compute_cluster_belong();
	uint count = 0;
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			VertexID cluster_id = i;
			set<VertexID> neighbors;
			compute_neighbor_cluster(cluster_id, neighbors);
			if (!neighbors.empty())
			{
				if (neighbors.size() == 1)
				{
					// This sub-cluster only has 1 neighbor cluster, so it is inside this neighbor cluster.
					// Just merge it into its neighbor cluster.
					set<VertexID>::iterator iter0 = neighbors.begin();
					merge_cluster(cluster_id, *iter0);
					valid_verts--;
					//cout<<"Cluster "<<cluster_id<<" will be merged into cluster "<<*iter0<<endl;
				}
			}
			neighbors.clear();
		}
	}
}

void ACVT::merge_cluster(VertexID from_cluster_id, VertexID to_cluster_id)
{
	VertexQuadrics *vtxFrom = vtxQuadricsList->at(from_cluster_id);
	VertexQuadrics *vtxTo = vtxQuadricsList->at(to_cluster_id);
	for (set<VertexID>::iterator iter = vtxFrom->cluster->begin(); 
		iter != vtxFrom->cluster->end(); ++iter)
	{
		vtxTo->add_vertex_to_cluster(*iter);
		Quadrics(to_cluster_id) += InitialQuadrics(*iter);
	}
	compute_target(to_cluster_id);
	vtxFrom->clearAll();
}

bool ACVT::check_cluster_inside_another(vector<VertexID> &vecTemp)
{
	// Get all neighbor clusters of each vertex in this sub-cluster
	if (!vecTemp.empty())
	{
		set<VertexID> neighbor_cluster_IDs;
		for (uint j = 0; j != vecTemp.size(); ++j)
		{
			VertexID vid = vecTemp[j];
			VertexQuadrics *vid_Quad = vtxQuadricsList->at(vid);
			for (set<VertexID>::iterator iter = vid_Quad->neighbors->begin();
				iter != vid_Quad->neighbors->end(); ++iter)
			{
				neighbor_cluster_IDs.insert(cluster_belong->at(*iter));
			}
		}
		VertexID cluster_id = cluster_belong->at(vecTemp[0]);
		if (vecTemp.size() == 1 && neighbor_cluster_IDs.size() == 1)
		{
			// This sub-cluster has only 1 neighbor cluster, then this is a floating vertex
			// inside another large cluster.
			set<VertexID>::iterator iter = neighbor_cluster_IDs.begin();
			if (cluster_id != *iter)
			{
				return true;
			}
		}
		if (neighbor_cluster_IDs.size() == 2)
		{
			neighbor_cluster_IDs.erase(cluster_id);
			if (neighbor_cluster_IDs.size() == 1)
			{
				return true;
			}
		}
		neighbor_cluster_IDs.clear();
		ClearSet(neighbor_cluster_IDs);
		return false;
	}
	else
	{
		return false;
	}	
}

void ACVT::get_empty_cluster_list(vector<VertexID> &cluster_list)
{
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (vtxQuad->cluster->empty())
			cluster_list.push_back(i);
	}
}

void ACVT::create_new_clusters(vector<VertexID> &empty_cluster_list, vector<vector<VertexID>> &newclusters)
{
	if (empty_cluster_list.size() > newclusters.size())
	{
		int i = 0;
		for (vector<VertexID>::iterator iter = empty_cluster_list.begin(); 
			iter != empty_cluster_list.end(); ++iter)
		{
			VertexID clusterId = *iter;
			VertexQuadrics *vtxQuad = vtxQuadricsList->at(clusterId);
			for (uint j = 0; j != newclusters[i].size(); ++j)
			{
				vtxQuad->add_vertex_to_cluster(newclusters[i][j]);
				Quadrics(clusterId) += InitialQuadrics(newclusters[i][j]);
			}
			compute_target(clusterId);
			//if (newclusters.size() < 100)
			//{
			//	cout<<"Newly created clusters:"<<clusterId<<endl;
			//}
			i++;
			valid_verts++;
			vtxQuad->set_cluster_color_flag(false);
			if (i == newclusters.size())
			{
				break;
			}
		}
	}
	else
	{
		cout<<"No empty cluster positions for creating new clusters!"<<endl;
	}
}

void ACVT::compute_target(VertexID i)
{
	MxQuadric Q = Quadrics(i);
	if (!Q.optimize(TargetVertex(i)))
	{
		// if Q is invertible, use the mass center as target vertex
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		VectorType v(dim());
		for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
			iter != vtxQuad->cluster->end(); ++iter)
		{
			VectorType v_temp(dim());
			getVertexWithExtraData(*iter, v_temp);
			v += v_temp;
		}
		v /= vtxQuad->cluster->size();
		TargetVertex(i) = v;
		//cout<<"Warning: the quadrics for cluster "<<i<<" is invertible"<<endl;
	}	
}

void ACVT::compute_new_edgelist_heap()
{
	// Clear related variables
	// clear heap
	uint heapSize = heap.size();
	for (uint i = 0; i != heapSize; ++i)
	{
		EdgeInfo *info = (EdgeInfo*) heap.extract();
	}
	assert(heap.size() == 0);
	// clear global edge list
	if (!edgeList->empty())
	{
		for (vector<EdgeInfo*>::iterator iter = edgeList->begin();
			iter != edgeList->end(); ++iter)
		{
			delete *iter;
		}
	}
	edgeList->clear();
	ClearVector(*edgeList);
	// clear edges and faces for each vertex
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		vtxQuad->clear_all_edges();
		vtxQuad->clear_all_faces();
	}

	// Compute new edgelist and add all related edges into the maximum heap
	edgeList = new vector<EdgeInfo*>;
	
	// Compute new neighbor clusters and edge pairs
	compute_all_neighbor_clusters();
	compute_new_edge_pairs();

	// Create new edge list
	uint count = 0;
	for (map<pair<uint,uint>, uint>::iterator iter = edgePairs_new.begin(); 
		iter != edgePairs_new.end(); ++iter)
	{
		// Create this edge
		EdgeInfo *e = new EdgeInfo((*iter).first.first, (*iter).first.second, count);
		count++;

		// Compute the edge information including the quadrics error for this edge
		compute_edge_info(e);

		// Put the edge key inside the maximum heap as the negative quadrics error
		// Note: The QEM needs to get the smallest error each time, but since we use a maximum heap
		// here for this pop out operation, so here we use the negative error in the heap as the key.
		e->heap_key(e->quadError);

		// Save the edge
		edgeList->push_back(e);
		
		// Push it into the heap
		heap.insert(e);

		// Save this edge into the edgelist of related vertex
		vtxQuadricsList->at(e->v1)->addEdge(e);
		vtxQuadricsList->at(e->v2)->addEdge(e);
	}
}

void ACVT::save_cluster_in_file(char *filename)
{
	ofstream writeout;
	writeout.open(filename, ios::trunc);
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		writeout<<vtxQuad->cluster->size();
		if (!vtxQuad->cluster->empty())
		{
			for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
				iter != vtxQuad->cluster->end(); ++iter)
			{
				writeout<<" "<<*iter;
			}
		}
		writeout<<endl;
	}
	writeout.close();
}

void ACVT::open_cluster_file(char *filename, vector<vector<VertexID>*> *cluster_all)
{
	ifstream readin;
	readin.open(filename, ios::in);
	VertexID num = 0, val = 0;
	for (uint i = 0; i != vertex_number; ++i)
	{
		vector<VertexID> *vecTemp = new vector<VertexID>;
		readin>>num;
		if (num != 0)
		{
			for (uint j = 0; j != num; ++j)
			{
				readin>>val;
				vecTemp->push_back(val);
			}
		}
		cluster_all->push_back(vecTemp);
	}
	readin.close();
}

void ACVT::set_clusters_from_input(vector<vector<VertexID>*> *cluster_all)
{
	for (uint i = 0; i != cluster_all->size(); ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		vtxQuad->clearAll();
		vector<VertexID>* vecTemp = cluster_all->at(i);
		if (vecTemp->size() != 0)
		{
			for (uint j = 0; j != vecTemp->size(); ++j)
			{
				vtxQuad->add_vertex_to_cluster(vecTemp->at(j));
				Quadrics(i) += InitialQuadrics(vecTemp->at(j));
			}
			compute_target(i);
		}
	}
}

bool ACVT::split_nonmanifold_clusters()
{
	set<VertexID> clusters_nonmanifold;
	check_nonmanifold_clusters(clusters_nonmanifold);
	cout<<"\tNon-manifold clusters: "<<clusters_nonmanifold.size()<<endl;
	bool flag = true;
	if (clusters_nonmanifold.size() > 0)
	{
		vector<vector<VertexID>> subclusters;
		uint count = 0;// Save the number of clusters to be split
		for (set<VertexID>::iterator iter = clusters_nonmanifold.begin();
			iter != clusters_nonmanifold.end(); ++iter)
		{
			VertexQuadrics *vtxQuad = vtxQuadricsList->at(*iter);
			// Only split clusters with more than 1 vertex
			if (vtxQuad->cluster->size() != 1)
			{
				count++;
				// Split the cluster and then delete it
				split_cluster_into_two(*iter, subclusters);
				//if (clusters_nonmanifold.size() < 50)
				//{
				//	cout<<"clusters to be split:"<<*iter<<endl;
				//}
				vtxQuad->clearAll();
				valid_verts--;
			}
		}
		if (count != 0)
		{
			cout<<"\tClusters to be split: "<<count<<endl;
			cout<<"\tSub-clusters to be created: "<<subclusters.size()<<endl;
			//if (subclusters.size() != 2*count)
			//{
			//	cout<<"Warning: there must be at least one cluster not split well."<<endl;
			//}
			//if (count < 50)
			//{
			//	cout<<"clusters to be split:"<<endl;
			//}
			vector<VertexID> empty_cluster_list;
			get_empty_cluster_list(empty_cluster_list);
			create_new_clusters(empty_cluster_list, subclusters);
			empty_cluster_list.clear();
			ClearVector(empty_cluster_list);
		}
		else
		{
			/* 
			* This is the case that there are still non-manifold clusters, but each of
			* these clusters only contains 1 vertex, since the original input model contains 
			* non-manifold edges.
			*/
			cout<<"\tClusters to be split: 0"<<endl;
			flag = false;
		}
		for (uint i = 0; i != subclusters.size(); ++i)
		{
			subclusters[i].clear();
			ClearVector(subclusters[i]);
		}
		subclusters.clear();
		ClearVector(subclusters);
	}
	else
	{
		// All clusters are manifold
		flag = false;
	}
	clusters_nonmanifold.clear();
	ClearSet(clusters_nonmanifold);
	return flag;
}

void ACVT::check_nonmanifold_clusters(set<VertexID>& cluster_ids)
{
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *v0_Quad = vtxQuadricsList->at(i);
		if (!v0_Quad->cluster->empty())
		{
			for (set<VertexID>::iterator iter = neighbor_clusters->at(i).begin();
				iter != neighbor_clusters->at(i).end(); ++iter)
			{				
				if (cluster_ids.find(i) != cluster_ids.end() && cluster_ids.find(*iter) != cluster_ids.end())
				{
					// These two clusters has already been checked, continue
					continue;
				}

				vector<VertexID> common_neighbors;
				uint a = neighbor_clusters->at(i).size();// number of neighbors of cluster i
				uint b = neighbor_clusters->at(*iter).size();// number of neighbors of (*iter), one of i's neighbors
				// The number of common neighbors cannot be larger than max(a,b)
				uint len = (a > b) ? a : b;
				common_neighbors.resize(len);
				// Get the intersection of i's and (*iter)'s neighbors to get their common neighbors
				// Note: here "it" denotes the pointer to the last position of common elements in "common_neighbors" 
				vector<VertexID>::iterator it = set_intersection(neighbor_clusters->at(i).begin(), neighbor_clusters->at(i).end(),
					neighbor_clusters->at(*iter).begin(), neighbor_clusters->at(*iter).end(), common_neighbors.begin());
				common_neighbors.resize(it-common_neighbors.begin());// resize "common_neighbors" only to keep the common neighbors
				if (common_neighbors.size() > 2)
				{
					// These two clusters have more than 2 common neighbors, then they denote a non-manifold edge.
					cluster_ids.insert(i);
					cluster_ids.insert(*iter);
					// NOTE: here the common neighbor clusters of the two non manifold clusters are also saved and to be split.
					// If you do not want this, please comment the following 3 lines.
					//for (uint j = 0; j != common_neighbors.size(); ++j)
					//{
					//	cluster_ids.insert(common_neighbors[j]);
					//}
				}
				common_neighbors.clear();
				ClearVector(common_neighbors);
			}
		}
	}
}

void ACVT::split_cluster_into_two(VertexID cluster_id, vector<vector<VertexID>> &subclusters)
{
	VertexQuadrics *vtxQuad = vtxQuadricsList->at(cluster_id);
	if (!vtxQuad->cluster->empty())
	{
		/* Compute center and first principal axis of this cluster first */
		DualQEM Q_cor;
		//Vector3d v_center;
		for (set<VertexID>::iterator iter = vtxQuad->cluster->begin();
			iter != vtxQuad->cluster->end(); ++iter)
		{
			VectorType v(dim());
			getVertexWithExtraData(*iter,v);
			DualQEM Q(v[0], v[1], v[2]);
			//v_center[0] += v[0];
			//v_center[1] += v[1];
			//v_center[2] += v[2];
			Q_cor += Q;
		}
		Vector3d v_axis;
		Q_cor.optimize(v_axis);// the first PCA principal axis
		Vector3d center_cluster = Q_cor.center;
		//v_center[0] /= vtxQuad->cluster->size();
		//v_center[1] /= vtxQuad->cluster->size();
		//v_center[2] /= vtxQuad->cluster->size();
		//Vector3d center_cluster = v_center;
		/* 
		* Split the cluster into two sub-clusters like this: let v be one point
		* in this cluster, e be the PCA principal axis, c be the center of this
		* cluster. If the angle between v-c and e is smaller than 90 degree,
		* then v belongs to the 1st sub-cluster; otherwise, v belongs to 
		* the 2nd sub-cluster.
		*/
		vector<VertexID> cluster1, cluster2;
		for (set<VertexID>::iterator iter = vtxQuad->cluster->begin();
			iter != vtxQuad->cluster->end(); ++iter)
		{
			VectorType v(dim());
			getVertexWithExtraData(*iter,v);
			Vector3d vtx;
			for (uint j = 0; j != 3; ++j)
			{
				vtx[j] = v[j] - center_cluster[j];
			}
			// If dot product between v-c and e is larger than 0, then v belongs to one cluster;
			// Otherwise, v belongs to another one.
			double dot_result = vtx[0]*v_axis[0]+vtx[1]*v_axis[1]+vtx[2]*v_axis[2];
			if (dot_result < 1e-10)
			{
				cluster1.push_back(*iter);
			}
			else
			{
				cluster2.push_back(*iter);
			}
		}
		if (!cluster1.empty() && !cluster2.empty())
		{
			subclusters.push_back(cluster1);
			subclusters.push_back(cluster2);
		}
		else
		{
			/*
			* Here the cluster cannot be split into two sub-clusters successfully. The 
			* most possible reason is that all the points in this cluster have almost "exactly" 
			* the same coordinates. Then, the PCA principal axis and the center coordinates
			* are not very accurate because of numerical error. Here I just split the cluster
			* equally based on the value of each point's "dot_result".
			*/
			//if (vtxQuad->cluster->size() == 2 || vtxQuad->cluster->size() == 3)
			//{
			//	cout<<"Cluster "<<cluster_id<<" has "<<vtxQuad->cluster->size()<<" points"<<endl;
			//	uint count_v = 0;
			//	for (set<VertexID>::iterator iter = vtxQuad->cluster->begin();
			//		iter != vtxQuad->cluster->end(); ++iter)
			//	{
			//		VectorType v(dim());
			//		getVertexWithExtraData(*iter,v);
			//		cout<<"v"<<count_v++<<":"<<v[0]<<","<<v[1]<<","<<v[2]<<endl;
			//	}
			//	cout<<"Center:"<<center_cluster[0]<<","<<center_cluster[1]<<","<<center_cluster[2]<<endl;
			//	cout<<"Principal Axis:"<<v_axis[0]<<","<<v_axis[1]<<","<<v_axis[2]<<endl;
			//}
			/* Get dot product (v-c)*e for each point v */
			cluster1.clear();
			cluster2.clear();
			set<pair<double,uint>> dot_results;
			for (set<VertexID>::iterator iter = vtxQuad->cluster->begin();
				iter != vtxQuad->cluster->end(); ++iter)
			{
				VectorType v(dim());
				getVertexWithExtraData(*iter,v);
				Vector3d vtx;
				for (uint j = 0; j != 3; ++j)
				{
					vtx[j] = v[j] - center_cluster[j];
				}
				double temp_val = vtx[0]*v_axis[0]+vtx[1]*v_axis[1]+vtx[2]*v_axis[2];
				dot_results.insert(make_pair(temp_val,*iter));
			}
			/* Split this cluster into two sub-clusters equally based on dot products*/
			uint size_cluster1 = vtxQuad->cluster->size()/2;
			uint count = 0;
			for (set<pair<double,uint>>::iterator iter = dot_results.begin();
				iter != dot_results.end(); ++iter)
			{
				if (count < size_cluster1)
				{
					cluster1.push_back((*iter).second);
				}
				else
				{
					cluster2.push_back((*iter).second);
				}
				count++;
			}
			if (cluster1.empty() || cluster2.empty())
			{
				cout<<"Warning: cluster "<<cluster_id<<" cannot be split."<<endl;
			}
			subclusters.push_back(cluster1);
			subclusters.push_back(cluster2);
			dot_results.clear();
			ClearSet(dot_results);
		}
		uint size_total = vtxQuad->cluster->size();
		uint size1 = cluster1.size(), size2 = cluster2.size();
		if (size_total != size1+size2)
		{
			cout<<"Warning: cluster "<<cluster_id<<" is not split very well"<<endl;
		}
		if (cluster_id == 75967)
		{
			uint count_v = 0;
			for (set<VertexID>::iterator iter = vtxQuad->cluster->begin();
				iter != vtxQuad->cluster->end(); ++iter)
			{
				VectorType v(dim());
				getVertexWithExtraData(*iter,v);
				cout<<"v"<<count_v++<<":"<<v[0]<<","<<v[1]<<","<<v[2]<<endl;
			}
			cout<<"Center:"<<center_cluster[0]<<","<<center_cluster[1]<<","<<center_cluster[2]<<endl;
			cout<<"Principal Axis:"<<v_axis[0]<<","<<v_axis[1]<<","<<v_axis[2]<<endl;
			for (uint j = 0; j != cluster1.size(); ++j)
			{
				cout<<"Sub-Cluster 1:"<<cluster1[j]<<endl;
			}
			for (uint j = 0; j != cluster2.size(); ++j)
			{
				cout<<"Sub-Cluster 2:"<<cluster2[j]<<endl;
			}
		}
	}
}

void ACVT::compute_final_vertex_quality()
{
	for (uint i = 0; i != vertex_number; ++i)
	{
		//VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		//vtxQuad->initialQuadrics->clear(0);
		InitialQuadrics(i).clear(0);
	}

	for (uint i = 0; i != face_number; ++i)
	{
		int *t = mesh->getTriangle(i);
		for (uint j = 0; j != 3; ++j)
		{
			VectorType v0(dim());
			VectorType v1(dim());
			VectorType v2(dim());
			v0[0] = mesh->getVtxExtraData(t[j],0);
			v1[0] = mesh->getVtxExtraData(t[(j+1)%3],0);
			v2[0] = mesh->getVtxExtraData(t[(j+2)%3],0);
			// The forth parameter r is a weight constraint to the face
			double r = 1.0;
			MxQuadric Q_p(v0, v1, v2, r);
			InitialQuadrics(t[j]) += Q_p;
		}
	}
	
	for (uint i = 0; i != vertex_number; ++i)
	{
		VertexQuadrics *vtxQuad = vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			MxQuadric Q(dim());
			for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
				iter != vtxQuad->cluster->end(); ++iter)
			{
				Q += InitialQuadrics(*iter);
			}
			VectorType v_quality(dim());
			if (Q.optimize(v_quality))
			{
				vtxQuad->vertexQuality = Q(v_quality);
			}
			else
			{
				// Use the average quality value of this cluster as the final
				// vertex quality value of this cluster vertex
				double sum = 0;
				for (set<VertexID>::iterator iter = vtxQuad->cluster->begin(); 
					iter != vtxQuad->cluster->end(); ++iter)
				{
					sum += mesh->getVtxExtraData(*iter);
				}
				vtxQuad->vertexQuality = sum/vtxQuad->cluster->size();
			}			
		}
	}
}