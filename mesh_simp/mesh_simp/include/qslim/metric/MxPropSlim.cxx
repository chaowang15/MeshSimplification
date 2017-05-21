/************************************************************************

  MxPropSlim

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxPropSlim.cxx,v 1.9 2000/11/20 20:36:38 garland Exp $

 ************************************************************************/

#include "stdmix.h"
#include "MxPropSlim.h"
#include "MxGeom3D.h"
#include "domain.h"
#include <set>
#include <vector>

extern domain * vertex_domains;
extern MxStdModel * init_model;
extern vector<float> * trajectory;
typedef MxQuadric Quadric;

MxPropSlim::MxPropSlim(MxStdModel *m0)
    : MxStdSlim(m0),
      __quadrics(m0->vert_count()),
	  initial_quadrics(m0->vert_count()),
      edge_links(m0->vert_count())
{
    consider_color();
    consider_texture();
    consider_normals();

    D = compute_dimension(m);

    will_decouple_quadrics = false;
}

void MxPropSlim::consider_color(bool will)
{
    use_color = will && (m->color_binding() == MX_PERVERTEX);
    D = compute_dimension(m);
}

void MxPropSlim::consider_texture(bool will)
{
    use_texture = will && (m->texcoord_binding() == MX_PERVERTEX);
    D = compute_dimension(m);
}

void MxPropSlim::consider_normals(bool will)
{
    use_normals = will && (m->normal_binding() == MX_PERVERTEX);
    D = compute_dimension(m);
}

uint MxPropSlim::compute_dimension(MxStdModel *m)
{
    uint d = 3;

    if( use_color )  d += 3;
    if( use_texture )  d += 2;
    if( use_normals )  d += 3;

    return d;
}

void MxPropSlim::pack_to_vector(MxVertexID id, MxVector& v)
{
    SanityCheck( v.dim() == D );
    SanityCheck( id < m->vert_count() );

    v[0] = m->vertex(id)[0];
    v[1] = m->vertex(id)[1];
    v[2] = m->vertex(id)[2];

    uint i = 3;
    if( use_color )
    {
	v[i++] = m->color(id).R();
	v[i++] = m->color(id).G();
	v[i++] = m->color(id).B();
    }
    if( use_texture )
    {
	v[i++] = m->texcoord(id)[0];
	v[i++] = m->texcoord(id)[1];
    }
    if( use_normals )
    {
	v[i++] = m->normal(id)[0];
	v[i++] = m->normal(id)[1];
	v[i++] = m->normal(id)[2];
    }
}

void MxPropSlim::pack_prop_to_vector(MxVertexID id, MxVector& v, uint target)
{
    if( target == 0 )
    {
	v[0] = m->vertex(id)[0];
	v[1] = m->vertex(id)[1];
	v[2] = m->vertex(id)[2];
	return;
    }

    uint i = 3;
    target--;

    if( use_color )
    {
	if( target == 0 )
	{
	    v[i]   = m->color(id).R();
	    v[i+1] = m->color(id).G();
	    v[i+2] = m->color(id).B();
	    return;
	}
	i += 3;
	target--;
    }
    if( use_texture )
    {
	if( target == 0 )
	{
	    v[i]   = m->texcoord(id)[0];
	    v[i+1] = m->texcoord(id)[1];
	    return;
	}
	i += 2;
	target--;
    }
    if( use_normals )
    {
	if( target == 0 )
	{
	    v[i]   = m->normal(id)[0];
	    v[i+1] = m->normal(id)[1];
	    v[i+2] = m->normal(id)[2];
	    return;
	}
    }
}

static inline void CLAMP(double& v, double lo, double hi)
{
    if( v<lo ) v = lo;
    if( v>hi ) v = hi;
}

void MxPropSlim::unpack_from_vector(MxVertexID id, MxVector& v)
{
    SanityCheck( v.dim() == D );
    SanityCheck( id < m->vert_count() );

    m->vertex(id)[0] = v[0];
    m->vertex(id)[1] = v[1];
    m->vertex(id)[2] = v[2];

    uint i = 3;
    if( use_color )
    {
	//CLAMP(v[i], 0, 1);
	//CLAMP(v[i+1], 0, 1);
	//CLAMP(v[i+2], 0, 1);
	m->color(id).set(v[i], v[i+1], v[i+2]);
	i += 3;
    }
    if( use_texture )
    {
	m->texcoord(id)[0] = v[i++];
        m->texcoord(id)[1] = v[i++];
    }
    if( use_normals )
    {
	float n[3];  n[0]=v[i++];  n[1]=v[i++];  n[2]=v[i++];
	mxv_unitize(n, 3);
	m->normal(id).set(n[0], n[1], n[2]);
    }
}

void MxPropSlim::unpack_prop_from_vector(MxVertexID id,MxVector& v,uint target)
{
    if( target == 0 )
    {
	m->vertex(id)[0] = v[0];
	m->vertex(id)[1] = v[1];
	m->vertex(id)[2] = v[2];
	return;
    }

    uint i=3;
    target--;

    if( use_color )
    {
	if( target == 0 )
	{
	    m->color(id).set(v[i], v[i+1], v[i+2]);
	    return;
	}
	i+=3;
	target--;
    }
    if( use_texture )
    {
	if( target == 0 )
	{
	    m->texcoord(id)[0] = v[i];
	    m->texcoord(id)[1] = v[i+1];
	    return;
	}
	i += 2;
	target--;
    }
    if( use_normals )
    {
	if( target == 0 )
	{
	    float n[3];  n[0]=v[i];  n[1]=v[i+1];  n[2]=v[i+2];
	    mxv_unitize(n, 3);
	    m->normal(id).set(n[0], n[1], n[2]);
	    return;
	}
    }
}


uint MxPropSlim::prop_count()
{
    uint i = 1;

    if( use_color ) i++;
    if( use_texture) i++;
    if( use_normals ) i++;

    return i;
}

void MxPropSlim::compute_face_quadric(MxFaceID i, MxQuadric& Q)
{
    MxFace& f = m->face(i);

    MxVector v1(dim());
    MxVector v2(dim());
    MxVector v3(dim());

    if( will_decouple_quadrics )
    {
	Q.clear();

	for(uint p=0; p<prop_count(); p++)
	{
	    v1=0.0;  v2=0.0;  v3=0.0;

	    pack_prop_to_vector(f[0], v1, p);
	    pack_prop_to_vector(f[1], v2, p);
	    pack_prop_to_vector(f[2], v3, p);

	    // !!BUG: Will count area multiple times (once per property)
	    MxQuadric Q_p(v1, v2, v3, m->compute_face_area(i));

	    // !!BUG: Need to only extract the relevant block of the matrix.
	    //        Adding the whole thing gives us extraneous stuff.
	    Q += Q_p;
	}
    }
    else
    {
	pack_to_vector(f[0], v1);
	pack_to_vector(f[1], v2);
	pack_to_vector(f[2], v3);

	Q = MxQuadric(v1, v2, v3, m->compute_face_area(i));
    }
}

void MxPropSlim::collect_quadrics()
{
    for(uint j=0; j<quadric_count(); j++)
	{
		__quadrics[j] = new MxQuadric(dim());
		initial_quadrics[j] = new MxQuadric(dim());
	}
	

    for(MxFaceID i=0; i<m->face_count(); i++)
    {
	MxFace& f = m->face(i);

	MxQuadric Q(dim());
	compute_face_quadric(i, Q);

// 	if( weight_by_area )
// 	    Q *= Q.area();

	quadric(f[0]) += Q;
	quadric(f[1]) += Q;
	quadric(f[2]) += Q;
    }
}

void MxPropSlim::initialize()
{
	vertex_domains = new domain(valid_verts);
	//neighbor_node = new set<int>[valid_verts];
	//for (int i=0; i<valid_verts; i++)
	//	neighbor_node[i].clear();
	
    collect_quadrics();

    if( boundary_weight > 0.0 )
 	//constrain_boundaries();

	//make a copy of the quadrics
	for(uint j=0; j<quadric_count(); j++)
	{
		//initial_quadrics[j] = __quadrics[j];
		init_quadric(j) = quadric(j);
	}

    collect_edges();

    is_initialized = true;
}

void MxPropSlim::compute_target_placement(edge_info *info)
{
    MxVertexID i=info->v1, j=info->v2;

    const MxQuadric &Qi=quadric(i), &Qj=quadric(j);
    MxQuadric Q=Qi;  Q+=Qj;

    double err;

    if( Q.optimize(info->target) )
    {
	err = Q(info->target);
    }
    else
    {
	// Fall back only on endpoints

	MxVector v_i(dim()), v_j(dim());

	pack_to_vector(i, v_i);
	pack_to_vector(j, v_j);

	double e_i = Q(v_i);
	double e_j = Q(v_j);

	if( e_i<=e_j )
	{
	    info->target = v_i;
	    err = e_i;
	}
	else
	{
	    info->target = v_j;
	    err = e_j;
	}
    }
	//if (err < 0) err += 1e200; // The punishment of negative error, make it a great positive error 
	//err += 1e-10 * Q.area(); // Make sure that area works even the error is zero
//     if( weight_by_area )
// 	err / Q.area();
    info->heap_key(-err);
}

bool MxPropSlim::decimate(uint target)
{
    MxPairContraction conx;

	while( valid_verts > target )
	{
	edge_info *info = (edge_info *)heap.extract();
	if( !info )  return false;

	MxVertexID v1=info->v1, v2=info->v2;

	if( m->vertex_is_valid(v1) && m->vertex_is_valid(v2) )
	{
	    m->compute_contraction(v1, v2, &conx);

	    conx.dv1[X] = info->target[X] - m->vertex(v1)[X];
	    conx.dv1[Y] = info->target[Y] - m->vertex(v1)[Y];
	    conx.dv1[Z] = info->target[Z] - m->vertex(v1)[Z];
	    conx.dv2[X] = info->target[X] - m->vertex(v2)[X];
	    conx.dv2[Y] = info->target[Y] - m->vertex(v2)[Y];
	    conx.dv2[Z] = info->target[Z] - m->vertex(v2)[Z];

	    apply_contraction(conx, info);
	}	
	delete info;
    }
    return true;
}

bool MxPropSlim::swap_per_pixel()
{
	int swap_count = 0;
	while(pixel_heap.size() > 0)
	{
		pixel_info * info = (pixel_info *)pixel_heap.extract();
		if( !info )  return false;

		apply_swap_per_pixel(info);
		swap_count++;

		delete info;
	}
	cout << "pixel swap done! swap " << swap_count << " times.";
	return true;	
}

bool MxPropSlim::swap_pixel_set()
{
	int swap_count = 0;
	int iteration_count = 0;
	set<pixel_info *> info_set;
	info_set.clear();
	while(pixel_heap.size() > 0)
	{
		for (int i=0; i<pixel_heap.size() ; i++)
		{
			pixel_info * info = (pixel_info *)pixel_heap.extract();
			assert(info != NULL);
			info_set.insert(info);
		}
		apply_swap_set(info_set);
		for (set<pixel_info *>::iterator it = info_set.begin(); it != info_set.end(); it++)
			delete *it;
		info_set.clear();

		iteration_count++;
		swap_count += pixel_heap.size() ;
	}
	cout << "pixel swap done! swap " << iteration_count << " iterations with " << swap_count << " times";
	return true;	
}

bool MxPropSlim::swap_pixel_lloyd()
{	
	int swap_count = 0;
	int iteration_count = 0;
	set<pixel_info *> info_set;
	info_set.clear();
	while(pixel_heap.size() > 0)
	{
		for (int i=0; i<pixel_heap.size() ; i++)
		{
			pixel_info * info = (pixel_info *)pixel_heap.extract();
			assert(info != NULL);
			info_set.insert(info);
		}
		lloyd_swap_set(info_set);
		for (set<pixel_info *>::iterator it = info_set.begin(); it != info_set.end(); it++)
			delete *it;
		info_set.clear();

		iteration_count++;
		swap_count += pixel_heap.size() ;
	}
	cout << "lloyd done! swap " << iteration_count << " iterations with " << swap_count << " times";
	return true;	
}

bool MxPropSlim::swap_pixel_lloyd_once()
{	
	set<pixel_info *> info_set;
	info_set.clear();
	if(pixel_heap.size() > 0)
	{
		for (int i=0; i<pixel_heap.size() ; i++)
		{
			pixel_info * info = (pixel_info *)pixel_heap.extract();
			assert(info != NULL);
			info_set.insert(info);
		}
		lloyd_swap_set(info_set);
		for (set<pixel_info *>::iterator it = info_set.begin(); it != info_set.end(); it++)
			delete *it;
		info_set.clear();
		cout << "perform lloyd swap one iteration" << endl;;
		return true;	
	}
	else
		return false;
}

bool MxPropSlim::swap_per_pixel_once()
{
	if(pixel_heap.size() > 0)
	{
		pixel_info * info = (pixel_info *)pixel_heap.extract();
		if( !info )  return false;

		apply_swap_per_pixel(info);

		delete info;
		return true;
	}
	else
		return false;
		
}

bool MxPropSlim::swap_pixel_set_once()
{
	set<pixel_info *> info_set;
	info_set.clear();
	if(pixel_heap.size() > 0)
	{
		for (int i=0; i<pixel_heap.size() ; i++)
		{
			pixel_info * info = (pixel_info *)pixel_heap.extract();
			assert(info != NULL);
			info_set.insert(info);
		}
		apply_swap_set(info_set);
		for (set<pixel_info *>::iterator it = info_set.begin(); it != info_set.end(); it++)
			delete *it;
		info_set.clear();

		cout << "perform pixel swap one iteration" << endl;;
		return true;
	}
	else
		return false;
}

void MxPropSlim::decimateOnce()
{
	MxPairContraction conx;

	edge_info *info = (edge_info *)heap.extract();
	if( !info )  return;

	MxVertexID v1=info->v1, v2=info->v2;

	if( m->vertex_is_valid(v1) && m->vertex_is_valid(v2) )
	{
	    m->compute_contraction(v1, v2, &conx);

	    conx.dv1[X] = info->target[X] - m->vertex(v1)[X];
	    conx.dv1[Y] = info->target[Y] - m->vertex(v1)[Y];
	    conx.dv1[Z] = info->target[Z] - m->vertex(v1)[Z];
	    conx.dv2[X] = info->target[X] - m->vertex(v2)[X];
	    conx.dv2[Y] = info->target[Y] - m->vertex(v2)[Y];
	    conx.dv2[Z] = info->target[Z] - m->vertex(v2)[Z];

	    apply_contraction(conx, info);
	}

	delete info;
}

void MxPropSlim::getTopQslimEdge(int &v1, int &v2)
 {
	 edge_info *info = (edge_info *)heap.top();
	 v1=info->v1;
	 v2=info->v2;
 }

 bool MxPropSlim::getTopPswapPixel(int &from, int &to, int &id)
 {
	 if (pixel_heap.size())
	 {
		 pixel_info * info = (pixel_info *)pixel_heap.top();
		 from = info->from_id;
		 to = info->to_id;
		 id = info->pixel_id;
		 return true;
	 }
	 else
		 return false;
	 
 }

 float MxPropSlim::vertex_dist(MxVertexID init_model_vert, MxVertexID contracted_model_vert)
 {
	 MxVertex init_vert = init_model->vertex(init_model_vert);
	 MxVertex contracted_vert = m->vertex(contracted_model_vert);
	 float dist =0;
	for (int j=0; j<3; j++)//coordinates
	{		 
		dist += (init_vert.as.pos[j] - contracted_vert.as.pos[j])*(init_vert.as.pos[j] - contracted_vert.as.pos[j]);
	}		 
	 
	 return dist;
 }

////////////////////////////////////////////////////////////////////////
//
// This is *very* close to the code in MxEdgeQSlim

void MxPropSlim::create_edge(MxVertexID i, MxVertexID j)
{
    edge_info *info = new edge_info(dim());

    edge_links(i).add(info);
    edge_links(j).add(info);

    info->v1 = i;
    info->v2 = j;

    compute_edge_info(info);
}

void MxPropSlim::discontinuity_constraint(MxVertexID i, MxVertexID j,
					  const MxFaceList& faces)
{
    for(uint f=0; f<faces.length(); f++)
    {
	Vec3 org(m->vertex(i)), dest(m->vertex(j));
	Vec3 e = dest - org;

	Vec3 v1(m->vertex(m->face(faces(f))(0)));
	Vec3 v2(m->vertex(m->face(faces(f))(1)));
	Vec3 v3(m->vertex(m->face(faces(f))(2)));
	Vec3 n = triangle_normal(v1,v2,v3);

	Vec3 n2 = e ^ n;
	unitize(n2);

	MxQuadric3 Q3(n2, -(n2*org));
	Q3 *= boundary_weight;

	MxQuadric Q(Q3, dim());

	quadric(i) += Q;
	quadric(j) += Q;
    }
}

void MxPropSlim::apply_contraction(const MxPairContraction& conx,
				   edge_info *info)
{
	vertex_domains->mergeDomain(conx.v1, conx.v2);
	
    valid_verts--;
    valid_faces -= conx.dead_faces.length();
    quadric(conx.v1) += quadric(conx.v2);

    update_pre_contract(conx);

    m->apply_contraction(conx);

    unpack_from_vector(conx.v1, info->target);

    // Must update edge_info here so that the meshing penalties
    // will be computed with respect to the new mesh rather than the old
    for(uint i=0; i<edge_links(conx.v1).length(); i++)
        compute_edge_info(edge_links(conx.v1)[i]);
}
void MxPropSlim::collect_update_domain(set<MxVertexID> & update_domain, pixel_info * info)
{
	update_domain.clear();
	update_domain.insert(info->from_id);
	update_domain.insert(info->to_id);

	MxVertexList star;
	star.reset();
	m->collect_vertex_star(info->from_id, star);
	for(uint j=0; j<star.length(); j++)
	{
		update_domain.insert(star(j));
	}
	star.reset();
	m->collect_vertex_star(info->to_id, star);
	for(uint j=0; j<star.length(); j++)
	{
		update_domain.insert(star(j));
	}
}

void MxPropSlim::apply_swap_per_pixel(pixel_info * info)
{
	cout << "key value : " << info->heap_key() << endl;
	assert(vertex_domains->vertex_domain[info->from_id].erase(info->pixel_id));
	vertex_domains->vertex_domain[info->to_id].insert(info->pixel_id);
	//update vertex position of v1 and v2 
	update_vertex_position(info);
	//delete link information
	assert(pixel_from_set[info->from_id].erase(info));
	assert(pixel_to_set[info->to_id].erase(info));
	pixel_id_pointer[info->pixel_id] = NULL;

	////aux
	//vertex_domains->red[info->pixel_id] = vertex_domains->red[info->to_id];
	//vertex_domains->green[info->pixel_id] = vertex_domains->green[info->to_id];
	//vertex_domains->blue[info->pixel_id] = vertex_domains->blue[info->to_id];

	//all the pixels belong to v1(v2) and its neighbor should update its information
	set<MxVertexID> update_domain;
	collect_update_domain(update_domain, info);

	for (set<MxVertexID>::iterator domain_it = update_domain.begin(); domain_it != update_domain.end(); domain_it++)
	{
		set<int>::iterator pixel_it;
		MxVertexID closest;
		float energy_saved_closest;
		MxVector from_update(dim());
		MxVector to_update(dim());

		//collect all the neighbor nodes
		MxVertexList star;
		star.reset();
		m->collect_vertex_star(*domain_it, star);

		for (pixel_it = vertex_domains->vertex_domain[*domain_it].begin(); pixel_it != vertex_domains->vertex_domain[*domain_it].end(); pixel_it++)
		{
			closest = *domain_it;
			energy_saved_closest=0;
			pixel_info *info = new pixel_info(dim());

			//compute the closest vertex
			for(uint j=0; j<star.length(); j++)
			{
				//float dist = vertex_dist(*pixel_it, star(j));
				float energy_saved = calc_energy_delta(*domain_it,star(j),*pixel_it, info);
				if ( energy_saved> energy_saved_closest)
				{
					from_update = info->from_temp_update;
					to_update = info->to_temp_update;
					closest = star(j);
					energy_saved_closest = energy_saved;
				}
			}
			if(energy_saved_closest >1e-5)
			{			
				if (pixel_id_pointer[*pixel_it] == NULL)
					//insert to heap
					heap_pixel_swap(energy_saved_closest,*domain_it,closest,*pixel_it, from_update,to_update, info);
				else
				{
					//modify the existing one
					pixel_info *exist_info = pixel_id_pointer[*pixel_it];
					assert(exist_info->is_in_heap());
					//from, id pointer should be the same
					assert(*domain_it == exist_info->from_id);
					pixel_to_set[exist_info->to_id].erase(exist_info);
					exist_info->to_id = closest;
					pixel_to_set[exist_info->to_id].insert(exist_info);
					exist_info->from_update = from_update;
					exist_info->to_update = to_update;
					exist_info->heap_key(energy_saved_closest);
					pixel_heap.update(exist_info);
					delete info;
				}
			}	
			else
			{
				//delete the existing one
				if (pixel_id_pointer[*pixel_it] != NULL)
				{
					pixel_info *exist_info = pixel_id_pointer[*pixel_it];				
					pixel_id_pointer[*pixel_it] = NULL;
					assert(pixel_from_set[exist_info->from_id].erase(exist_info));
					assert(pixel_to_set[exist_info->to_id].erase(exist_info));
					pixel_heap.remove(exist_info);
				}
				delete info;
			}
		}
	}
}

void MxPropSlim::apply_swap_set(set<pixel_info *> & info_set)
{
	set<MxVertexID> position_domain;// the vertex set that need position update
	position_domain.clear();
	for (set<pixel_info *>::iterator info_it = info_set.begin(); info_it != info_set.end(); info_it++)
	{
		pixel_info * info = *info_it;
		//cout << "key value : " << info->heap_key() << endl;
		assert(vertex_domains->vertex_domain[info->from_id].erase(info->pixel_id));
		vertex_domains->vertex_domain[info->to_id].insert(info->pixel_id);
		
		//delete link information
		assert(pixel_from_set[info->from_id].erase(info));
		assert(pixel_to_set[info->to_id].erase(info));
		pixel_id_pointer[info->pixel_id] = NULL;

		////aux
		//vertex_domains->red[info->pixel_id] = vertex_domains->red[info->to_id];
		//vertex_domains->green[info->pixel_id] = vertex_domains->green[info->to_id];
		//vertex_domains->blue[info->pixel_id] = vertex_domains->blue[info->to_id];

		//update quadric of v1 and v2 
		quadric(info->from_id) -= init_quadric(info->pixel_id);
		quadric(info->to_id) += init_quadric(info->pixel_id);

		//collect the vertexs that need position update
		position_domain.insert(info->from_id);
		position_domain.insert(info->to_id);
	}
	//update vertex position
	pixel_info * info = new pixel_info(dim());// for temp use, no meaning
	for(set<MxVertexID>::iterator vertex_it = position_domain.begin(); vertex_it != position_domain.end(); vertex_it++)
	{
		const MxQuadric &Q =quadric(*vertex_it);
		if(Q.optimize(info->to_update))
		{
			m->vertex(*vertex_it).as.pos[0] = info->to_update[0];
			m->vertex(*vertex_it).as.pos[1] = info->to_update[1];
			m->vertex(*vertex_it).as.pos[2] = info->to_update[2];

			trajectory[*vertex_it].push_back(m->vertex(*vertex_it).as.pos[0]);
			trajectory[*vertex_it].push_back(m->vertex(*vertex_it).as.pos[1]);
		}
		else
		{
			//do nothing
			//cout << "why vertex no position!" << endl;
		}
	}
	delete info;
	//all the pixels belong to vertex set and/or its neighbor should update its information
	set<MxVertexID> update_domain;//the domains need update
	update_domain.clear();
	MxVertexList star;
	for(set<MxVertexID>::iterator vertex_it = position_domain.begin(); vertex_it != position_domain.end(); vertex_it++)
	{
		star.reset();
		m->collect_vertex_star(*vertex_it, star);
		for(uint j=0; j<star.length(); j++)
		{
			update_domain.insert(star(j));
		}
	}
	for (set<MxVertexID>::iterator domain_it = update_domain.begin(); domain_it != update_domain.end(); domain_it++)
	{
		set<int>::iterator pixel_it;
		MxVertexID closest;
		float energy_saved_closest;
		MxVector from_update(dim());
		MxVector to_update(dim());

		//collect all the neighbor nodes
		MxVertexList star;
		star.reset();
		m->collect_vertex_star(*domain_it, star);

		for (pixel_it = vertex_domains->vertex_domain[*domain_it].begin(); pixel_it != vertex_domains->vertex_domain[*domain_it].end(); pixel_it++)
		{
			closest = *domain_it;
			energy_saved_closest=0;
			pixel_info *info = new pixel_info(dim());

			//compute the closest vertex
			for(uint j=0; j<star.length(); j++)
			{
				//float dist = vertex_dist(*pixel_it, star(j));
				float energy_saved = calc_energy_delta(*domain_it,star(j),*pixel_it, info);
				if ( energy_saved> energy_saved_closest)
				{
					from_update = info->from_temp_update;
					to_update = info->to_temp_update;
					closest = star(j);
					energy_saved_closest = energy_saved;
				}
			}
			if(energy_saved_closest >1e-5)
			{			
				if (pixel_id_pointer[*pixel_it] == NULL)
					//insert to heap
					heap_pixel_swap(energy_saved_closest,*domain_it,closest,*pixel_it, from_update,to_update, info);
				else
				{
					//modify the existing one
					pixel_info *exist_info = pixel_id_pointer[*pixel_it];
					assert(exist_info->is_in_heap());
					//from, id pointer should be the same
					assert(*domain_it == exist_info->from_id);
					pixel_to_set[exist_info->to_id].erase(exist_info);
					exist_info->to_id = closest;
					pixel_to_set[exist_info->to_id].insert(exist_info);
					exist_info->from_update = from_update;
					exist_info->to_update = to_update;
					exist_info->heap_key(energy_saved_closest);
					pixel_heap.update(exist_info);
					delete info;
				}
			}	
			else
			{
				//delete the existing one
				if (pixel_id_pointer[*pixel_it] != NULL)
				{
					pixel_info *exist_info = pixel_id_pointer[*pixel_it];				
					pixel_id_pointer[*pixel_it] = NULL;
					assert(pixel_from_set[exist_info->from_id].erase(exist_info));
					assert(pixel_to_set[exist_info->to_id].erase(exist_info));
					pixel_heap.remove(exist_info);
				}
				delete info;
			}
		}
	}
}

void MxPropSlim::update_vertex_position(pixel_info * info)
{
	quadric(info->from_id) -= init_quadric(info->pixel_id);
	quadric(info->to_id) += init_quadric(info->pixel_id);

	 m->vertex(info->from_id).as.pos[0] = info->from_update[0];
	 m->vertex(info->from_id).as.pos[1] = info->from_update[1];
	 m->vertex(info->from_id).as.pos[2] = info->from_update[2];

	 m->vertex(info->to_id).as.pos[0] = info->to_update[0];
	 m->vertex(info->to_id).as.pos[1] = info->to_update[1];
	 m->vertex(info->to_id).as.pos[2] = info->to_update[2];

	 trajectory[info->from_id].push_back(m->vertex(info->from_id).as.pos[0]);
	 trajectory[info->from_id].push_back(m->vertex(info->from_id).as.pos[1]);
	 trajectory[info->to_id].push_back(m->vertex(info->to_id).as.pos[0]);
	 trajectory[info->to_id].push_back(m->vertex(info->to_id).as.pos[1]);
}


void MxPropSlim::lloyd_swap_set(set<pixel_info *> & info_set)
{
	set<MxVertexID> position_domain;// the vertex set that need position update
	position_domain.clear();
	for (set<pixel_info *>::iterator info_it = info_set.begin(); info_it != info_set.end(); info_it++)
	{
		pixel_info * info = *info_it;
		assert(vertex_domains->vertex_domain[info->from_id].erase(info->pixel_id));
		vertex_domains->vertex_domain[info->to_id].insert(info->pixel_id);

		//delete link information
		assert(pixel_from_set[info->from_id].erase(info));
		assert(pixel_to_set[info->to_id].erase(info));
		pixel_id_pointer[info->pixel_id] = NULL;

		////aux
		//vertex_domains->red[info->pixel_id] = vertex_domains->red[info->to_id];
		//vertex_domains->green[info->pixel_id] = vertex_domains->green[info->to_id];
		//vertex_domains->blue[info->pixel_id] = vertex_domains->blue[info->to_id];

		//update quadric of v1 and v2 
		quadric(info->from_id) -= init_quadric(info->pixel_id);
		quadric(info->to_id) += init_quadric(info->pixel_id);

		//collect the vertexs that need position update
		position_domain.insert(info->from_id);
		position_domain.insert(info->to_id);
	}
	//update vertex position
	pixel_info * info = new pixel_info(dim());// for temp use, no meaning
	for(set<MxVertexID>::iterator vertex_it = position_domain.begin(); vertex_it != position_domain.end(); vertex_it++)
	{
		const MxQuadric &Q =quadric(*vertex_it);
		if(Q.optimize(info->to_update))
		{
			m->vertex(*vertex_it).as.pos[0] = info->to_update[0];
			m->vertex(*vertex_it).as.pos[1] = info->to_update[1];
			m->vertex(*vertex_it).as.pos[2] = info->to_update[2];

			trajectory[*vertex_it].push_back(m->vertex(*vertex_it).as.pos[0]);
			trajectory[*vertex_it].push_back(m->vertex(*vertex_it).as.pos[1]);
		}
		else
		{
			//do nothing
		}
	}
	delete info;
	//all the pixels belong to vertex set and/or its neighbor should update its information
	set<MxVertexID> update_domain;//the domains need update
	update_domain.clear();
	MxVertexList star;
	for(set<MxVertexID>::iterator vertex_it = position_domain.begin(); vertex_it != position_domain.end(); vertex_it++)
	{
		update_domain.insert(*vertex_it);
		star.reset();
		m->collect_vertex_star(*vertex_it, star);
		for(uint j=0; j<star.length(); j++)
		{
			update_domain.insert(star(j));
		}
	}
	for (set<MxVertexID>::iterator domain_it = update_domain.begin(); domain_it != update_domain.end(); domain_it++)
	{
		set<int>::iterator pixel_it;
		MxVertexID closest;
		float distance_closest;

		//collect all the neighbor nodes
		MxVertexList star;
		star.reset();
		m->collect_vertex_star(*domain_it, star);

		for (pixel_it = vertex_domains->vertex_domain[*domain_it].begin(); pixel_it != vertex_domains->vertex_domain[*domain_it].end(); pixel_it++)
		{
			closest = *domain_it;
			//initially the closest one is the domain itself
			distance_closest= vertex_dist(*pixel_it, *domain_it);
			pixel_info *info = new pixel_info(dim());

			//compute the closest vertex
			for(uint j=0; j<star.length(); j++)
			{
				float dist = vertex_dist(*pixel_it, star(j));
				if ( dist< distance_closest)
				{
					closest = star(j);
					distance_closest = dist;
				}
			}
			if(closest != *domain_it)
			{			
				if (pixel_id_pointer[*pixel_it] == NULL)
					//insert to heap
					heap_lloyd_swap(vertex_dist(*pixel_it, *domain_it) - distance_closest,*domain_it,closest,*pixel_it, info);
				else
				{
					//modify the existing one
					pixel_info *exist_info = pixel_id_pointer[*pixel_it];
					assert(exist_info->is_in_heap());
					//from, id pointer should be the same
					assert(*domain_it == exist_info->from_id);
					pixel_to_set[exist_info->to_id].erase(exist_info);
					exist_info->to_id = closest;
					pixel_to_set[exist_info->to_id].insert(exist_info);
					exist_info->heap_key(vertex_dist(*pixel_it, *domain_it)-distance_closest);
					pixel_heap.update(exist_info);
					delete info;
				}
			}	
			else
			{
				//delete the existing one
				if (pixel_id_pointer[*pixel_it] != NULL)
				{
					pixel_info *exist_info = pixel_id_pointer[*pixel_it];				
					pixel_id_pointer[*pixel_it] = NULL;
					assert(pixel_from_set[exist_info->from_id].erase(exist_info));
					assert(pixel_to_set[exist_info->to_id].erase(exist_info));
					pixel_heap.remove(exist_info);
				}
				delete info;
			}
		}
	}
}

void MxPropSlim::heap_lloyd_swap(float distance_saved,MxVertexID from, MxVertexID to, MxVertexID init_vertex_id, pixel_info *info)
{

	info->from_id = from;
	info->to_id = to;
	info->pixel_id = init_vertex_id;
	info->heap_key(distance_saved);

	pixel_from_set[from].insert(info);
	pixel_to_set[to].insert(info);
	pixel_id_pointer[init_vertex_id] = info;

	pixel_heap.insert(info);
}

float MxPropSlim::calc_energy_delta(MxVertexID from, MxVertexID to, MxVertexID init_vertex_id, pixel_info *info)
{
	float energy_from;
	float energy_to;
	//only care about the energy before swapping, no need to use the original position
	const MxQuadric &Qi=quadric(from), &Qj=quadric(to), &Q_init = init_quadric(init_vertex_id);
	MxQuadric Q_from = Qi;
	MxQuadric Q_to = Qj;
	if(Q_from.optimize(info->from_temp_update))
		energy_from = Q_from.getmin(info->from_temp_update);
	else //if rank(A) = 0, exist an x satisfy (x-p)tA(x-p)=0;
		energy_from = 0;

	if(Q_to.optimize(info->to_temp_update))
		energy_to = Q_to.getmin(info->to_temp_update);
	else
		energy_to = 0;
	float energy_before = energy_from + energy_to;

	Q_from -= Q_init;
	Q_to += Q_init;

	if(Q_from.optimize(info->from_temp_update))
		energy_from = Q_from.getmin(info->from_temp_update);
	else //if rank(A) = 0, exist an x satisfy (x-p)tA(x-p)=0;
		energy_from = 0;

	if(Q_to.optimize(info->to_temp_update))
		energy_to = Q_to.getmin(info->to_temp_update);
	else
		energy_to = 0;
	info->from_temp_id = from;
	info->to_temp_id = to;
	info->pixel_id = init_vertex_id;
	float energy_after = energy_from + energy_to;
	return energy_before - energy_after;
}

////////////////////////////////////////////////////////////////////////
//
// These were copied *unmodified* from MxEdgeQSlim
// (with some unsupported features commented out).
//

void MxPropSlim::collect_edges()
{
    MxVertexList star;

    for(MxVertexID i=0; i<m->vert_count(); i++)
    {
        star.reset();
        m->collect_vertex_star(i, star);

        for(uint j=0; j<star.length(); j++)
            if( i < star(j) )  // Only add particular edge once
			{
                create_edge(i, star(j));
				/*neighbor_node[i].insert(star(j));
				neighbor_node[star(j)].insert(i);*/
			}
    }
}

void MxPropSlim::collect_domain_pixel(MxVertexID i)
{
	set<int>::iterator pixel_it;
	MxVertexID closest;
	float energy_saved_closest;
	MxVector from_update(dim());
	MxVector to_update(dim());

	//collect all the neighbor nodes
	MxVertexList star;
	star.reset();
	m->collect_vertex_star(i, star);

	for (pixel_it = vertex_domains->vertex_domain[i].begin(); pixel_it != vertex_domains->vertex_domain[i].end(); pixel_it++)
	{
		closest = i;
		energy_saved_closest=0;
		pixel_info *info = new pixel_info(dim());

		//compute the closest vertex
		for(uint j=0; j<star.length(); j++)
		{
			//float dist = vertex_dist(*pixel_it, star(j));
			float energy_saved = calc_energy_delta(i,star(j),*pixel_it, info);
			if ( energy_saved> energy_saved_closest)
			{
				from_update = info->from_temp_update;
				to_update = info->to_temp_update;
				closest = star(j);
				energy_saved_closest = energy_saved;
			}
		}
		// insert
		if(energy_saved_closest > 0)
			heap_pixel_swap(energy_saved_closest,i,closest,*pixel_it, from_update, to_update, info);
		//no need to insert
		else
			delete info;
	}
}

void MxPropSlim::init_pswap()
{
	trajectory = new vector<float>[m->vert_count()];
	pixel_from_set = new pixel_set[m->vert_count()];
	pixel_to_set = new pixel_set[m->vert_count()];
	pixel_id_pointer = new pixel_info*[m->vert_count()];
	for (int i=0; i<m->vert_count(); i++)
	{
		trajectory[i].clear();
		trajectory[i].push_back(m->vertex(i).as.pos[0]);
		trajectory[i].push_back(m->vertex(i).as.pos[1]);
		pixel_from_set[i].clear();
		pixel_to_set[i].clear();
		pixel_id_pointer[i]= NULL;
	}

	for(MxVertexID i=0; i<m->vert_count(); i++)
	{
		if (!vertex_domains->vertex_domain[i].empty())
		{
			collect_domain_pixel(i);
		} 
	}

	////here is a potential special case to distinguish this method from lloyd
	//spec_ex_lloyd();
}

void MxPropSlim::spec_ex_lloyd()
{
	decimate(4);
	int count = 0;
	int valid_v[4];
	for(MxVertexID i=0; i<m->vert_count(); i++)
	{		
		if (!vertex_domains->vertex_domain[i].empty())
		{
			if (count == 0)
				valid_v[count] = i;
			else if (count == 1)
				valid_v[count]= i;
			else if (count == 2)
				valid_v[count] = i;

			else if (count == 3)
				valid_v[count] = i;
			count ++;
		} 
	}
	//no disturb
	for (int i=0; i<4; i++)
	{
		vertex_domains->vertex_domain[valid_v[i]].clear();
		quadric(valid_v[i]) = init_quadric(valid_v[i]);
	}
	for(MxVertexID i=0; i<m->vert_count(); i++)
	{
		if (init_model->vertex(i).as.pos[0] <= 127 && init_model->vertex(i).as.pos[1] <= 65)
		{
			vertex_domains->vertex_domain[valid_v[0]].insert(i);
			quadric(valid_v[0]) += init_quadric(i);
		} 
		else if (init_model->vertex(i).as.pos[0] > 127 && init_model->vertex(i).as.pos[1] <= 65)
		{
			vertex_domains->vertex_domain[valid_v[1]].insert(i);
			quadric(valid_v[1]) += init_quadric(i);
		} 
		else if (init_model->vertex(i).as.pos[0] <= 127 && init_model->vertex(i).as.pos[1] > 65)
		{
			vertex_domains->vertex_domain[valid_v[2]].insert(i);
			quadric(valid_v[2]) += init_quadric(i);
		} 
		else if (init_model->vertex(i).as.pos[0] > 127 && init_model->vertex(i).as.pos[1] > 65)
		{
			vertex_domains->vertex_domain[valid_v[3]].insert(i);
			quadric(valid_v[3]) += init_quadric(i);
		} 
	}

	////try disturb somehow
	//float angle = 10/180.0 * 3.14159;
	//for (int i=0; i<4; i++)
	//{
	//	vertex_domains->vertex_domain[valid_v[i]].clear();
	//	quadric(valid_v[i]) = init_quadric(valid_v[i]);
	//}
	//for(MxVertexID i=0; i<m->vert_count(); i++)
	//{
	//	if ((init_model->vertex(i).as.pos[1] - 65) < tan(angle) *(init_model->vertex(i).as.pos[0] - 127) && (init_model->vertex(i).as.pos[1] - 65) < -1/tan(angle) *(init_model->vertex(i).as.pos[0] - 127) )
	//	{
	//		vertex_domains->vertex_domain[valid_v[0]].insert(i);
	//		quadric(valid_v[0]) += init_quadric(i);
	//	} 
	//	else if ((init_model->vertex(i).as.pos[1] - 65) >= tan(angle) *(init_model->vertex(i).as.pos[0] - 127) && (init_model->vertex(i).as.pos[1] - 65) < -1/tan(angle) *(init_model->vertex(i).as.pos[0] - 127) )
	//	{
	//		vertex_domains->vertex_domain[valid_v[1]].insert(i);
	//		quadric(valid_v[1]) += init_quadric(i);
	//	} 
	//	else if ((init_model->vertex(i).as.pos[1] - 65) < tan(angle) *(init_model->vertex(i).as.pos[0] - 127) && (init_model->vertex(i).as.pos[1] - 65) >= -1/tan(angle) *(init_model->vertex(i).as.pos[0] - 127) )
	//	{
	//		vertex_domains->vertex_domain[valid_v[2]].insert(i);
	//		quadric(valid_v[2]) += init_quadric(i);
	//	} 
	//	else if ((init_model->vertex(i).as.pos[1] - 65) >= tan(angle) *(init_model->vertex(i).as.pos[0] - 127) && (init_model->vertex(i).as.pos[1] - 65) >= -1/tan(angle) *(init_model->vertex(i).as.pos[0] - 127) )
	//	{
	//		vertex_domains->vertex_domain[valid_v[3]].insert(i);
	//		quadric(valid_v[3]) += init_quadric(i);
	//	} 
	//}

	for (int i=0; i<4; i++)
	{
		const MxQuadric &Q=quadric(valid_v[i]);
		MxVector v(dim());
		assert(Q.optimize(v));
		m->vertex(valid_v[i]).as.pos[0] = v[0];
		m->vertex(valid_v[i]).as.pos[1] = v[1];
		
	}
	trajectory = new vector<float>[m->vert_count()];
	pixel_from_set = new pixel_set[m->vert_count()];
	pixel_to_set = new pixel_set[m->vert_count()];
	pixel_id_pointer = new pixel_info*[m->vert_count()];
	for (int i=0; i<m->vert_count(); i++)
	{
		trajectory[i].clear();
		trajectory[i].push_back(m->vertex(i).as.pos[0]);
		trajectory[i].push_back(m->vertex(i).as.pos[1]);
		pixel_from_set[i].clear();
		pixel_to_set[i].clear();
		pixel_id_pointer[i]= NULL;
	}

	for(MxVertexID i=0; i<m->vert_count(); i++)
	{
		if (!vertex_domains->vertex_domain[i].empty())
		{
			collect_domain_pixel(i);
		} 
	}
}

void MxPropSlim::init_lloyd()
{
	trajectory = new vector<float>[m->vert_count()];
	pixel_from_set = new pixel_set[m->vert_count()];
	pixel_to_set = new pixel_set[m->vert_count()];
	pixel_id_pointer = new pixel_info*[m->vert_count()];
	for (int i=0; i<m->vert_count(); i++)
	{
		trajectory[i].clear();
		trajectory[i].push_back(m->vertex(i).as.pos[0]);
		trajectory[i].push_back(m->vertex(i).as.pos[1]);
		pixel_from_set[i].clear();
		pixel_to_set[i].clear();
		pixel_id_pointer[i]= NULL;
	}

	for(MxVertexID i=0; i<m->vert_count(); i++)
	{
		if (!vertex_domains->vertex_domain[i].empty())
		{
			collect_lloyd_pixel(i);
		} 
	}
}

void MxPropSlim::collect_lloyd_pixel(MxVertexID i)
{
	set<int>::iterator pixel_it;
	MxVertexID closest;
	float distance_closest;

	//collect all the neighbor nodes
	MxVertexList star;
	star.reset();
	m->collect_vertex_star(i, star);
	
	for (pixel_it = vertex_domains->vertex_domain[i].begin(); pixel_it != vertex_domains->vertex_domain[i].end(); pixel_it++)
	{
		closest = i;
		distance_closest=vertex_dist(*pixel_it, i);
		pixel_info *info = new pixel_info(dim());

		//compute the closest vertex
		for(uint j=0; j<star.length(); j++)
		{
			float dist = vertex_dist(*pixel_it, star(j));
			if ( dist < distance_closest)
			{
				closest = star(j);
				distance_closest = dist;
			}
		}
		// insert
		if(distance_closest < vertex_dist(*pixel_it, i))
			heap_lloyd_swap(vertex_dist(*pixel_it, i) - distance_closest,i,closest,*pixel_it, info);
		//no need to insert
		else
			delete info;
	}
}

void MxPropSlim::heap_pixel_swap(float energy_saved,MxVertexID from, MxVertexID to, MxVertexID init_vertex_id, MxVector from_update, MxVector to_update, pixel_info *info)
{
	
	info->from_id = from;
	info->to_id = to;
	info->pixel_id = init_vertex_id;
	info->from_update = from_update;
	info->to_update = to_update;
	info->heap_key(energy_saved);

	pixel_from_set[from].insert(info);
	pixel_to_set[to].insert(info);
	pixel_id_pointer[init_vertex_id] = info;
	
	pixel_heap.insert(info);
}

void MxPropSlim::print_trajectory_length()
{
	float last_x;
	float last_y;
	float x;
	float y;
	float total_length = 0;
	for (int i=0; i<m->vert_count(); i++)
	{
		if (!vertex_domains->vertex_domain[i].empty())
		{
			vector<float>::iterator iter = trajectory[i].begin();
			int count = 0;
			while (iter != trajectory[i].end())
			{
				if (count == 0)
				{
					last_x = x = *iter;
					iter++;
					last_y = y = *iter;
					iter++;
				}
				else
				{
					last_x = x;
					last_y = y;
					x = *iter;
					iter++;
					y = *iter;
					iter++;
					total_length += sqrt((y-last_y)*(y-last_y)+(x-last_x)*(x-last_x));
				}
				count++;
			}
		}				
	}
	cout << "Final trajectory length is " << total_length << endl;
}

void MxPropSlim::print_total_energy()
{
	float energy = 0;
	MxVector temp_v(dim());
	for(MxVertexID i=0; i<m->vert_count(); i++)
	{
		if (!vertex_domains->vertex_domain[i].empty())
		{
			//the energy is defined via the center
			pack_to_vector(i,temp_v);
			for (set<int>::iterator it = vertex_domains->vertex_domain[i].begin(); it != vertex_domains->vertex_domain[i].end(); it++)
			{
				const MxQuadric &Qi=init_quadric(*it);
				energy += Qi(temp_v);
			}
		} 
	}
	cout << "energy is " << energy << endl;
}

void MxPropSlim::constrain_boundaries()
{
    MxVertexList star;
    MxFaceList faces;

    for(MxVertexID i=0; i<m->vert_count(); i++)
    {
	star.reset();
	m->collect_vertex_star(i, star);

	for(uint j=0; j<star.length(); j++)
	    if( i < star(j) )
	    {
		faces.reset();
		m->collect_edge_neighbors(i, star(j), faces);
		if( faces.length() == 1 )
		    discontinuity_constraint(i, star(j), faces);
	    }
    }
}

void MxPropSlim::compute_edge_info(edge_info *info)
{
    compute_target_placement(info);

//     if( will_normalize_error )
//     {
//         double e_max = Q_max(info->vnew);
//         if( weight_by_area )
//             e_max /= Q_max.area();

//         info->heap_key(info->heap_key() / e_max);
//     }

    finalize_edge_update(info);
}

double MxPropSlim::check_local_compactness(uint v1, uint/*v2*/,
	const float *vnew)
{
	const MxFaceList& N1 = m->neighbors(v1);
	double c_min = 1.0;

	for(uint i=0; i<N1.length(); i++)
		if( m->face_mark(N1[i]) == 1 )
		{
			const MxFace& f = m->face(N1[i]);
			Vec3 f_after[3];
			for(uint j=0; j<3; j++)
				f_after[j] = (f[j]==v1)?Vec3(vnew):Vec3(m->vertex(f[j]));

			double c=triangle_compactness(f_after[0], f_after[1], f_after[2]);

			if( c < c_min ) c_min = c;
		}

		return c_min;
}

uint MxPropSlim::check_local_validity(uint v1, uint /*v2*/, const float *vnew)

{
	const MxFaceList& N1 = m->neighbors(v1);
	uint nfailed = 0;
	uint i;

	for(i=0; i<N1.length(); i++)
		if( m->face_mark(N1[i]) == 1 )
		{
			MxFace& f = m->face(N1[i]);
			uint k = f.find_vertex(v1);
			uint x = f[(k+1)%3];
			uint y = f[(k+2)%3];

			float d_yx[3], d_vx[3], d_vnew[3], f_n[3], n[3];
			mxv_sub(d_yx, m->vertex(y),  m->vertex(x), 3);   // d_yx = y-x
			mxv_sub(d_vx, m->vertex(v1), m->vertex(x), 3);   // d_vx = v-x
			mxv_sub(d_vnew, vnew, m->vertex(x), 3);          // d_vnew = vnew-x

			mxv_cross3(f_n, d_yx, d_vx);
			mxv_cross3(n, f_n, d_yx);     // n = ((y-x)^(v-x))^(y-x)
			mxv_unitize(n, 3);

			// assert( mxv_dot(d_vx, n, 3) > -FEQ_EPS );
			if(mxv_dot(d_vnew,n,3)<local_validity_threshold*mxv_dot(d_vx,n,3))
				nfailed++;
		}

		return nfailed;
}

void MxPropSlim::apply_mesh_penalties(edge_info *info)
{
	uint i;

	const MxFaceList& N1 = m->neighbors(info->v1);
	const MxFaceList& N2 = m->neighbors(info->v2);

	// Set up the face marks as the check_xxx() functions expect.
	//
	for(i=0; i<N2.length(); i++) m->face_mark(N2[i], 0);
	for(i=0; i<N1.length(); i++) m->face_mark(N1[i], 1);
	for(i=0; i<N2.length(); i++) m->face_mark(N2[i], m->face_mark(N2[i])+1);

	double base_error = info->heap_key();
	double bias = 0.0;

	// Check for excess over degree bounds.
	//
	uint max_degree = MAX(N1.length(), N2.length());
	if( max_degree > vertex_degree_limit )
		bias += (max_degree-vertex_degree_limit) * meshing_penalty * 0.001;

#if ALTERNATE_DEGREE_BIAS
	// ??BUG:  This code was supposed to be a slight improvement over
	//         the earlier version above.  But it performs worse.
	//         Should check into why sometime.
	//
	uint degree_excess = check_local_degree(info->v1, info->v2, info->vnew);
	if( degree_excess )
		bias += degree_excess * meshing_penalty;
#endif

	// Local validity checks
	//
	float vnew[3] = {info->target[0],info->target[1],info->target[2]};
	uint nfailed = check_local_validity(info->v1, info->v2, vnew);
	nfailed += check_local_validity(info->v2, info->v1, vnew);
	if( nfailed )
		bias += nfailed*meshing_penalty;

	if( compactness_ratio > 0.0 )
	{
		double c1_min=check_local_compactness(info->v1, info->v2, vnew);
		double c2_min=check_local_compactness(info->v2, info->v1, vnew);
		double c_min = MIN(c1_min, c2_min);

		// !!BUG: There's a small problem with this: it ignores the scale
		//        of the errors when adding the bias.  For instance, enabling
		//        this on the cow produces bad results.  I also tried
		//        += (base_error + FEQ_EPS) * (2-c_min), but that works
		//        poorly on the flat planar thing.  A better solution is
		//        clearly needed.
		//
		//  NOTE: The prior heuristic was
		//        if( ratio*cmin_before > cmin_after ) apply penalty;
		//
		if( c_min < compactness_ratio )
			bias += (1-c_min);
	}

#if USE_OLD_INVERSION_CHECK
	double Nmin1 = check_local_inversion(info->v1, info->v2, info->vnew);
	double Nmin2 = check_local_inversion(info->v2, info->v1, info->vnew);
	if( MIN(Nmin1, Nmin2) < 0.0 )
		bias += meshing_penalty;
#endif

	info->heap_key(base_error - bias);
}


void MxPropSlim::finalize_edge_update(edge_info *info)
{
     if( meshing_penalty > 1.0 )
         apply_mesh_penalties(info);

    if( info->is_in_heap() )
        heap.update(info);
    else
        heap.insert(info);
}

void MxPropSlim::update_pre_contract(const MxPairContraction& conx)
{
    MxVertexID v1=conx.v1, v2=conx.v2;
    uint i, j;

    star.reset();
    m->collect_vertex_star(v1, star);

    for(i=0; i<edge_links(v2).length(); i++)
    {
        edge_info *e = edge_links(v2)(i);
        MxVertexID u = (e->v1==v2)?e->v2:e->v1;
        SanityCheck( e->v1==v2 || e->v2==v2 );
        SanityCheck( u!=v2 );

        if( u==v1 || varray_find(star, u) )
        {
            // This is a useless link --- kill it
            bool found = varray_find(edge_links(u), e, &j);
            assert( found );
            edge_links(u).remove(j);
            heap.remove(e);
            if( u!=v1 ) delete e; // (v1,v2) will be deleted later
        }
        else
        {
            // Relink this to v1
            e->v1 = v1;
            e->v2 = u;
            edge_links(v1).add(e);
        }
    }

    edge_links(v2).reset();
}
