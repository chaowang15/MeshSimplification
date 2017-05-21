#ifndef MXPROPSLIM_INCLUDED // -*- C++ -*-
#define MXPROPSLIM_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  MxPropSlim

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxPropSlim.h,v 1.6 1998/10/26 21:09:15 garland Exp $

 ************************************************************************/

#include "MxStdSlim.h"
#include "MxQMetric.h"
#include <set>
class MxPropSlim : public MxStdSlim
{
private:
    uint D;

    bool use_color;
    bool use_texture;
    bool use_normals;

	class edge_info : public MxHeapable
	{
	public:
		MxVertexID v1, v2;
		MxVector target;

		edge_info(uint D) : target(D) { }
	};
    typedef MxSizedDynBlock<edge_info*, 6> edge_list;


    MxBlock<edge_list> edge_links;	// 1 per vertex
    MxBlock<MxQuadric*> __quadrics;	// 1 per vertex
	MxBlock<MxQuadric*> initial_quadrics;	// 1 per vertex

    //
    // Temporary variables used by methods
    MxVertexList star, star2;
    MxPairContraction conx_tmp;

	//pixel re-arrangement
	class pixel_info : public MxHeapable
	{
	public:
		MxVertexID from_id, to_id;
		MxVertexID pixel_id;
		//float energy_before;

		MxVector from_update;
		MxVector to_update;
		//used for temp
		MxVertexID from_temp_id;
		MxVertexID to_temp_id;
		MxVector from_temp_update;
		MxVector to_temp_update;

		pixel_info(uint D) : from_update(D),to_update(D),from_temp_update(D),to_temp_update(D) { }
	};
	typedef set<pixel_info*> pixel_set;
	pixel_set* pixel_from_set;
	pixel_set* pixel_to_set;
	pixel_info** pixel_id_pointer;
	MxHeap pixel_heap;

protected:
	void collect_domain_pixel(MxVertexID i);
	void collect_lloyd_pixel(MxVertexID i);
	float calc_energy_delta(MxVertexID from, MxVertexID to, MxVertexID init_vertex_id, pixel_info *info);
	void heap_pixel_swap(float energy_saved,MxVertexID from, MxVertexID to, MxVertexID init_vertex_id, MxVector from_update, MxVector to_update, pixel_info *info);
	void heap_lloyd_swap(float distance_saved,MxVertexID from, MxVertexID to, MxVertexID init_vertex_id, pixel_info *info);

	void apply_swap_per_pixel(pixel_info * info);
	void apply_swap_set(set<pixel_info *> & info_set);
	void lloyd_swap_set(set<pixel_info *> & info_set);
	void collect_update_domain(set<MxVertexID> & update_domain, pixel_info * info);

	void update_vertex_position(pixel_info * info);
	//void compute_pixel_info(pixel_info * info);

protected:
    uint compute_dimension(MxStdModel *);
    void pack_to_vector(MxVertexID, MxVector&);
    void unpack_from_vector(MxVertexID, MxVector&);
    uint prop_count();
    void pack_prop_to_vector(MxVertexID, MxVector&, uint);
    void unpack_prop_from_vector(MxVertexID, MxVector&, uint);

    void compute_face_quadric(MxFaceID, MxQuadric&);
    void collect_quadrics();

    void create_edge(MxVertexID, MxVertexID);
    void collect_edges();
    void constrain_boundaries();
    void discontinuity_constraint(MxVertexID, MxVertexID, const MxFaceList&);
    void compute_edge_info(edge_info *);
    void finalize_edge_update(edge_info *);
    void compute_target_placement(edge_info *);

    void apply_contraction(const MxPairContraction&, edge_info *);
    void update_pre_contract(const MxPairContraction&);
	void apply_mesh_penalties(edge_info *info);
	uint check_local_validity(uint v1, uint /*v2*/, const float *vnew);
	double check_local_compactness(uint v1, uint/*v2*/,
		const float *vnew);

public:
    bool will_decouple_quadrics;

public:
    MxPropSlim(MxStdModel *);

    uint dim() const { return D; }

    void consider_color(bool will=true);
    void consider_texture(bool will=true);
    void consider_normals(bool will=true);

    uint quadric_count() const { return __quadrics.length(); }
    MxQuadric&       quadric(uint i)       { return *(__quadrics(i)); }
    const MxQuadric& quadric(uint i) const { return *(__quadrics(i)); }

	MxQuadric&       init_quadric(uint i)       { return *(initial_quadrics(i)); }
	const MxQuadric& init_quadric(uint i) const { return *(initial_quadrics(i)); }


    void initialize();
    bool decimate(uint);

	void decimateOnce();
	void getTopQslimEdge(int &v1, int &v2);
	bool getTopPswapPixel(int &from, int &to, int &id);
	float vertex_dist(MxVertexID init_model_vert, MxVertexID contracted_model_vert);
	//void updateMassCenter();
	void init_pswap();
	void init_lloyd();
	bool swap_per_pixel();
	bool swap_per_pixel_once();
	bool swap_pixel_set();
	bool swap_pixel_set_once();
	bool swap_pixel_lloyd();
	bool swap_pixel_lloyd_once();

	//print out criteria
	void print_trajectory_length();
	void print_total_energy();
	//test functions
	void spec_ex_lloyd();
};

// MXPROPSLIM_INCLUDED
#endif
