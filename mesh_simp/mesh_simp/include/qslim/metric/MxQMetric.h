#ifndef MXQMETRIC_INCLUDED // -*- C++ -*-
#define MXQMETRIC_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  n-D Quadric Error Metric

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxQMetric.h,v 1.5 1998/10/26 21:09:17 garland Exp $

 ************************************************************************/

#include "MxQMetric3.h"
#include "MxVector.h"
#include "MxMatrix.h"

class MxQuadric
{
private:
    MxMatrix A;
    MxVector b;
    double c;

    double r;

public:
	MxQuadric(const float x_coor, const float y_coor, const float z_coor, int N=3);
	MxQuadric(double x_coor, double y_coor, double z_coor, double nx, double ny, double nz, int N=3);
	MxQuadric(float * x_coor, float * y_coor, float *z_coor, float* n_to_return, int np, int N=3);
    MxQuadric(unsigned int N) : A(N), b(N) { clear(); }
    MxQuadric(const MxVector& p1, const MxVector& p2, const MxVector& p3,
	      double area=1.0);
	MxQuadric(const MxVector& p1, const double coefficient);	// Compute an isotropic constraint quadrics
	MxQuadric(const MxVector& p1, const MxVector& p2);	// Used for adding boundary constraint (Added by Chao Wang)
    //MxQuadric(const MxQuadric3&, unsigned int N);
    MxQuadric(const MxQuadric& Q)
	: A(Q.A.dim()), b(Q.b.dim()) { *this = Q; }

    const MxMatrix& tensor() const { return A; }
    const MxVector& vector() const { return b; }
    double offset() const { return c; }
    double area() const { return r; }
    MxMatrix& homogeneous(MxMatrix& H) const;

    void clear(double val=0.0) { A=val; b=val; c=val; r=val; }
    MxQuadric& operator=(const MxQuadric& Q)
	{ A=Q.A; b=Q.b; c=Q.c; r=Q.r; return *this; }
    MxQuadric& operator+=(const MxQuadric& Q)
	{ A+=Q.A; b+=Q.b; c+=Q.c; r+=Q.r; return *this; }
    MxQuadric& operator-=(const MxQuadric& Q)
	{ A-=Q.A; b-=Q.b; c-=Q.c; r-=Q.r; return *this; }
    MxQuadric& operator*=(double s)
	{ A*=s; b*=s; c*=s; return *this; }

	// Compute the quadrics matrix with a vertex as input
    double evaluate(const MxVector& v) const;
	// Overloaded operator: Compute the quadrics matrix with a vertex as input
    double operator()(const MxVector& v) const { return evaluate(v); }
	// Compute the contract vertex v
	// If A^(-1), compute contract vertex v = A^(-1)*B and return true;
	// If A^(-1) does not exist, return false.
    bool optimize(MxVector& v) const;
	double getmin(MxVector& v) const;
	void isospace(MxVector& v);
	void decompose(MxVector& axis1, MxVector& axis2, MxVector& axis3, double& l1, double& l2, double& l3);
};

// MXQMETRIC_INCLUDED
#endif
