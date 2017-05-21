/************************************************************************

  n-D Quadric Error Metrics

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxQMetric.cxx,v 1.5 1998/10/26 21:09:17 garland Exp $

 ************************************************************************/

#include "stdmix.h"
#include "MxQMetric.h"
#include <set>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
using namespace Eigen;

static
void symmetric_subfrom(MxMatrix& A, const MxVector& a, const MxVector& b)
{
    for(int i=0; i != A.dim(); i++)  
		for(int j=0; j !=A.dim(); j++)
			A(i,j) -= a[i]*b[j];
}

double gassian(double d, double size)
{
	double para = -(d*d)/size*size;
	return exp(para);
}

MxQuadric::MxQuadric(double x_coor, double y_coor, double z_coor, double nx, double ny, double nz, int N):
A(N), b(N)
{
	Vector3d n;
	n(0) = nx;
	n(1) = ny;
	n(2) = nz;

	for(int i=0; i<N; i++)  
		for(int j=0; j<N; j++)
			A(i,j) = n(i)*n(j);

	for (int i=0; i< N; i++)
	{
		b(i) = -(A(i, 0) * x_coor + A(i, 1)* y_coor + A(i, 2)* z_coor);
	}
	c = - (b(0)*x_coor  + b(1)*y_coor + b(2)*z_coor);
	r = 1;
}

MxQuadric::MxQuadric(float * x_coor, float * y_coor, float *z_coor, float* n_to_return, int np, int N):
A(N), b(N)
{
	//Vector3d avg;
	//avg.setZero();
	//for (int i=0; i<np; i++)
	//{
	//	avg += Vector3d(x_coor[i], y_coor[i], z_coor[i]);
	//}
	//avg /= np;

	//Vector3d Y;
	//Matrix3d M;
	//M.setZero();
	//Vector3d n;
	//for (int i=0; i<np; i++)
	//{
	//	Y = Vector3d(x_coor[i], y_coor[i], z_coor[i]) - avg;
	//	M += Y * Y.transpose();
	//}
	//EigenSolver<Matrix3d> es(M);
	//int smallest = 0;
	//for (int i=1; i<3; i++)
	//{
	//	if (es.eigenvalues()[smallest].real() > es.eigenvalues()[i].real())
	//		smallest = i;
	//}
	//n = es.eigenvectors().col(smallest).real();
	//n = n.normalized();

	//n_to_return[0] = n(0);
	//n_to_return[1] = n(1);
	//n_to_return[2] = n(2);

	//for(int i=0; i<N; i++)  
	//	for(int j=0; j<N; j++)
	//		A(i,j) = n(i)*n(j);

	//for (int i=0; i< N; i++)
	//{
	//	b(i) = -(A(i, 0) * x_coor[0] + A(i, 1)* y_coor[0] + A(i, 2)* z_coor[0]);
	//}
	//c = - (b(0)*x_coor[0]  + b(1)*y_coor[0] + b(2)*z_coor[0]);
	//r = 1;

	Vector3d avg;
	avg.setZero();
	for (int i=0; i<np; i++)
	{
		avg += Vector3d(x_coor[i], y_coor[i], z_coor[i]);
	}
	avg /= np;

	double size;
	std::set<double> sizes;
	for (int i=0; i< np; i++)
	{
		Vector3d dist = Vector3d(x_coor[i], y_coor[i], z_coor[i]) - avg;
		size = dist.norm();
		sizes.insert(size);
	}
	std::set<double>::iterator it = sizes.begin();
	for(int i=0; i != sizes.size()/2; i++)
		it++;
	size = *(it);

	Vector3d Y;
	Matrix3d M;
	M.setZero();
	Vector3d n;
	for (int i=0; i<np; i++)
	{
		Y = Vector3d(x_coor[i], y_coor[i], z_coor[i]) - avg;
		M += gassian(Y.norm(), size) * Y * Y.transpose();
	}
	EigenSolver<Matrix3d> es(M);
	int smallest = 0;
	for (int i=1; i<3; i++)
	{
		if (es.eigenvalues()[smallest].real() > es.eigenvalues()[i].real())
			smallest = i;
	}
	n = es.eigenvectors().col(smallest).real();
	n = n.normalized();

	n_to_return[0] = (float)n(0);
	n_to_return[1] = (float)n(1);
	n_to_return[2] = (float)n(2);

	for(int i=0; i<N; i++)  
		for(int j=0; j<N; j++)
			A(i,j) = n(i)*n(j);

	for (int i=0; i< N; i++)
	{
		b(i) = -(A(i, 0) * x_coor[0] + A(i, 1)* y_coor[0] + A(i, 2)* z_coor[0]);
	}
	c = - (b(0)*x_coor[0]  + b(1)*y_coor[0] + b(2)*z_coor[0]);
	r = size*size;

	A *= r;
	b *= r;
	c *= r;
}

MxQuadric::MxQuadric(const float x_coor, const float y_coor, const float z_coor, int N):
A(N), b(N)
{
	//cvt case
	mxm_identity(A, A.dim());
	MxVector p(N);
	p[0] = x_coor;
	p[1] = y_coor;
	p[2] = z_coor;
	b -= p;
	c = p*p;
	r = 1;
}

MxQuadric::MxQuadric(const MxVector& p1,const MxVector& p2,const MxVector& p3,
		     double area)
    : A(p1.dim()), b(p1.dim())
{
	//original code
	AssertBound( p1.dim()==p2.dim() && p1.dim()==p3.dim() );

	MxVector e1=p2; e1-=p1; unitize(e1); // e1 = p2-p1; unitize
	MxVector e2=p3; e2-=p1;              // e2 = p3-p1
	MxVector t = e1;
	t*=e1*e2; e2-=t; unitize(e2);        // e2 = p3-p1-e1*(e1*(p3-p1)); unitize

	double p1e1 = p1*e1;
	double p1e2 = p1*e2;

	mxm_identity(A, A.dim());
	symmetric_subfrom(A, e1,e1);
	symmetric_subfrom(A, e2,e2);

	b=e1; b*=p1e1; t=e2; t*=p1e2; b += t; b -= p1;

	c = p1*p1 - p1e1*p1e1 - p1e2*p1e2;

	r = area;

	//add some identity part for isotropic constraint quadrics
	//double coefficient = 0.1;
	//for (int i=0; i< A.dim(); i++)
	//{
	//	A(i,i) += coefficient;
	//}
	//for (int i=0; i< A.dim(); i++)
	//{
	//	b(i) -= coefficient*p1(i);
	//}
	//c += p1*p1*coefficient;

	////cvt case
	//mxm_identity(A, A.dim());
	//b -= p1;
	//c = p1*p1;
	//r = 1;
}

// Compute an isotropic constraint quadrics
MxQuadric::MxQuadric(const MxVector& p1, const double coefficient): A(p1.dim()), b(p1.dim())
{
	for (int i=0; i< A.dim(); i++)
	{
		A(i,i) = coefficient;
	}
	for (int i=0; i< A.dim(); i++)
	{
		b(i) = -coefficient*p1(i);
	}
	c = p1*p1*coefficient;
}

MxQuadric::MxQuadric(const MxVector& p1, const MxVector& p2) : A(p1.dim()), b(p1.dim())
{
	AssertBound( p1.dim()==p2.dim());
	MxVector e1=p2; e1-=p1; unitize(e1); // e1 = p2-p1; unitize
	mxm_identity(A, A.dim());
	symmetric_subfrom(A, e1,e1);
	b = A*p1;	// use the first endpoint as p
	mxv_neg(b, b.dim());
	c = p1*(A*p1);
	r = 1.0;
}

//MxQuadric::MxQuadric(const MxQuadric3& Q3, unsigned int N)
//    : A(N), b(N)
//{
//    uint i, j;
//
//    clear();
//
//    Mat3 A3 = Q3.tensor();
//    Vec3 b3 = Q3.vector();
//
//    for(i=0; i<3; i++)
//    {
//	for(j=0; j<3; j++)
//	    A(i,j) = A3(i,j);
//
//	b[i] = b3[i];
//    }
//
//    c = Q3.offset();
//    r = Q3.area();
//}

MxMatrix& MxQuadric::homogeneous(MxMatrix& H) const
{
    AssertBound( H.dim() == A.dim()+1 );

    uint i, j;

    for(i=0; i<A.dim(); i++)  for(j=0; j<A.dim(); i++)
	H(i,j) = A(i,j);

    for(i=0; i<b.dim(); i++)
	H(i, b.dim()) = H(b.dim(), i) = b[i];

    H(b.dim(), b.dim()) = c;

    return H;
}

double MxQuadric::evaluate(const MxVector& v) const
{
    AssertBound( v.dim() == b.dim() );
    return v*(A*v) + 2*(b*v) + c;
}

double MxQuadric::getmin(MxVector& v) const
{
	AssertBound( v.dim() == b.dim() );
	return (b*v) + c;
}

bool MxQuadric::optimize(MxVector& v) const
{
	//original code
	MxMatrix Ainv(A.dim());
	// Compute inverse A
	double det = A.invert(Ainv);
	if( FEQ(det, 0.0, 1e-12) )
		return false;

	v = (Ainv * b);
	mxv_neg(v, v.dim());
	return true;

	////general case
	//MatrixXd left(A.dim(),A.dim());
	//VectorXd right(A.dim());
	//for (int i=0; i<A.dim(); i++)
	//{
	//	right(i)=b[i];
	//	for (int j=0; j<A.dim(); j++)
	//	{
	//		left(i,j) = A[j*A.dim()+i];
	//	}			
	//}

	//LLT<MatrixXd> llt(left);		
	//VectorXd result = llt.solve(-right);

	//for (int i=0; i<A.dim(); i++)
	//{
	//	v[i] = result(i);
	//}
	//return true;

	////////special case to test speed optimized by eigen
	//Matrix3d left;
	//Vector3d right;
	//for (int i=0; i<3; i++)
	//{
	//	right(i)=b[i];
	//	for (int j=0; j<3; j++)
	//	{
	//		left(i,j) = A[j*3+i];
	//	}			
	//}

	//LLT<Matrix3d> llt(left);		
	//Vector3d result = llt.solve(-right);

	//for (int i=0; i<3; i++)
	//{
	//	v[i] = result(i);
	//}
	//return true;

	//////special case to test speed optimized by eigen
	//Matrix3d left;
	//Vector3d right;
	//for (int i=0; i<3; i++)
	//{
	//	right(i)=b[i];
	//	for (int j=0; j<3; j++)
	//	{
	//		left(i,j) = A[j*3+i];
	//	}			
	//}
	//	
	//Vector3d result = left.ldlt().solve(right);
	//result = -result;
	//for (int i=0; i<3; i++)
	//{
	//	v[i] = result(i);
	//}
	//return true;
}

void MxQuadric::isospace(MxVector& v)
{
	Matrix3d left;
	Vector3d m_space;
	for (int i=0; i<3; i++)
	{
		m_space(i) = v[i];
		for (int j=0; j<3; j++)
		{
			left(i,j) = A[j*3+i];
		}			
	}

	LLT<Matrix3d> llt(left);		
	Vector3d iso_space = m_space * llt.matrixL();

	for (int i=0; i<3; i++)
	{
		v[i] = iso_space(i);
	}
}

void MxQuadric::decompose(MxVector& axis1, MxVector& axis2, MxVector& axis3, double& l1, double& l2, double& l3)
{
	Matrix3d M;
	M.setZero();

	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			M(i,j) = A(i,j);
		}
	}
	EigenSolver<Matrix3d> es(M);
	Vector3d n = es.eigenvectors().col(0).real();
	l1 = 1.0/sqrt(es.eigenvalues()[0].real());
	axis1(0) = n(0);
	axis1(1) = n(1);
	axis1(2) = n(2);

	n = es.eigenvectors().col(1).real();
	l2 = 1.0/sqrt(es.eigenvalues()[1].real());
	axis2(0) = n(0);
	axis2(1) = n(1);
	axis2(2) = n(2);

	n = es.eigenvectors().col(2).real();
	l3 = 1.0/sqrt(es.eigenvalues()[2].real());
	axis3(0) = n(0);
	axis3(1) = n(1);
	axis3(2) = n(2);
}
