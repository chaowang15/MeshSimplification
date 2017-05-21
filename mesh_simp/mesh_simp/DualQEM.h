/************************************************************************/  
/* 
* DualQEM.h
*	Define DualQEM, a quadrics object which can be used in QEM, but here 
*	I only use this class for computing PCA axis of a cluster, given this 
*	cluster's quadrics.
* 
* Author: Chao Wang
* Date: 2.25.2013
*/  
/************************************************************************/

#ifndef DUAL_QEM
#define DUAL_QEM
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
using namespace Eigen;
#include <iostream>
using namespace std;

//extern double d_ratio;
class DualQEM
{
public:
	Matrix3d co_var;
	Vector3d center;
	float n_p;
	double v1;
	double v2;
	double v3;
	Vector3d e1;
	Vector3d e2;
	Vector3d e3;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
	DualQEM();
	DualQEM(const float x_coor, const float y_coor, const float z_coor);
	double evaluate(const Vector3d& v) const;
	double operator()(const Vector3d& v) const { return evaluate(v); }
	double printratio(){cout << v1/v2 << endl;   return v1/v2;}
	//double printregularity() {cout << d_ratio * sqrt(v2 * v1); return d_ratio * sqrt(v2 * v1);}
	double printerror(){cout << v3; return v3;}
	DualQEM& operator=(const DualQEM& Q)
	{ 
		co_var = Q.co_var; center = Q.center; n_p = Q.n_p;
		v1 = Q.v1; v2 = Q.v2; v3 = Q.v3;
		return *this; 
	}

	DualQEM& operator+=(const DualQEM& Q)
	{ 
		Vector3d ncenter = (center*n_p + Q.center * Q.n_p)/(n_p+Q.n_p);
		Vector3d v = ncenter - center;
		Matrix3d m = v * v.transpose();
		co_var += Q.co_var; 
		co_var += n_p * m;
		v = ncenter - Q.center;
		m = v * v.transpose();
		co_var += Q.n_p * m;

		center = ncenter;
		n_p += Q.n_p; 

		invalidateev();
		return *this; 
	}

	DualQEM& operator-=(const DualQEM& Q)
	{ 
		Vector3d ncenter = (center*n_p - Q.center * Q.n_p)/(n_p-Q.n_p);
		Vector3d v = ncenter - center;
		Matrix3d m = v * v.transpose();
		co_var -= Q.co_var; 
		co_var -= (n_p-Q.n_p) * m;
		v = Q.center - center;
		m = v * v.transpose();
		co_var -= Q.n_p * m;

		center = ncenter;
		n_p -= Q.n_p; 
		assert(n_p >= 0);

		invalidateev();
		return *this; 
	}
	//DualQEM& operator*=(double s)
	//{ co_var*=s; return *this; }

	void clear() { co_var.setZero(); center.setZero(); n_p = 0; invalidateev();}
	void invalidateev() {v1 = v2 = v3 = DBL_MAX; e1.setZero(); e2.setZero(); e3.setZero();}

	bool optimize(Vector3d& v);
	//double getmin(Vector3d& v, int type) const;
	//void isospace(MxVector& v);
	//void decompose(MxVector& axis1, MxVector& axis2, MxVector& axis3, double& l1, double& l2, double& l3);
};

#endif
