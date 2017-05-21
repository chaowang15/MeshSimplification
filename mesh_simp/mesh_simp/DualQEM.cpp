#include "DualQEM.h"

DualQEM::DualQEM()
{
	clear();
}

DualQEM::DualQEM(const float x_coor, const float y_coor, const float z_coor)
{
	co_var.setZero();
	center = Vector3d(x_coor, y_coor, z_coor);
	n_p = 1;
	invalidateev();
}

double DualQEM::evaluate(const Vector3d& v) const
{
	return v.transpose() * co_var * v;
}

bool DualQEM::optimize(Vector3d& v)
{
	SelfAdjointEigenSolver<Matrix3d> es(co_var);
	int smallest = 0;
	int largest = 0;
	int medium = 0;
	for (int i=1; i<3; i++)
	{
		if (abs(es.eigenvalues()[i] > abs(es.eigenvalues()[largest])))
			largest = i;
		if (abs(es.eigenvalues()[i]) < abs(es.eigenvalues()[smallest]))
			smallest = i;
	}
	for (int i=0; i<3; i++)
	{
		if (i != smallest && i!= largest)
		{
			medium = i;
			break;
		}
	}

	//assign values
	v1 = es.eigenvalues()[largest];
	v2 = es.eigenvalues()[medium];
	v3 = es.eigenvalues()[smallest];
	e1 = es.eigenvectors().col(largest).real();
	e2 = es.eigenvectors().col(medium).real();
	e3 = es.eigenvectors().col(smallest).real();
	
	/////////////////////////////////////////////////////////////////
	// Return the eigenvector corresponding to the largest eigenvalue
	if (v3<0)
		v3 = -v3;
	v = e1;
	//if (v3<-1e-15)
	//{
	//	cout <<v1 << " " << v2 << " " << v3 << endl;
	//	cout << co_var << endl;
	//}
	return true;
}
//
//double DualQEM::getmin(Vector3d& v, int type) const
//{
//	//SHAPE
//	assert(v == e3);
//	if (type == 0)
//		return v3;
//	//COMBINATION
//	//isotropic metric is added
//	else if (type == 1)
//	{
//		//assert(v3>=-1e-15);
//		//assert(co_var(0,0) >= -1e-15);
//		//assert(co_var(1,1) >= -1e-15);
//		//assert(co_var(2,2) >= -1e-15);
//		return v3 + d_ratio*(co_var(0,0) + co_var(1,1) + co_var(2,2));
//	}
//
//	////COMBINATION
//	////anisotropic metric is added
//	//else if (type == 1)
//	//{
//	//	assert(v3>=-1e-15);
//	//	assert(co_var(0,0) >= -1e-15);
//	//	assert(co_var(1,1) >= -1e-15);
//	//	assert(co_var(2,2) >= -1e-15);
//
//	//	if(abs(v3) < 1e-10)
//	//	{
//	//		return v3;
//	//	}
//	//	else
//	//	{
//	//		return v3 + d_ratio * sqrt(v2 * v1);
//	//	}
//	//}		
//}
