/************************************************************************/  
/* 
* tools.h
*	Define a great deal of common global computation functions for use.
* 
* Author: Chao Wang
* Date: 2.25.2013
*/  
/************************************************************************/

#ifndef TOOL_H
#define TOOL_H

#include "math.h"
#include "GeometryElement.h"
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <QVector>
#include <QString>
#include "global.h"

using namespace std;

/********************************************************/
/* 
* Geometric Calculation
*/
// Compute the Euclidean distance between two points
double getLineLength(Vertex *v0, Vertex *v1);
double getLineLength(double ax,double ay,double az,double bx,double by,double bz);

// Cross Product
void crossProduct(Vertex *v0, Vertex *v1, Vertex *v2, double coord[3]);
void crossProduct(Vector *vec1, Vector *vec2, double coord[3]);
// Cross product
template <class T>
inline void crossProduct(T vec1[3], T vec2[3], T vec[3])
{
	T ax = vec1[0], ay = vec1[1], az = vec1[2];
	T bx = vec2[0], by = vec2[1], bz = vec2[2];
	vec[0] = ay * bz - az * by;
	vec[1] = az * bx - ax * bz;
	vec[2] = ax * by - ay * bx;	
}

// Dot Product
template <class T>
T dot3(T a[3], T b[3])
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double dot3(Vector *vec1, Vector *vec2);

// Normalize a vector into unit one
// Normalize a vertex (or a vector) into unit one
template <class T>
inline void normalize(T a[3])
{
	T n = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	if (n > 1e-8)
	{
		a[0] /= n;
		a[1] /= n;
		a[2] /= n;
	}
}

// Add a triangle's normal onto a vertex's normal
void addTriNorm2VtxNorm(Vertex *v, Triangle *t);

// Check if an edge (with two endpoints v0 and v1) is in a triangle
bool IsEdgeInTriangle(int v0, int v1, Triangle *t);

// Compute the length of a vector
double vectorLength(double v[3]);

// compute the RGB color of an input weight value in jet format
template <class T>
void computeRGBcolorOfWeight
(// input
	T minWgt,		// minimum weight 
	T maxWgt,		// maximum weight
	T inputWgt,		// input weight value needed to be colored
	// output
	T clr[4]		// output color RGB values (between 0.0 and 1.0)
)
{
	// Set the jet color first
	int i = 0;
	int r[64] = {0}, g[64] = {0}, b[64] = {0};
	for (i = 0; i <= 7; ++i)
	{
		r[i] = 0;
		g[i] = 0;
		b[i] = 16 *(i+9);
	}
	for (i = 1; i <= 16; ++i)
	{
		r[i+7] = 0;
		g[i+7] = g[i+6] + 16;
		b[i+7] = b[i+6];
	}
	for (i = 1; i <= 16; ++i)
	{
		r[i+23] = r[i+22] + 16;
		g[i+23] = g[i+22];
		b[i+23] = b[i+22] - 16;
	}
	for (i = 1; i <= 16; ++i)
	{
		r[i+39] = r[i+38];
		g[i+39] = g[i+38] - 16;
		b[i+39] = b[i+38];
	}
	for (i = 1; i <= 8; ++i)
	{
		r[i+55] = r[i+54] - 16;
		g[i+55] = 0;
		b[i+55] = 0;
	}

	// Get the color of input weight value
	int index = 0;
	if (maxWgt - minWgt > 1e-5)
	{
		index = (int)(64*(inputWgt-minWgt)/(maxWgt-minWgt));
		if (inputWgt == maxWgt)
			index--;
	}
	clr[0] = (T)r[index] / 255;
	clr[1] = (T)g[index] / 255;
	clr[2] = (T)b[index] / 255;
	clr[3] = 1.0;
}

// Compute the angle (in degree) between two input vectors
double computeAngleOfVectors(Vector *vec1, Vector *vec2);
double computeAngleOfVectors(double vec1[3], double vec2[3]);

// Compute the distance between an input point and an input
// line determined by two input VERTICES
double computeDisBtwnPointAndLine(double p[3], double v0[3], double v1[3]);

//	check if a is included in the three numbers v0, v1 and v2.
//	Input: four integers a, v0, v1, v2
//	Output: if a is not included in v0, v1 and v2, then return -1;
//	    if a is equal to v0, then return 0; or if a is equal to v1,
//	    then return 1; if a is equal to v2, then return 2.
//
template <class T>
int isIncluded(T a, T v0, T v1, T v2)
{
	if (a != v0 && a != v1 && a != v2)
		return -1;
	else if (a == v0)
		return 0;
	else if (a == v1)
		return 1;
	else
		return 2;
}

// Given two vertex indexes in a triangle, find the third vertex index
// If a and b are both in the triangle, then return
// the 3rd vertex, or return -1
int get3rdVtrOfTriangle(Triangle *t, int v0, int v1);

// Sort three unsigned integers
void sort3(unsigned int &v0, unsigned int &v1, unsigned int &v2);

// Sort two unsigned integers
template <class T>
inline void sort2(T &v0, T &v1)
{
	if (v0 > v1)
	{
		T temp = v1;
		v1 = v0;
		v0 = temp;
	}
}

/********************************************************/
/* 
* String Processing
*/
// Given the full path of a file, find its file name 
void getFileName(char *filePath, char *fileName);

// Given the full path of a file, find its suffix (with .)
void getFileSuffix(char *filePath, char *fileSuffix);

// Convert string to QString
QString s2q(const string &s);

// Convert QString to string
string q2s(const QString &s);

// Convert wstring to string
std::string ws2s(const std::wstring& ws);

// Convert string to wstring
std::wstring s2ws(const std::string& s);

// Convert int to QString
QString int2QString(const int &n);

// Convert QString to char * (C-stype string)
char* cstyleString(QString& str);


/********************************************************/
/* 
* Free memory functions
*/
// Really free memory of a vector
template <class T>
inline void ClearVector( vector<T>& vt ) 
{
	vector<T> vtTemp; 
	vtTemp.swap( vt );
}

// Really free memory of a set
template <class T>
inline void ClearSet( set<T>& vt ) 
{
	set<T> vtTemp; 
	vtTemp.swap( vt );
}

#endif