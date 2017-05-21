
#include <math.h>
#include "GeometryElement.h"

Vertex operator + (Vertex &v0, Vertex &v1)
{return Vertex(v0.coord[0]+v1.coord[0], v0.coord[1]+v1.coord[1], 
v0.coord[2]+v1.coord[2]);}

Vertex operator - (Vertex &v0, Vertex &v1)
{return Vertex(v0.coord[0]-v1.coord[0], v0.coord[1]-v1.coord[1], 
v0.coord[2]-v1.coord[2]);}

void Vertex::normalizeCoordinates()
{
	float sum = sqrt(coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2]);
	coord[0] /= sum;
	coord[1] /= sum;
	coord[2] /= sum;
}

void Vertex::normalizeNormal()
{
	float sum = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
	norm[0] /= sum;
	norm[1] /= sum;
	norm[2] /= sum;
}

void Triangle::normalizeNormal()
{
	float sum = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
	norm[0] /= sum;
	norm[1] /= sum;
	norm[2] /= sum;
}
