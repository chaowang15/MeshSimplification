/************************************************************************/
/* 
** GeometryElement.h
** Class Definition of Some Fundamental geometry elements
** This file includes the following classes:
** (1) Vertex
** (2) Triangle
** (3) Vector
** (4) Edge
*/
/************************************************************************/


#ifndef GEOMETRYELEMENT_H
#define GEOMETRYELEMENT_H

/*
**  Class Vertex
*/
class Vertex
{
public:
	Vertex()
	{
		for (int i = 0; i < 3; ++i)
		{
			coord[i] = 0;
			norm[i] = 0;
		}
		id = -1;
		txtVal = 0;
	}

	Vertex(float *cor)
	{
		for (int i = 0; i < 3; ++i)
		{
			coord[i] = cor[i];
			norm[i] = 0;
		}
		id = -1;
		txtVal = 0;
	}

	Vertex(Vertex *v)
	{
		for (int i = 0; i < 3; ++i)
		{
			coord[i] = v->coord[i];
			norm[i] = v->norm[i];
		}
		id = v->id;
		txtVal = 0;
	}

	Vertex(float cor0, float cor1, float cor2, int idNo = -1)
	{
		coord[0] = cor0;
		coord[1] = cor1;
		coord[2] = cor2;
		id = idNo;
		norm[0] = norm[1] = norm[2] = 0;
		txtVal = 0;
	}

	friend Vertex operator + (Vertex &v0, Vertex &v1);
	friend Vertex operator - (Vertex &v0, Vertex &v1);

	// normalize coordinates and normal vector
	void normalizeCoordinates();
	void normalizeNormal();

public:
	float coord[3];	// coordinates
	float norm[3];	// normal
	int id;			// id index of vertex
	float txtVal;	// texture value
};


/*
**  Class Triangle
*/
class Triangle
{
public:
	Triangle()
	{
		for (int i = 0; i < 3; ++i)
		{
			vtx[i] = 0;
			norm[i] = 0;
		}
		id = -1;
	}

	Triangle(int *t)
	{
		for (int i = 0; i < 3; ++i)
		{
			vtx[i] = t[i];
			norm[i] = 0;
		}
		id = -1;	
	}

	Triangle(int v0, int v1, int v2, int idNo = -1)
	{
		vtx[0] = v0;
		vtx[1] = v1;
		vtx[2] = v2;
		id = idNo;
		norm[0] = norm[1] = norm[2] = 0;
	}

	Triangle(Triangle *t)
	{
		for (int i = 0; i < 3; ++i)
		{
			vtx[i] = t->vtx[i];
			norm[i] = t->norm[i];
		}
		id = t->id;
	}

	void copyTriangle(Triangle *t)
	{
		for (int i = 0; i < 3; ++i)
		{
			vtx[i] = t->vtx[i];
			norm[i] = t->norm[i];
		}
		id = t->id;
	}

	// normalize normal vector
	void normalizeNormal();

public:
	int vtx[3];		// vertex indexes
	float norm[3];	// normal vector
	int id;			// id index of triangle
};


/*
**  Class Vector
*/
class Vector
{
public:
	Vector()
	{
		for(int i = 0; i < 3; ++i)
			coord[i] = 0;
	}

	Vector(float *v)
	{
		for(int i = 0; i < 3; ++i)
			coord[i] = v[i];
	}
	Vector(float *v0, float *v1)
	{
		for(int i = 0; i < 3; ++i)
			coord[i] = v1[i] - v0[i];
	}
	Vector(float v0, float v1, float v2)
	{
		coord[0] = v0; coord[1] = v1; coord[2] = v2;
	}
	Vector(Vertex *v0, Vertex *v1)
	{
		for(int i = 0; i < 3; ++i)
			coord[i] = v1->coord[i] - v0->coord[i];
	}
	~Vector(){}

public:
	float coord[3];		// vector coordinates
};

// Edge
class Edge
{
public:
	Edge(int i, int j, int id = -1){
		node[0]=i;
		node[1]=j;
		faceId=id;
	}
	Edge(const Edge &e){
		node[0] = e.node[0]; 
		node[1] = e.node[1];
		faceId = e.faceId;
	}
	Edge(Edge *e){
		node[0] = e->node[0]; 
		node[1] = e->node[1];
		faceId = e->faceId;
	}
	~Edge(){}

	// Compare edge
	bool compareEdge(Edge *e)
	{
		if ((node[0] == e->node[0] && node[1] == e->node[1]) || 
			(node[0] == e->node[1] && node[1] == e->node[0]))
			return true;
		else
			return false;
	}

public:
	int node[2];	// two nodes connecting the edge
	int faceId;		// index of the face that this edge belongs to
};



#endif