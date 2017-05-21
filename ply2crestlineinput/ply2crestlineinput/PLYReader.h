#ifndef PLYREADER_H
#define PLYREADER_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "GeoCommon.h"
#include <map>
#include <algorithm>
#include <vector>

using namespace std;

const char typeName[11][10] = {"char","uchar","short","ushort","int","uint","float","double", "uint8", "float32", "int32"} ;
const int typeSize[11] = {1,1,2,2,4,4,4,8,1,4,4} ;
const int totalTypes = 11 ;

class PLYReader
{
public:
	PLYReader( char* fname );
	~PLYReader();

	// Write a ply file (in ascii mode or binary mode)
	void writePLY(char *fname);

	// Write a ply2 file (only in ascii mode)
	void writePLY2(char *fname);

	// Write a gts file (only in ascii mode)
	void writeGTS(char *fname);

	// Write a txt file for the input of detect-crest-line code (only in ascii mode)
	void writeCrestLineInputTxt(char *fname);

	// Flip the 4 bits of x into the inverse order
	// Used for saving faces in binary_big_endian Mode
	void flipBits32 ( void *x );

	void printInfo( ){
		printf("  Vertices: %d Polygons: %d\n",this->numVerts, this->numTrians) ;
	}

	int getNumTriangles(){
		return this->numTrians ;
	}

	int getNumVertices(){
		return this->numVerts ;
	}

	float* getVertex(int ind){
		return &(this->vertarray[ 3 * ind ]) ;
	}

	int* getTriangle(int ind){
		return &(this->faceArray[3*ind]);
	}

public:
	FILE* fin;						// File handle
	int numTrians, curtrian;		// Number of triangles and current triangle
	int numVerts, curvert;			// Number of vertices and current vertex
	float* vertarray ;				// original Vertex array
	int* faceArray;					// original face array
	int countDataSize;				// Data size of count type for a face
	int extraVertex ;				// Number of extra vertex
	int mode ;						// PLY Mode: 0 for ascii, 1 for little-endian, 2 for big-endian (uncommon)
	int vmode ;						// Vertex Mode
	long offset ;					// Offset for reading faces
};

#endif