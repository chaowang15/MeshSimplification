/************************************************************************/  
/* 
* PLYReader.h
*	Define PLYReader for reading and writing ply models
* 
* Author: Chao Wang
* Date: 2.25.2013
*/  
/************************************************************************/

#ifndef PLYREADER_H
#define PLYREADER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kdtree.h"
#include <vector>
#include "PLYWriter.h"
#include <algorithm>
#include "ALGraph.h"
#include "tools.h"

const char typeName[11][10] = {"char","uchar","short","ushort","int","uint","float","double", "uint8", "float32", "int32"} ;
const int typeSize[11] = {1,1,2,2,4,4,4,8,1,4,4} ;
const int totalTypes = 11 ;

using namespace std;

// Compare function
inline bool compPair(const pair<int, int>& x, const pair<int, int>& y)  
{  
	return (x.second > y.second);
}  


class PLYReader
{
public:

	FILE* fin;						// File handle
	int numTrians, curtrian;		// Number of triangles and current triangle
	int numVerts, curvert;			// Number of VERTICES and current vertex
	float* vertArray ;				// original Vertex array
	int* faceArray;					// original face array
	float min[3], max[3] ;			// Bounding box
	float rawmin[3], rawmax[3] ;	// 
	float maxsize ;					// Size of box
	int countDataSize;				// Data size of count type for a face
	int extraVertex ;				// Number of extra vertex
	int temparray[64] ;				// Temporary array for non-triangle faces
	int tottri, curtri ;
	int mode ;						// PLY Mode: 0 for ascii, 1 for little-endian, 2 for big endian
	int vmode ;						// Vertex Mode
	long offset ;					// Offset for reading faces
	char *modelname;				// model name
	float ctr[3];					// center coordinate

	/* Variables for new mesh (i.e. result mesh) */
	vector<float*> *newVertices;			// New Vertex vector
	vector<int*> *newFaces;					// New face vector
	vector<bool> *vtxFlag;					// flag for each vertex
	vector<vector<float>> *vtxNewExtraData; // New Vertex Extra data

	/* Variables for original mesh */
	vector<vector<float>> vtxExtraData;		// original extra data for each vertex (such as vertex quality data, or RGB colors)
	vector<char*> vtxExtraNames;			// names of the extra data (used in saving meshes)
	vector<int> *vtxExtraDataIndex;			// extra data index in original mesh for each vertex in result mesh
											// (Used in transferring original vertex extra data to result mesh)
	
	int extraTypeNumber;					// number of types for extra data for each vertex
											// 0 means there is no extra data  for each vertex (only x,y,z coordinates); 
											// larger than 0 means that there are extra data for each vertex with type number this 
											// parameter. For instance, if there is a vertex quality (color) data for each vertex, 
											// this parameter is 1.
	
public:
	PLYReader( char* fname, float scl );
	~PLYReader();

	// get center
	inline float* getCenter(){return ctr;}

	// Clear the newly created vertex vector and face vector
	void clearNewVertexFaceVector();

	// Delete the duplicate VERTICES in the model
	void DeleteDuplicateVertices();

	// Write a ply
	void writeOriginalPLY(char *fname);

	// Write newly created ply file
	void writeNewPLY(char *fname);

	// Write ply with extra data for each vertex
	void writePLYWithExtraData(char *fname, vector<char*> vecExtraNames, vector<vector<float>> vecExtraData);

	// Write ply with extra data for each vertex
	void writeNewPLYWithExtraData(char *fname);

	// Flip the 4 bits of x into the inverse order
	// Used for saving faces in binary_big_endian Mode
	void flipBits32 ( void *x );

	// Remove the largest difference between two meshes
	void RemoveLargestDifferenceFrom(PLYReader *inputPly);

	// Get the difference (VERTICES and faces) with another input mesh and save the difference in vecVtxIndex and vecFace
	void GetDifferenceWithAnotherMesh(PLYReader *inputPly, vector<int> *vecResultFaceIndex);

	// Get connected components from input faces with indexes in [startIndex,endIndex) (the order is from largest to smallest)
	// Input: 
	//	vecFaceIndex: input face indexes, each face is denoted by three integers
	//	startIndex: the start index of all the components you want to get
	//	endIndex: the end index of all the components you want to get
	// Output:
	//	vecResultVtxIndex: result face indexes containing all the desired connected components
	// Example: GetConnectedComponents(vecFaceIndex, 0, 5, vecResultVtxIndex) means that the connected components with indexes 
	//	from 0 to 4 (from large to small) in the vecFaceIndex will be obtained and output into the vecResultVtxIndex.
	void GetConnectedComponents(vector<int> *vecFaceIndex, int startIndex, int endIndex, vector<int> *vecResultVtxIndex);

	// Get connected components from input faces with number of VERTICES smaller than vertexThreshold
	// Input: 
	//	vecFaceIndex: input face indexes, each face is denoted by three integers
	//	vertexThreshold: vertex number threshold
	// Output:
	//	vecResultVtxIndex: result face indexes containing all the desired connected components
	// Example: GetConnectedComponents(vecFaceIndex, 100, vecResultVtxIndex) means that the connected components with point 
	//	number smaller than 100
	void GetSmallConnectedComponents(vector<int> *vecFaceIndex, int vertexThreshold, vector<int> *vecResultVtxIndex);

	// Get floating VERTICES
	void GetfloatingVertices(vector<int> *vecVtx);

	// Delete VERTICES saved in vecDeletedVtx and save the rest VERTICES and faces in newVertices and newFaces
	void DeleteVertices(vector<int> *vecDeletedVtx);

	// Compute Euclidean distance between each point in original model and the nearest point in the input model
	void computeNearestEuclideanDistance(PLYReader *inputPly, vector<double> *vecDis);

	// Get the extra data information for each vertex in the new mesh
	void getExtraDataForNewMesh(PLYReader *inputPly, vector<double> *vecDis);

	// Compute center
	void computeCenter();

	void printNewModelInfo()
	{
		printf("Output Model:\n VERTICES: %d Polygons: %d\n", newVertices->size(),newFaces->size()) ;
	}

	void reset( )
	{
		curtrian = 0 ;
		curvert = 0 ;
		fseek( fin, offset, SEEK_SET ) ;
	}

	void printInfo ( )
	{
		printf("  VERTICES: %d Polygons: %d\n",this->numVerts, this->numTrians) ;
	}

	// Get bounding box
	float getBoundingBox ( float origin[3] )
	{
		for ( int i = 0 ; i < 3 ; i ++ )
			origin[i] = min[i] ;
		return maxsize ;
	};

	inline float getRange(){return maxsize;}

	void getRawBoundingBox ( float low[3], float high[3] )
	{
		for ( int i = 0 ; i < 3 ; i ++ )
		{
			low[i] = rawmin[i] ;
			high[i] = rawmax[i] ;
		}
	}

	// Get storage size
	int getMemory () {
		return sizeof ( class PLYReader ) + this->numVerts * sizeof( float ) * 3 ;	
	};

	void printBits32 ( void *x ){
		unsigned char *temp = (unsigned char *)x;
		printf("%2x %2x %2x %2x\n", temp[0], temp[1], temp[2], temp[3]) ;
	}

	int getNumTriangles(){
		return this->numTrians ;
	}

	int getNumVertices(){
		return this->numVerts ;
	}

	float* getVertex( int ind ){
		return &(this->vertArray[ 3 * ind ]) ;
	}

	int* getTriangle(int ind){
		return &(this->faceArray[3*ind]);
	}

	int getMode(){
		return this->mode ;
	}

	// Get vertex extra data (i.e. vertex quality data)
	// Note: extraDim = 0 by default
	float getVtxExtraData(int vtxIndex, int extraDim=0){
		return vtxExtraData[extraDim][vtxIndex];}

	// Get maximum and minimum value of extra data
	// Note: extraDim is the dimension of extra data (usually 0)
	void getVtxMaxMinExtraData(int extraDim, float &maxVal, float &minVal){
		maxVal = *max_element(vtxExtraData[extraDim].begin(), vtxExtraData[extraDim].end());
		minVal = *min_element(vtxExtraData[extraDim].begin(), vtxExtraData[extraDim].end());
	}
};


#endif