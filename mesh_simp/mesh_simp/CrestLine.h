/************************************************************************/
/* 
*	CrestLine.h
*		Define a class to save the crest line data.
*
*/
/************************************************************************/

#ifndef CRESTLINE_H
#define CRESTLINE_H

//#ifdef _WIN64
//typedef char TCHARS1;
//#else
//typedef char TCHARS2;
//#endif

#include <vector>
#include "GeometryElement.h"
#include "global.h"

using namespace std;

typedef pair<pair<VertexID,VertexID>, VertexID> FeatureEdgeType;

class CrestLine
{
public:
	CrestLine();
	~CrestLine();

	void clearAll();

	bool readCrestLineFile(char *infilename);
	
	// Set thresholds
	void set_ridgeness_threshold(double thredVal){ridgeness_threshold = thredVal;}
	void set_sphericalness_threshold(double thredVal){sphericalness_threshold = thredVal;}
	void set_cyclideness_threshold(double thredVal){cyclideness_threshold = thredVal;}

	// Update crest lines based on the updated thresholds
	void update_crestlines();

	inline bool is_feature(VertexID id){return is_feature_flag->at(id);}
	inline float* getVertex(VertexID id){return vertices->at(id);}
	inline float* getNormal(VertexID id){return edgeNormals->at(id);}

public:
	vector<float*> *vertices;	// crest line VERTICES
	vector<FeatureEdgeType> *edges;		// crest line edge
	vector<float*> *edgeNormals;	// normals for each edge
	vector<int> *componentID;	// crest line connected component id for each vertex
	vector<float> *ridgeness;	// ridgeness value for each connected component
	vector<float> *sphericalness;// sphericalness value for each connected component
	vector<float> *cyclideness;	// cyclideness value for each connected component
	vector<bool> *is_feature_flag;	// flag for if it is a feature or not
	int vertexNumber;			// vertex number
	int edgeNumber;				// edge number
	int componentNumber;		// crest line connected component number
	double ridgeness_threshold;
	double sphericalness_threshold;
	double cyclideness_threshold;
};

#endif