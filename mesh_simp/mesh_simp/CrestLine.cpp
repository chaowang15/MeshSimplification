#include "CrestLine.h"
#include "tools.h"
#include <fstream>

CrestLine::CrestLine()
{
	vertexNumber = edgeNumber = componentNumber = 0;
	vertices = new vector<float*>;
	edges = new vector<FeatureEdgeType>;
	componentID = new vector<int>;
	ridgeness = new vector<float>;
	sphericalness = new vector<float>;
	cyclideness = new vector<float>;
	is_feature_flag = new vector<bool>;
	edgeNormals = new vector<float*>;
	ridgeness_threshold = sphericalness_threshold = cyclideness_threshold = 0;
}

CrestLine::~CrestLine()
{
	clearAll();
}

void CrestLine::clearAll()
{
	if (!vertices->empty())
	{
		for (vector<float*>::iterator iter = vertices->begin(); iter != vertices->end(); ++iter)
		{
			delete [](*iter);
		}
		vertices->clear();
		ClearVector(*vertices);
	}
	if (!edges->empty())
	{
		edges->clear();
		ClearVector(*edges);
	}
	if (!edgeNormals->empty())
	{
		for (vector<float*>::iterator iter = edgeNormals->begin(); iter != edgeNormals->end(); ++iter)
		{
			delete [](*iter);
		}
		edgeNormals->clear();
		ClearVector(*edgeNormals);
	}
	componentID->clear();
	ClearVector(*componentID);
	ridgeness->clear();
	ClearVector(*ridgeness);
	sphericalness->clear();
	ClearVector(*sphericalness);
	cyclideness->clear();
	ClearVector(*cyclideness);
	is_feature_flag->clear();
	ClearVector(*is_feature_flag);
}

bool CrestLine::readCrestLineFile(char *infilename)
{
	ifstream infile;
	infile.open(infilename, ios::in);
	if (!infile.good())
	{
		cout<<"Can not open the file \""<<infilename<<"\""<<endl;
		return false;
	}
	infile>>vertexNumber>>edgeNumber>>componentNumber;
	
	// Read each crest line vertex
	float tempVal = 0.0;
	int id = 0, edgeNode0 = 0, edgeNode1 = 0;
	for (int i = 0; i != vertexNumber; ++i)
	{
		float *v = new float[3];
		for (int j = 0; j != 3; ++j)
		{
			infile>>v[j];
		}
		vertices->push_back(v);
		infile>>id;
		componentID->push_back(id);
	}
	for (int i = 0; i != componentNumber; ++i)
	{
		infile>>tempVal;
		ridgeness->push_back(tempVal);
		infile>>tempVal;
		sphericalness->push_back(tempVal);
		infile>>tempVal;
		cyclideness->push_back(tempVal);
	}
	for (int i = 0; i != edgeNumber; ++i)
	{
		infile>>edgeNode0>>edgeNode1>>id;
		edges->push_back(FeatureEdgeType(pair<VertexID,VertexID>(edgeNode0, edgeNode1), id));
	}
	infile.close();
	return true;
}

void CrestLine::update_crestlines()
{
	is_feature_flag->clear();
	is_feature_flag->resize(vertexNumber, true);
	for (int i = 0; i != vertexNumber; ++i)
	{
		int cid = componentID->at(i);
		if ((ridgeness->at(cid) < ridgeness_threshold)
			|| (sphericalness->at(cid) < sphericalness_threshold)
			|| (cyclideness->at(cid) < cyclideness_threshold))
		{
			is_feature_flag->at(i) = false;
		}
	}
}