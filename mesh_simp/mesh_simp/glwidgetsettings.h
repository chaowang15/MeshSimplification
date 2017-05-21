/************************************************************************/  
/* 
* glwidgetsettings.h
*	Define GLWidgetSettings containing related parameters for rendering models
*	in Qt, such as rendering colors, line sizes, feature values, etc.
* 
* Author: Chao Wang
* Date: 2.25.2013
*/  
/************************************************************************/

#ifndef GLWIDGETSETTINGS_H
#define GLWIDGETSETTINGS_H

#include "global.h"

class GLWidgetSettings
{
public:
	GLWidgetSettings();
	~GLWidgetSettings();

	void init();
	void clear();
	void reset();
	inline void setMeshDisplayMode(MeshDisplayType type){ meshDisplayMode = type;}

public:
	float point_ambient[4], point_diffuse[4], point_shinness;
	float line_ambient[4], line_diffuse[4], line_shinness;
	float face_ambient[4], face_diffuse[4], face_shinnes;
	float border_ambient[4], border_diffuse[4], border_shinnes;
	float other_ambient[4], other_diffuse[4], other_shinnes;
	float visibility_ambient[4], visibility_diffuse[4], visibility_shinnes;
	float feature_ambient[4], feature_diffuse[4], feature_shinness;
	float feature_line_size;

	float transMatrix[4][4];
	float xShift, yShift, zShift;
	float clipping;

	MeshDisplayType meshDisplayMode; 	// Mesh Viewmode
	float lightBrightness;
	float pointSize, lineSize;

	bool enableLighting, enableBoundingBox; // Flags
	bool textureFlag;
	bool enableOriWgtDescription, enableRstWgtDescription;
};

#endif