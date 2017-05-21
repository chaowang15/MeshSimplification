#include "glwidgetsettings.h"


//////////////////////////////////////////////////////////////////////////
// Class GLWidgetSettings
GLWidgetSettings::GLWidgetSettings()
{
	init();
}

GLWidgetSettings::~GLWidgetSettings()
{
	clear();
}

void GLWidgetSettings::init()
{
	meshDisplayMode = VERTICES;
	pointSize = 5.0;
	lineSize = lightBrightness = 1.5;
	enableBoundingBox = false;
	enableLighting = true;
	textureFlag = false;
	enableOriWgtDescription = true;
	enableRstWgtDescription = false;

	point_ambient[0] = 0.2; point_ambient[1] = 0.2; point_ambient[2] = 0.2; point_ambient[3] = 1.0;
	point_diffuse[0] = 0.8; point_diffuse[1] = 0.8; point_diffuse[2] = 0.8 ; point_diffuse[3] = 1.0;

	line_ambient[0] = 0.2; line_ambient[1] = 0.2; line_ambient[2] = 0.2; line_ambient[3] = 1.0;
	line_diffuse[0] = 0.0; line_diffuse[1] = 0.0; line_diffuse[2] = 0.0; line_diffuse[3] = 1.0; 

	feature_ambient[0] = 1.0; feature_ambient[1] = 0.0; feature_ambient[2] = 0.0; feature_ambient[3] = 1.0; 
	feature_diffuse[0] = 1.0; feature_diffuse[1] = 0.0; feature_diffuse[2] = 0.0; feature_diffuse[3] = 1.0;
	face_shinnes = 30;
	feature_line_size = 10;

	border_ambient[0] = 1.0; border_ambient[1] = 1.0 ; border_ambient[2] = 0.0; border_ambient[3] = 1.0;
	border_diffuse[0] = 1.0; border_diffuse[1] = 1.0; border_diffuse[2] = 0.0 ; border_diffuse[3] = 1.0;

	other_ambient[0] = 0.5; other_ambient[1] = 0.2 ; other_ambient[2] = 0.5; other_ambient[3] = 1.0;
	other_diffuse[0] = 0.4; other_diffuse[1] = 0.5; other_diffuse[2] = 0.6 ; other_diffuse[3] = 1.0;

	visibility_ambient[0] = 0.0; visibility_ambient[1] = 1.0 ; visibility_ambient[2] = 1.0; visibility_ambient[3] = 1.0;
	visibility_diffuse[0] = 0.0; visibility_diffuse[1] = 1.0; visibility_diffuse[2] = 1.0 ; visibility_diffuse[3] = 1.0;

	face_ambient[0] = 0.2; face_ambient[1] = 0.2 ; face_ambient[2] = 0.2 ; face_ambient[3] = 1.0;
	face_diffuse[0] = 0.5; face_diffuse[1] = 1.0; face_diffuse[2] = 0.0; face_diffuse[3] = 1.0;
	point_shinness = 20; line_shinness = 20; face_shinnes = 15; border_shinnes = 20; other_shinnes = 5; visibility_shinnes = 4;

	reset();
}

void GLWidgetSettings::clear()
{
	init();
}

void GLWidgetSettings::reset()
{
	xShift = 0.0;
	yShift = 0.0;
	zShift = 0.0; 
	clipping = 0.0;

	// create identity matrix
	for (int i=0; i < 4; i++)
		for (int j=0; j < 4; j++)
			transMatrix[i][j]= (i == j) ? 1.0 : 0.0;
}