/************************************************************************/  
/* 
* global.h
*	Define global parameters such as global constants, type-definitions, etc
* 
* Author: Chao Wang
* Date: 2.25.2013
*/  
/************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

enum MouseAtionMode {NONE = 0, SELECT, ROTATE, TRANSLATE};
enum MeshDisplayType {VERTICES=0, EDGES, HIDDENLINES, FLAT, SMOOTH, TEXTURES, NOTHING};
enum DistortionType {MSE = 0, KGERROR, RIBBON, HOSDORFF, STEP};

const double PI = 3.141592653589793238462643;
const double penalty_factor = 100000;

typedef unsigned int uint;
typedef uint VertexID;
typedef uint EdgeID;
typedef uint FaceID;
typedef int VertexType;

//#ifndef NORMAL
//#define NORMAL
//#endif


/************************************************************************/
/* 
* Crest line parameters
*/
#ifndef CRESTLINE
#define CRESTLINE
#endif

#ifndef BORDERS
#define BORDERS
#endif
// parameter for computing the crest-line weight
// const uint yita = 5;
// Crest line threshold parameters
//extern double ridgeness_threshold;
//extern double sphericalness_threshold;
//extern double cyclideness_threshold;

// parameter for the vertex quality weight
// const uint garma = 5;

// Original QEM Face list 
#ifndef QEMCONNECTIVITY
#define QEMCONNECTIVITY
#endif

/************************************************************************/
/* 
* Cluster parameters
*/
const uint color_number = 17;
const float colors[color_number][3] = {
	//{1.0, 0.0, 0.0},
	{1.0, 0.4, 1.0},
	{0.0, 0.1, 0.0},
	{0.0, 0.0, 1.0},
	{0.0, 1.0, 1.0},
	{1.0, 1.0, 0.0},	
	{0.0, 0.5, 0.0},
	{0.6, 0.0, 1.0},
	{1.0, 1.0, 0.5},
	{1.0, 0.0, 1.0},	
	{0.6, 0.3, 0.0},
	{1.0, 0.0, 0.6},
	{0.5, 0.5, 0.5},
	{1.0, 0.5, 0.5},
	{0.2, 0.6, 0.3},
	{1.0, 0.5, 0.0},
	{0.0, 0.4, 1.0},
	{0.7, 1.0, 1.0}
};

// Used for checking if a normal direction is reversed
const float thresVal_normal_direction = 1e-6;

// Isotropic coefficient for each quadrics
const double isotropic_coefficient = 0.1;

// epsilon: a constant parameter for computing crest line weight or vertex quality weight
const double epsilon_const = 0.1;

// eta: a constant parameter for computing crest line weight
//const double eta_const = 3;

// gamma: a constant parameter for computing vertex quality weight
//const double gamma_const = 7;

// A weight threshold for computing vertex quality weight
// const double weight_threshold_const = 0.1;

/************************************************************************/
#endif