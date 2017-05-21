/************************************************************************/  
/* 
* glwidget.h
*	Define GLWidget, an inherit class from QGLWidget for rendering models in Qt
* 
* Author: Chao Wang
* Date: 2.25.2013
*/  
/************************************************************************/

#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QtGui/QtGui>
#include <QFileDialog>
#include "global.h"
#include "glwidgetsettings.h"
#include "acvt_mesh.h"

enum ProgramType {BLANK = 0, QEM, QEMOPTIMIZATION, SPLITCLUSTERS, REQEM, MANIFOLDTOPOLOGY, NEWCONNECTIVITY};

class GLWidget : public QGLWidget
{
	Q_OBJECT

public:
	GLWidgetSettings *settings;
	int frameIndex;
	bool flag_isOpened;
	QPoint mouseLastPos;
	QPoint mousecurrentPos;
	float lastMousePos3F[3];
	MouseAtionMode mouseMode;
	float range;
	ACVT *acvt;
	ProgramType programType;
	uint final_vertex_number;

	/* flags for rendering */
	bool show_crest_line_flag;
	bool show_clusters_flag;
	uint connectivity_flag;	/* 0 denotes original mesh connectivity, 
							   1 denotes QEM connectivity, 
							   2 denotes the new cluster connectivity */
	bool show_vertex_weight_flag;

public:
	GLWidget(QWidget *parent = 0);
	~GLWidget();

	void newBlank();
	void clearAll();
	void renderModel();
	inline void show_crest_line(bool flag){show_crest_line_flag = flag;updateGL();}
	inline void show_clusters(bool flag){show_clusters_flag = flag;updateGL();}
	inline void show_vertex_weights(bool flag){show_vertex_weight_flag = flag;updateGL();}
	void run_QEM(uint vertexNum, double reductionRatio = 0.0);
	void run_QEM_Optimization();
	void run_QEM_Optimization_SwapOnce();
	void run_split_clusters();
	void run_create_connectivity();
	void run_ReQEM_after_split(uint vertexNum, double reductionRatio);
	void run_split_nonmanifold_clusters();
	void run_reverse_normals(){acvt->reverse_current_normals(); updateGL();}
	void run_compute_vertex_weight();
	void run_all_steps(uint vertexNum, double reductionRatio);
	void run_compute_vertex_quality();

public slots:
	void openFile();
	void saveFile();
	void setRenderModePoints(){settings->meshDisplayMode = VERTICES;updateGL();}
	void setRenderModeEdges(){settings->meshDisplayMode = EDGES;updateGL();}
	void setRenderModeHiddenLines(){settings->meshDisplayMode = HIDDENLINES;updateGL();}
	void setRenderModeFlat(){settings->meshDisplayMode = FLAT;updateGL();}
	void setRenderModeSmooth(){settings->meshDisplayMode = SMOOTH;updateGL();}
	void setRenderModeTextures(){settings->meshDisplayMode = TEXTURES;updateGL();}
	void setRenderModeNothing(){settings->meshDisplayMode = NOTHING;updateGL();}
	void set_connectivity_flag(uint flag){connectivity_flag = flag;updateGL();}
	void setRidgenessThreshold(int t){
		if (show_crest_line_flag)
		{
			double threval = ((double)(t))/10;
			acvt->ridge->set_ridgeness_threshold(threval);
			acvt->ridge->update_crestlines();
			acvt->ravine->set_ridgeness_threshold(threval);
			acvt->ravine->update_crestlines();
			cout<<"Set Ridgeness Threshold = "<<threval<<endl;
			updateGL();
		}
	}
	void setSphericalnessThreshold(int t){
		if (show_crest_line_flag)
		{
			double threval = ((double)(t))/10;
			acvt->ridge->set_sphericalness_threshold(threval);
			acvt->ridge->update_crestlines();
			acvt->ravine->set_sphericalness_threshold(threval);
			acvt->ravine->update_crestlines();
			cout<<"Set Sphericalness Threshold = "<<threval<<endl;
			updateGL();
		}
		
	}
	void setCyclidenessThreshold(int t){
		if (show_crest_line_flag)
		{
			double threval = ((double)(t))/10;
			acvt->ridge->set_cyclideness_threshold(threval);
			acvt->ridge->update_crestlines();
			acvt->ravine->set_cyclideness_threshold(threval);
			acvt->ravine->update_crestlines();
			cout<<"Set Cyclideness Threshold = "<<threval<<endl;
			updateGL();
		}	
	}

protected:
	void initializeGL();
	void paintGL();
	void resizeGL(int width, int height);
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void resizeEvent(QResizeEvent* event);
	void mouseReleaseEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);

private:
	void initLighting();
	// Convert a 2D coordinate in the image into a 3D coordinate in 3D space
	void point2Triple(int x, int y, float triplef[3]);
	void displayPoints();
	void displayEdges();
	void displayHiddenLinesRemoval();
	void displayTriangles(GLenum disType);
	void displayTextures();
	void displayClusters();
	void displayResultPoints();
	void displayResultEdges_QEM();
	void displayResultEdges_Clusters();
	void displayResultTriangles_QEM(GLenum disType);
	void displayResultTriangles_Clusters(GLenum disType);
	void displayResultClustersInTriangles(GLenum disType);
	void displayResultTriangles_Clusters_ColorInverseNormalFace(GLenum disType);
	void displayHalfEdges();
	void displayVertexWeightDistributions();
//	void displayTextureTriangles(Mesh *inputMesh, GLenum disType = GL_SMOOTH);
};

#endif