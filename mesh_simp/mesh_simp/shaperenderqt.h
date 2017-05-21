/************************************************************************/  
/* 
* shaperenderqt.h
*	Define ShapeRenderQt, the main Qt class for the interface window
* 
* Author: Chao Wang
* Date: 2.25.2013
*/  
/************************************************************************/

#ifndef SHAPERENDERQT_H
#define SHAPERENDERQT_H

#include <QtWidgets/QMainWindow>
#include "ui_shaperenderqt.h"
#include "glwidget.h"
#include <QScrollArea>
#include <QGridLayout>
//#include <QDockWidget>

class QScrollArea;
class GLWidget;
class QGridLayout;
//class QDockWidget;

class ShapeRenderQt : public QMainWindow
{
	Q_OBJECT

public:
	ShapeRenderQt(QWidget *parent = 0);
	~ShapeRenderQt();


public slots:
	void closeWindow();
	void newBlank();
	void checkbox_showcrestlines_StateChanged(int a);
	void checkbox_showclusters_StateChanged(int a);
	void pushbutton_QEM();
	void pushbutton_QEM_Optimization();
	//void pushbutton_QEM_Optimization_SwapOnce();
	void pushbutton_SplitClusters();
	void pushbutton_CreateConnectivity();
	void pushbutton_ReQEM();
	void pushbutton_RunReverseNormals();
	//void comboBox_Connectivity_Changed(int a);
	void pushbutton_SplitNonManifold();
	void checkbox_usecrestlines_StateChanged(int a);
	void checkbox_usevertexquality_StateChanged(int a);
	void pushbutton_ComputeVertexWeight();
	void pushbutton_AllSteps();
	void pushbutton_ComputeVertexQuality();

private:
	Ui::ShapeRenderQtClass ui;
	GLWidget *originalGLWidget;
	QScrollArea *originalGLWidgetArea;// to bound originalGLWidget
	QGridLayout *centralLayout;
	//QDockWidget *qemDockWidget;
};

#endif // SHAPERENDERQT_H
