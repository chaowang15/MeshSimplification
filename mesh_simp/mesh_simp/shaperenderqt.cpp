#include "shaperenderqt.h"

ShapeRenderQt::ShapeRenderQt(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	originalGLWidget = new GLWidget;
	originalGLWidgetArea = new QScrollArea(parent);
	originalGLWidgetArea->setWidget(originalGLWidget);
	originalGLWidgetArea->setWidgetResizable(true);
	originalGLWidgetArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	originalGLWidgetArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	originalGLWidgetArea->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	originalGLWidgetArea->setMinimumSize(50, 50);

	centralLayout = new QGridLayout;
	centralLayout->addWidget(originalGLWidgetArea, 0, 0);
	ui.centralWidget->setLayout(centralLayout);

	ui.dockWidget_QEM->showMaximized();

	//qemDockWidget = new QDockWidget(tr("QEM Operations"),this);

	ui.lineEdit_eta->setEnabled(false);
	ui.lineEdit_gamma->setEnabled(false);

	// Connect actions
	connect(ui.action_Exit, SIGNAL(triggered()), this, SLOT(closeWindow()));
	connect(ui.action_New, SIGNAL(triggered()), this, SLOT(newBlank()));
	connect(ui.action_Open, SIGNAL(triggered()), originalGLWidget, SLOT(openFile()));
	connect(ui.action_Save, SIGNAL(triggered()), originalGLWidget, SLOT(saveFile()));
	//connect(ui.action_QEM, SIGNAL(triggered()), originalGLWidget, SLOT(run_QEM()));
	//connect(ui.action_QEM_Optimization, SIGNAL(triggered()), originalGLWidget, SLOT(run_QEM_Optimization()));
	connect(ui.action_Points, SIGNAL(triggered()), originalGLWidget, SLOT(setRenderModePoints()));
	connect(ui.action_Edges, SIGNAL(triggered()), originalGLWidget, SLOT(setRenderModeEdges()));
	connect(ui.action_Hidden_Lines, SIGNAL(triggered()), originalGLWidget, SLOT(setRenderModeHiddenLines()));
	connect(ui.action_Flat, SIGNAL(triggered()), originalGLWidget, SLOT(setRenderModeFlat()));
	connect(ui.action_Smooth, SIGNAL(triggered()), originalGLWidget, SLOT(setRenderModeSmooth()));
	connect(ui.action_Textures, SIGNAL(triggered()), originalGLWidget, SLOT(setRenderModeTextures()));
	connect(ui.checkBox_ShowCrestLine, SIGNAL(stateChanged(int)), this, SLOT(checkbox_showcrestlines_StateChanged(int)));
	connect(ui.checkBox_ShowClusters, SIGNAL(stateChanged(int)), this, SLOT(checkbox_showclusters_StateChanged(int)));
	connect(ui.checkBox_UseCrestLine, SIGNAL(stateChanged(int)), this, SLOT(checkbox_usecrestlines_StateChanged(int)));
	connect(ui.checkBox_UseVertexQuality,  SIGNAL(stateChanged(int)), this, SLOT(checkbox_usevertexquality_StateChanged(int)));
	connect(ui.horizontalScrollBar_Ridgeness, SIGNAL(valueChanged(int)), originalGLWidget, SLOT(setRidgenessThreshold(int)));
	connect(ui.horizontalScrollBar_Sphericalness, SIGNAL(valueChanged(int)), originalGLWidget, SLOT(setSphericalnessThreshold(int)));
	connect(ui.horizontalScrollBar_Cyclideness, SIGNAL(valueChanged(int)), originalGLWidget, SLOT(setCyclidenessThreshold(int)));
	connect(ui.pushButton_QEM, SIGNAL(clicked()), this, SLOT(pushbutton_QEM()));
	connect(ui.pushButton_QEM_Optimization, SIGNAL(clicked()), this,  SLOT(pushbutton_QEM_Optimization()));
	//connect(ui.pushButton_QEM_Optimization_SwapOnce, SIGNAL(clicked()), this,  SLOT(pushbutton_QEM_Optimization_SwapOnce()));
	connect(ui.pushButton_SplitClusters, SIGNAL(clicked()), this,  SLOT(pushbutton_SplitClusters()));
	connect(ui.pushButton_CreateConnectivity, SIGNAL(clicked()), this,  SLOT(pushbutton_CreateConnectivity()));
	connect(ui.pushButton_ReQEM, SIGNAL(clicked()), this,  SLOT(pushbutton_ReQEM()));
	connect(ui.pushButton_SplitNonManifold, SIGNAL(clicked()), this,  SLOT(pushbutton_SplitNonManifold()));
	connect(ui.pushButton_ReverseNormals, SIGNAL(clicked()), this,  SLOT(pushbutton_RunReverseNormals()));
	connect(ui.pushButton_ComputeWeights, SIGNAL(clicked()), this,  SLOT(pushbutton_ComputeVertexWeight()));
	//connect(ui.comboBox_Connectivity, SIGNAL(currentIndexChanged(int)), this, SLOT(comboBox_Connectivity_Changed(int)));
	connect(ui.pushButton_AllSteps, SIGNAL(clicked()), this,  SLOT(pushbutton_AllSteps()));
	//connect(ui.pushButton_ComputeVertexQuality, SIGNAL(clicked()), this,  SLOT(pushbutton_ComputeVertexQuality()));

}

ShapeRenderQt::~ShapeRenderQt()
{
	delete originalGLWidget;
	delete originalGLWidgetArea;
	delete centralLayout;
}

//void ShapeRenderQt::comboBox_Connectivity_Changed(int a)
//{
//	originalGLWidget->set_connectivity_flag(a);
//}

void ShapeRenderQt::checkbox_showcrestlines_StateChanged(int a)
{
	if (a > 0)
	{
		ui.horizontalScrollBar_Cyclideness->setEnabled(true);
		ui.horizontalScrollBar_Ridgeness->setEnabled(true);
		ui.horizontalScrollBar_Sphericalness->setEnabled(true);
		originalGLWidget->show_crest_line(true);
	}
	else
	{
		ui.horizontalScrollBar_Cyclideness->setEnabled(false);
		ui.horizontalScrollBar_Ridgeness->setEnabled(false);
		ui.horizontalScrollBar_Sphericalness->setEnabled(false);
		originalGLWidget->show_crest_line(false);
	}
}

void ShapeRenderQt::checkbox_showclusters_StateChanged(int a)
{
	if (a > 0)
	{
		originalGLWidget->show_clusters(true);
	}
	else
	{
		originalGLWidget->show_clusters(false);
	}
}

void ShapeRenderQt::checkbox_usecrestlines_StateChanged(int a)
{
	if (a > 0)
	{
		ui.lineEdit_eta->setEnabled(true);
	}
	else
	{
		ui.lineEdit_eta->setEnabled(false);
	}
}

void ShapeRenderQt::checkbox_usevertexquality_StateChanged(int a)
{
	if (a > 0)
	{
		ui.lineEdit_gamma->setEnabled(true);
	}
	else
	{
		ui.lineEdit_gamma->setEnabled(false);
	}
}

void ShapeRenderQt::pushbutton_QEM()
{
	QString ratio = ui.lineEdit_ReductionRatio->text();
	QString finalNum = ui.lineEdit_FinalVertexNumber->text();
	originalGLWidget->run_QEM(finalNum.toFloat(), ratio.toFloat());
}

void ShapeRenderQt::pushbutton_CreateConnectivity()
{
	originalGLWidget->run_create_connectivity();
}

void ShapeRenderQt::pushbutton_ReQEM()
{
	QString ratio = ui.lineEdit_ReductionRatio->text();
	QString finalNum = ui.lineEdit_FinalVertexNumber->text();
	originalGLWidget->run_ReQEM_after_split(finalNum.toFloat(), ratio.toFloat());
}

void ShapeRenderQt::pushbutton_QEM_Optimization()
{
	originalGLWidget->run_QEM_Optimization();
}

void ShapeRenderQt::pushbutton_SplitNonManifold()
{
	originalGLWidget->run_split_nonmanifold_clusters();
}

//void ShapeRenderQt::pushbutton_QEM_Optimization_SwapOnce()
//{
//	originalGLWidget->run_QEM_Optimization_SwapOnce();
//}

void ShapeRenderQt::pushbutton_SplitClusters()
{
	originalGLWidget->run_split_clusters();
}

void ShapeRenderQt::pushbutton_RunReverseNormals()
{
	originalGLWidget->run_reverse_normals();
}

void ShapeRenderQt::pushbutton_ComputeVertexWeight()
{
	if (ui.checkBox_UseCrestLine->isChecked())
	{
		originalGLWidget->acvt->setCrestLineFlag(true);
		originalGLWidget->acvt->set_eta(ui.lineEdit_eta->text().toDouble());
	}
	else
	{
		originalGLWidget->acvt->setCrestLineFlag(false);
	}
	if (ui.checkBox_UseVertexQuality->isChecked())
	{
		originalGLWidget->acvt->setVertexQualityFlag(true);
		originalGLWidget->acvt->set_gamma(ui.lineEdit_gamma->text().toDouble());
		originalGLWidget->acvt->set_vertex_weight_threshold(ui.lineEdit_VertexQualityThreshold->text().toDouble());
	}
	else
	{
		originalGLWidget->acvt->setVertexQualityFlag(false);
	}

	originalGLWidget->run_compute_vertex_weight();
}

void ShapeRenderQt::closeWindow()
{
	this->close();
}

void ShapeRenderQt::newBlank()
{
	ui.horizontalScrollBar_Cyclideness->setEnabled(false);
	ui.horizontalScrollBar_Ridgeness->setEnabled(false);
	ui.horizontalScrollBar_Sphericalness->setEnabled(false);
	originalGLWidget->show_crest_line(false);
	ui.checkBox_ShowCrestLine->setChecked(false);
	ui.checkBox_UseCrestLine->setChecked(false);
	ui.checkBox_UseVertexQuality->setChecked(false);
	ui.checkBox_ShowClusters->setChecked(true);
	ui.lineEdit_eta->setEnabled(false);
	ui.lineEdit_gamma->setEnabled(false);
	originalGLWidget->newBlank();
}

void ShapeRenderQt::pushbutton_AllSteps()
{
	QString ratio = ui.lineEdit_ReductionRatio->text();
	QString finalNum = ui.lineEdit_FinalVertexNumber->text();
	originalGLWidget->run_all_steps(finalNum.toFloat(), ratio.toFloat());
}

void ShapeRenderQt::pushbutton_ComputeVertexQuality()
{
	originalGLWidget->run_compute_vertex_quality();
}