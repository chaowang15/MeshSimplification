#include "glwidget.h"
#include "tools.h"

//////////////////////////////////////////////////////////////////////////
/* 
* Fundamental functions
*/

GLWidget::GLWidget(QWidget *parent) 
	: QGLWidget(parent)
{
	settings = new GLWidgetSettings;
	flag_isOpened = false;
	range = 1.0;
	programType = BLANK;
	show_clusters_flag = show_crest_line_flag = show_vertex_weight_flag = false;
	connectivity_flag = 0;
	final_vertex_number = 0;
}

GLWidget::~GLWidget()
{
	clearAll();
}

void GLWidget::clearAll()
{
	if (settings != NULL)
	{
		delete settings;
	}
	if (flag_isOpened)
	{
		delete acvt;
	}
}

void GLWidget::newBlank()
{
	clearAll();
	settings = new GLWidgetSettings;
	flag_isOpened = false;
	show_clusters_flag = show_crest_line_flag = show_vertex_weight_flag = false;
	connectivity_flag = 0;
	programType = BLANK;
	final_vertex_number = 0;
	cout<<"New blank"<<endl;
	updateGL();
}

//////////////////////////////////////////////////////////////////////////
/*
* OpenGL related functions : protected functions inherited from QGLWidget
*/

void GLWidget::initLighting()
{
	const GLfloat diffuseLight[] = {1.0, 1.0, 1.0, 1.0};
	const GLfloat ambientLight[] = {0.1, 0.1, 0.1, 1.0};
	const GLfloat lightPositions[4][4] = 
	{
		{1.0, 1.0,  0.0, 0.0},
		{1.0, 1.0, -1.0, 0.0},
		{-0.1, -0.1,  1.0, 0.0},
		{ 0.1,  0.1, -1.0, 0.0}
	};

	//glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);
	//glLightfv(GL_LIGHT0, GL_POSITION, lightPositions[0]);
	glEnable(GL_LIGHT0);

	// Enable lighting on both front and back sides
	// Note: you must draw the back sides with normals if you want to show it
	//glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);

	// const GLfloat diffuseLight1[] = {1.0, 1.0, 1.0, 1.0};
	//glLightfv(GL_LIGHT1, GL_DIFFUSE, lightPositions[1]);
	//glLightfv(GL_LIGHT1, GL_POSITION, lightPositions[1]);
	//glEnable(GL_LIGHT1);

	//glLightfv(GL_LIGHT2, GL_DIFFUSE   , diffuseLight);
	//glLightfv(GL_LIGHT2, GL_POSITION, lightPositions[2]);
	//glEnable(GL_LIGHT2);

	//glLightfv(GL_LIGHT3, GL_DIFFUSE   , diffuseLight);
	//glLightfv(GL_LIGHT3, GL_POSITION, lightPositions[3]);
	//glEnable(GL_LIGHT3);
}

void GLWidget::initializeGL()
{
	initLighting();
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glLineWidth(settings->lineSize);
	glPointSize(settings->pointSize);
	//glEnable(GL_POINT_SMOOTH);
	setCursor(QCursor(Qt::ArrowCursor));
	//glClearColor(1.0, 1.0, 1.0, 1.0);

	//this->format().setDoubleBuffer(true);
}

void GLWidget::paintGL()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBlendFunc(GL_ONE, GL_ZERO);
	glDisable(GL_ALPHA_TEST);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);

	glEnable(GL_LIGHTING);

	glPushMatrix();
	glMultMatrixf((GLfloat*)settings->transMatrix);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	float residue = 20.0;

	// orthogonal projection
	// Note: the last parameter zFar should not be very large, or the front of the model
	// is very likely to be transparent. Here I use 100.
	glOrtho((-residue - settings->xShift - settings->zShift) / residue * range,	// left
		(residue - settings->xShift + settings->zShift) / residue * range,		// right
		(-residue - settings->yShift - settings->zShift) / residue * range,		// bottom
		(residue - settings->yShift + settings->zShift) / residue * range,		// top
		(-2*residue - settings->zShift - settings->clipping) / (residue) * range,		// zNear
		100);																// zFar

	glMatrixMode(GL_MODELVIEW);
	renderModel();
	glPopMatrix();
	glFlush();
}

void GLWidget::resizeGL(int width, int height)
{
	int max = width > height ? width : height;
	glViewport((width-max)/2, (height-max)/2, max, max);
	//   int side = qMin(width, height);
	//   //glViewport((width - side) / 2, (height - side) / 2, side, side);
	//glViewport(0,0, width, height);
	/*QRect rect = QApplication::desktop()->availableGeometry();*/
	//gluPerspective(90, (GLdouble)(width/height), 5, 5+viewSpeed);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-1.0, +1.0, -1.0, 1.0, 5.0, 60.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslated(0.0, 0.0, -40.0);
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
	mouseLastPos = event->pos();

	const int bufSize = 32; 
	GLuint selectBuf[bufSize] = {0};
	GLint hits, viewport[4];

	if(event->button() == Qt::MidButton)
	{
		mouseMode = TRANSLATE;
		setCursor(QCursor(Qt::SizeAllCursor));
	}
	if (event->button() == Qt::LeftButton)
	{
		mouseMode = ROTATE;
		point2Triple(mouseLastPos.x(),mouseLastPos.y(),lastMousePos3F);
		setCursor( QCursor(Qt::CrossCursor) );
	}
}

void GLWidget::point2Triple(int x, int y, float triplef[3])
{
	float df, a;
	triplef[0] = (2.0 * x - width()) / width();
	triplef[1] = (height() - 2.0 * y) / height();
	df = sqrt(triplef[0] * triplef[0] + triplef[1] * triplef[1]);
	triplef[2] = cos((PI / 2) * ((df < 1.0) ? df : 1.0));
	a = 1.0 / sqrt(triplef[0] * triplef[0] + triplef[1] * triplef[1] + triplef[2] * triplef[2]);
	triplef[0] *= a;
	triplef[1] *= a;
	triplef[2] *= a;
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
	int x= event->x(), y=event->y();

	// --------------------------------------------------
	// Show the 3D coordinates in real time
	//coordTransformation(x,y,lastPos3D);
	//mouse2DPosStr = tr("2D Coordinates: (%1, %2)").arg(x).arg(y);
	//mouse3DPosStr = tr("OpenGL Coordinates: (%1, %2, %3)").arg(lastPos3D[0])
	//	.arg(lastPos3D[1]).arg(lastPos3D[2]);
	//updateGL();	// Must be added in order to flush 
	// ---------------------------------------------------

	switch(mouseMode)
	{
	case ROTATE:
		float currentPosition[3], angle, dx, dy, dz;

		point2Triple(x, y, currentPosition);

		dx = currentPosition[0] - lastMousePos3F[0];
		dy = currentPosition[1] - lastMousePos3F[1];
		dz = currentPosition[2] - lastMousePos3F[2];

		angle = 180.0 * sqrt(dx * dx + dy * dy + dz * dz);

		float rotateAxis[3];
		rotateAxis[0] = lastMousePos3F[1] * currentPosition[2] -
			lastMousePos3F[2] * currentPosition[1];
		rotateAxis[1] = lastMousePos3F[2] * currentPosition[0] -
			lastMousePos3F[0] * currentPosition[2];
		rotateAxis[2] = lastMousePos3F[0] * currentPosition[1] -
			lastMousePos3F[1] * currentPosition[0];
		lastMousePos3F[0] = currentPosition[0];
		lastMousePos3F[1] = currentPosition[1];
		lastMousePos3F[2] = currentPosition[2];

		glPushMatrix();
		glLoadIdentity();
		glRotatef(angle, rotateAxis[0], rotateAxis[1], rotateAxis[2]);
		glMultMatrixf((GLfloat*) settings->transMatrix);
		glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *)settings->transMatrix);
		glPopMatrix();
		updateGL();
		break;

	case TRANSLATE:
		settings->xShift += (x-mouseLastPos.x())/(float)width()*8;
		settings->yShift -= (y-mouseLastPos.y())/(float)height()*8;
		mouseLastPos = event->pos();
		updateGL();
		break;

	default:
		//lastPos = event->pos();
		updateGL();
		break;
	}
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
	mouseMode = NONE;
	setCursor(QCursor(Qt::ArrowCursor));
	updateGL();
}

void GLWidget::resizeEvent(QResizeEvent* event)
{
	QSize size = event->size();
	resizeGL(size.width(),size.height());
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
	int numDegrees = event->delta() / 120;
	settings->zShift += numDegrees * 0.8;
	updateGL();
}


/**************************************************************************/
/*
* Public Slots
*/

void GLWidget::openFile()
{
	// Get the folder path of original animation frames
	QString qstr_filename = QFileDialog::getOpenFileName(this,tr("Open Model"),"",tr("PLY format (*.ply)\nCluster File (*.cluster)"));
	if(qstr_filename == "")
		return;
	char *filename = cstyleString(qstr_filename);
	if( strstr(filename,".ply")!=NULL )
	{
		acvt = new ACVT;
		cout<<"Loading mesh..."<<endl;
		double load_t = get_cpu_time();
		acvt->readMesh(filename);
		load_t = get_cpu_time() - load_t;
		cout<<"\tTime: "<<load_t<<endl;
		acvt->mesh->printInfo();
		cout<<"Done"<<endl;
		flag_isOpened = true;
		range = acvt->mesh->getRange();
	}
	else if (strstr(filename,".cluster")!=NULL)
	{
		// Read a cluster file
		if (flag_isOpened == false || acvt->vertex_number == 0)
		{
			cout<<"You must open a mesh first!"<<endl;
			return;
		}
		cout<<"Reading Cluster File..."<<endl;
		vector<vector<VertexID>*> *cluster_all = new vector<vector<VertexID>*>;
		acvt->open_cluster_file(filename, cluster_all);

		cout<<"Initializing..."<<endl;
		acvt->init_qslim();

		//cout<<"Adding border constraint to the related clusters..."<<endl;
		//acvt->update_border_constraint();

		// Set clusters infos such as quadrics and target coordinates
		cout<<"Setting clusters based on input cluster file..."<<endl;
		acvt->set_clusters_from_input(cluster_all);

		cout<<"Clearing cluster vector..."<<endl;
		for (uint i = 0; i != cluster_all->size(); ++i)
		{
			cluster_all->at(i)->clear();
			ClearVector(*(cluster_all->at(i)));
		}
		cluster_all->clear();
		ClearVector(*cluster_all);

		//acvt->merge_all_clusters_inside_another();

		cout<<"Vertex number in clusters:"<<acvt->get_vertex_number_in_clusters()<<endl;

		cout<<"Assigning cluster colors..."<<endl;
		acvt->assign_cluster_colors();
		show_clusters(true);
		set_connectivity_flag(0);
		cout<<"Done"<<endl;

		programType = QEM;
		updateGL();
	}
	else
	{
		cout<<"Unknown model format. Quiting..."<<endl;
	}
}

void GLWidget::run_QEM(uint vertexNum, double reductionRatio)
{
	if (flag_isOpened)
	{
		// Initialization of QSlim method
		cout<<"Running QSlim..."<<endl;
		double t = get_cpu_time();
		cout<<"Initializing..."<<endl;
		acvt->init_qslim();
		cout<<"\tTime: "<<get_cpu_time() - t<<endl;

		acvt->compute_all_neighbor_clusters();
		set<VertexID> cluster_ids;
		acvt->check_nonmanifold_clusters(cluster_ids);
		cout<<"Non manifold clusters:"<<cluster_ids.size()<<endl;
		cluster_ids.clear();
		ClearSet(cluster_ids);

		// Simplification
		final_vertex_number = vertexNum;
		if (reductionRatio != 0)
			final_vertex_number = (uint)(acvt->vertex_number * reductionRatio);

		cout<<"Simplifying..."<<endl;
		cout<<"Initial Energy:"<<acvt->compute_total_energy()<<endl;
		t = get_cpu_time();
		while (acvt->valid_verts > final_vertex_number)
		{
			acvt->simplify();
		}
		cout<<"\tTime: "<<get_cpu_time() - t<<endl;
		cout<<"Result Energy:"<<acvt->compute_total_energy()<<endl;

		// Show the original QEM connectivity
		cout<<"Updating QEM normals for new mesh..."<<endl;
		acvt->compute_QEM_normals();

		set_connectivity_flag(1);
		cout<<"Assigning cluster colors..."<<endl;
		acvt->assign_cluster_colors();

		cout<<"Vertex Number before merging:"<<acvt->get_current_vertex_number()<<endl;
		cout<<acvt->valid_verts<<endl;
		cout<<"Merging clusters with only one neighbor..."<<endl;
		acvt->merge_all_clusters_inside_another();
		cout<<"Vertex Number after merging:"<<acvt->get_current_vertex_number()<<endl;
		cout<<acvt->valid_verts<<endl;

		//cout<<"Updating normals and center for new mesh..."<<endl;
		//acvt->update_normals();

		//cout<<"Adding border constraint to the related clusters..."<<endl;
		//acvt->update_border_constraint();

		show_clusters_flag = true;

		programType = QEM;
		updateGL();
		cout<<"Done"<<endl;
	}
	else
	{
		cout<<"Please open the model first!"<<endl;
	}
}

void GLWidget::saveFile()
{
	// Save Result Model
	cout<<"Saving..."<<endl;
	QString qstr_filename = QFileDialog::getSaveFileName(this,tr("Open Model"),"",tr("PLY format (*.ply)\nCluster File (*.cluster)"));
	if(qstr_filename == "")
		return;

	char *filename = cstyleString(qstr_filename);
	if( strstr(filename,".ply")!=NULL )
	{
		acvt->save_result_model(connectivity_flag);
		acvt->mesh->writeNewPLYWithExtraData(filename);
	}
	else if (strstr(filename,".cluster")!=NULL)
	{
		acvt->save_cluster_in_file(filename);
	}
	else
	{
		cout<<"Format Not Supported"<<endl;
	}
	cout<<"Done"<<endl;
}

void GLWidget::run_QEM_Optimization()
{
	if (programType == QEM || programType == QEMOPTIMIZATION)
	{
		// Firstly, show clusters
		if (!show_clusters_flag)
		{
			show_clusters(true);
			updateGL();
		}
		//acvt->compute_cluster_belong();
		acvt->compute_all_neighbor_clusters();
		cout<<"Running optimization Process..."<<endl;
		uint count = 0;
		double t = get_cpu_time();
		double energy_initial = 0;
		do
		{
			cout<<"\tTime:"<<get_cpu_time()-t<<endl;
			cout<<"Iteration "<<count<<": "<<endl;
			//cout<<"Energy = "<<acvt->compute_total_energy()<<endl;
			count++;
			if (count > 300)
			{
				break;
				//double energy_now = acvt->compute_total_energy();
				//if (fabs(energy_now - energy_initial) < 1e-2)
				//{
				//	break;
				//}
				//energy_initial = energy_now;
			}
			//cout<<"\tVertex Number:"<<acvt->get_vertex_number_in_clusters()<<endl;
			updateGL();
			t = get_cpu_time();
		}while (acvt->swap_once());
		cout<<"Time for swapping:"<<get_cpu_time()-t<<endl;

		cout<<"Vertex Number after swapping:"<<acvt->get_current_vertex_number()<<endl;

		// Merge some small clusters into their belonging "father" clusters
		cout<<"Merging clusters with only one neighbor..."<<endl;
		acvt->merge_all_clusters_inside_another();
		cout<<"Vertex Number after merging: "<<acvt->get_current_vertex_number()<<endl;

		// Show the original QEM connectivity
		cout<<"Updating QEM normals for new mesh..."<<endl;
		acvt->compute_QEM_normals();
		set_connectivity_flag(1);

		cout<<"Ending optimization..."<<endl;
		programType = QEMOPTIMIZATION;
		cout<<"Done"<<endl;
	}
	else
	{
		cout<<"Please run QEM process first!"<<endl;
	}
}

void GLWidget::run_QEM_Optimization_SwapOnce()
{
	if (programType == QEM || programType == QEMOPTIMIZATION)
	{
		cout<<"Swapping Once..."<<endl;
		if (!show_clusters_flag)
		{
			show_clusters(true);
			updateGL();
		}
		cout<<"Initial Energy = "<<acvt->compute_total_energy()<<endl;
		acvt->swap_once();
		cout<<"Result Energy = "<<acvt->compute_total_energy()<<endl;
		cout<<"Done"<<endl;
		programType = QEMOPTIMIZATION;
		updateGL();
	}
}

void GLWidget::run_split_clusters()
{
	cout<<"Splitting Clusters into separate connected components..."<<endl;

	if (!show_clusters_flag)
	{
		show_clusters(true);
		updateGL();
	}
	double t = get_cpu_time();
	acvt->split_clusters();
	cout<<"Time for splitting clusters:"<<get_cpu_time()-t<<endl;
	cout<<"Merging clusters with only one neighbor..."<<endl;
	acvt->merge_all_clusters_inside_another();
	cout<<"Assign cluster colors..."<<endl;
	acvt->assign_cluster_colors();

	cout<<"Done"<<endl;
	programType = SPLITCLUSTERS;
	updateGL();
}

void GLWidget::run_create_connectivity()
{
	if (programType != BLANK)
	{
		set_connectivity_flag(2);
		cout<<"Compute new connectivity..."<<endl;
		double t = get_cpu_time();
		acvt->compute_new_faces_nonrecursive();
		//acvt->compute_new_faces();
		cout<<"Time="<<get_cpu_time()-t<<endl;
		

		programType = NEWCONNECTIVITY;
		cout<<"Done"<<endl;
	}
	else
	{
		cout<<"You should run QEM or QEM Optimization or Split Clusters first!"<<endl;
	}
}

void GLWidget::run_ReQEM_after_split(uint vertexNum, double reductionRatio)
{
	if (programType == SPLITCLUSTERS)
	{
		cout<<"Re-running QEM after splitting clusters..."<<endl;
		
		cout<<"Initializing..."<<endl;
		acvt->compute_new_edgelist_heap();

		cout<<"Run QEM again..."<<endl;
		cout<<"Initial Energy:"<<acvt->compute_total_energy()<<endl;
		double t = get_cpu_time();
		final_vertex_number = vertexNum;
		if (reductionRatio != 0)
			final_vertex_number = (uint)(acvt->vertex_number * reductionRatio);
		while (acvt->valid_verts > final_vertex_number)
		{
			acvt->simplify();
		}		
		cout<<"\tTime: "<<get_cpu_time() - t<<endl;
		cout<<"Result Energy:"<<acvt->compute_total_energy()<<endl;

		show_clusters_flag = true;
		programType = REQEM;
		updateGL();
		cout<<"Done"<<endl;
	}
	else
	{
		cout<<"You should run 'Spliting Clusters' Step first!"<<endl;
	}
}

void GLWidget::run_split_nonmanifold_clusters()
{
	cout<<"Splitting non-manifold clusters..."<<endl;
	uint vtxNum = acvt->get_current_vertex_number();
	cout<<"Vertex Number now:"<<vtxNum<<endl;
	uint count = 0;
	do 
	{
		acvt->compute_all_neighbor_clusters();
		if (count > 300)
		{
			break;
		}
		cout<<"Iteration "<<++count<<":"<<endl;
	} while (acvt->split_nonmanifold_clusters());

	//cout<<"Vertex Number in clusters:"<<acvt->get_vertex_number_in_clusters();

	//acvt->compute_all_neighbor_clusters();
	//acvt->split_nonmanifold_clusters();

	uint vtxNumNow = acvt->get_current_vertex_number();
	cout<<"Vertex Number in the Result:"<<vtxNumNow<<endl;
	cout<<"Number of vertices added:"<<vtxNumNow-vtxNum<<endl;
	//cout<<"\tTime: "<<get_cpu_time() - t<<endl;

	acvt->assign_cluster_colors();
	show_clusters_flag = true;
	programType = MANIFOLDTOPOLOGY;
	updateGL();
	cout<<"Done"<<endl;
}

void GLWidget::run_compute_vertex_weight()
{
	cout<<"Computing vertex weights..."<<endl;
	acvt->compute_vertex_weight();
	show_vertex_weights(true);
	updateGL();
	cout<<"Done"<<endl;
}

void GLWidget::run_all_steps(uint vertexNum, double reductionRatio)
{
	if (flag_isOpened)
	{
		/*
		* Step 0: Initialization of QEM
		*/
		cout<<"Step 0: Initialization of QEM step..."<<endl;
		//double t = get_cpu_time();
		acvt->init_qslim();
		//cout<<"\tTime 0: "<<get_cpu_time() - t<<endl;

		/*
		* Step 1: QEM step
		*/
		cout<<"Step 1: Running QEM Simplification..."<<endl;
		final_vertex_number = vertexNum;
		if (reductionRatio != 0)
			final_vertex_number = (uint)(acvt->vertex_number * reductionRatio);
		cout<<"  Initial Energy:"<<acvt->compute_total_energy()<<endl;
		while (acvt->valid_verts > final_vertex_number)
		{
			acvt->simplify();
		}
		cout<<"  Result Energy:"<<acvt->compute_total_energy()<<endl;

		cout<<"Merging clusters with only one neighbor..."<<endl;
		acvt->merge_all_clusters_inside_another();
		cout<<"Vertex Number after merging: "<<acvt->get_current_vertex_number()<<endl;

		cout<<"Updating QEM normals for new mesh..."<<endl;
		acvt->compute_QEM_normals();

		/*
		* Step 2: Optimization step
		*/
		cout<<"Step 2: Running optimization after QEM..."<<endl;
		acvt->compute_all_neighbor_clusters();
		uint count = 0;
		do
		{
			cout<<"Iteration "<<count<<": "<<endl;
			count++;
			if (count > 300)
			{
				break;
			}
		}while (acvt->swap_once());
		// Merge some small clusters into their belonging "father" clusters
		cout<<"Merging clusters with only one neighbor..."<<endl;
		acvt->merge_all_clusters_inside_another();
		cout<<"Vertex Number after merging: "<<acvt->get_current_vertex_number()<<endl;
		// Show the original QEM connectivity
		cout<<"Updating QEM normals for new mesh..."<<endl;
		acvt->compute_QEM_normals();
		set_connectivity_flag(1);

		/*
		* Step 3: Split clusters into separate connected components
		*/
		cout<<"Step 3: Splitting Clusters into separate connected components..."<<endl;
		acvt->split_clusters();
		cout<<"Merging clusters with only one neighbor..."<<endl;
		acvt->merge_all_clusters_inside_another();
		
		/*
		* Step 4: Re-QEM after split
		*/
		cout<<"Step 4: Re-running QEM after splitting clusters..."<<endl;
		cout<<"Initializing..."<<endl;
		acvt->compute_new_edgelist_heap();
		cout<<"Run QEM again..."<<endl;
		cout<<"Initial Energy:"<<acvt->compute_total_energy()<<endl;
		//final_vertex_number = vertexNum;
		//if (reductionRatio != 0)
		//	final_vertex_number = (uint)(acvt->vertex_number * reductionRatio);
		while (acvt->valid_verts > final_vertex_number)
		{
			acvt->simplify();
		}
		cout<<"Result Energy:"<<acvt->compute_total_energy()<<endl;

		/*
		* Step 5: Split non-manifold clusters
		*/
		cout<<"Step 5: Splitting non-manifold clusters..."<<endl;
		uint vtxNum = acvt->get_current_vertex_number();
		cout<<"Vertex Number now:"<<vtxNum<<endl;
		count = 0;
		do 
		{
			acvt->compute_all_neighbor_clusters();
			cout<<"Iteration "<<++count<<":"<<endl;
			if (count > 300)
			{
				break;
			}
		} while (acvt->split_nonmanifold_clusters());
		uint vtxNumNow = acvt->get_current_vertex_number();
		cout<<"Vertex Number in the Result:"<<vtxNumNow<<endl;
		cout<<"Number of vertices added:"<<vtxNumNow-vtxNum<<endl;

		/*
		* Step 6: Compute New Connectivity
		*/
		cout<<"Step 6: Compute new connectivity..."<<endl;
		//acvt->compute_new_faces();
		acvt->compute_new_faces_nonrecursive();
		/*
		* Final set up
		*/
		cout<<"Assigning cluster colors..."<<endl;
		acvt->assign_cluster_colors();
		set_connectivity_flag(2);
		show_clusters_flag = true;
		programType = NEWCONNECTIVITY;
		
		updateGL();
		cout<<"All Done!"<<endl;
	}
	else
	{
		cout<<"Please open the model first!"<<endl;
	}
}

void GLWidget::run_compute_vertex_quality()
{
	if (flag_isOpened && programType != BLANK)
	{
		cout<<"Computing vertex quality..."<<endl;
		acvt->compute_final_vertex_quality();
		cout<<"Done"<<endl;
	}
	else
	{
		cout<<"Please open the model first!"<<endl;
	}
}


/************************************************************************/
/* 
* Rendering functions
*/

void GLWidget::renderModel()
{
	if (flag_isOpened)
	{
		if (programType == BLANK)
		{
			switch (settings->meshDisplayMode)
			{
			case VERTICES:
				if (show_vertex_weight_flag)
				{
					displayVertexWeightDistributions();
				}
				else
				{
					displayPoints();
				}
				break;
			case EDGES:
				displayEdges();
				break;	
			case HIDDENLINES:
				displayHiddenLinesRemoval();
				break;

			case FLAT:
				displayTriangles(GL_FLAT);
				break;

			case SMOOTH:
				displayTriangles(GL_SMOOTH);
				break;

			case TEXTURES:
				displayTextures();
				break;

			case NOTHING:
				break;

			default:
				break;
			}
		}
		else
		{
			if (show_clusters_flag)
			{
				switch (settings->meshDisplayMode)
				{
					case VERTICES:
						displayClusters();
						break;
					case FLAT:
						displayResultClustersInTriangles(GL_FLAT);
						displayHalfEdges();
						break;
					case SMOOTH:
						displayResultClustersInTriangles(GL_SMOOTH);
						break;
					default:
						break;
				}
			}
			else
			{
				switch (settings->meshDisplayMode)
				{
				case VERTICES:
					displayResultPoints();
					break;
				case EDGES:
					if (connectivity_flag == 0)
						displayEdges();
					else if (connectivity_flag == 1)
						displayResultEdges_QEM();
					else 
						displayResultEdges_Clusters();
					break;
				case FLAT:
					displayResultTriangles_Clusters_ColorInverseNormalFace(GL_FLAT);
					displayResultEdges_Clusters();
					break;
				case SMOOTH:
					if (connectivity_flag == 0)
						displayEdges();
					else if (connectivity_flag == 1)
						displayResultEdges_QEM();
					else 
						displayHalfEdges();
					break;
				default:
					break;
				}
			}
		}
	}
}

void GLWidget::displayPoints()
{
	uint vertex_number = acvt->vertex_number;
	float *ctr = acvt->mesh->getCenter();

	if (show_crest_line_flag)
	{
#ifdef CRESTLINE
		glPushMatrix();
		glPolygonMode (GL_FRONT, GL_LINE);
		glMaterialfv(GL_FRONT, GL_AMBIENT, settings->feature_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, settings->feature_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, settings->feature_shinness);
		glLineWidth(settings->feature_line_size);
		glBegin(GL_LINES);
		for(int i = 0 ; i != acvt->ridge->edgeNumber; i++)
		{
			int id1 = acvt->ridge->edges->at(i).first.first;
			int id2 = acvt->ridge->edges->at(i).first.second;
			if (acvt->ridge->is_feature(id1) && acvt->ridge->is_feature(id2))
			{
				float *v1 = acvt->ridge->getVertex(id1);
				float *v2 = acvt->ridge->getVertex(id2);

				//glNormal3fv(acvt->ridge->getNormal(id1));
				glVertex3f(v1[0]-ctr[0],v1[1]-ctr[1],v1[2]-ctr[2]);
				//glNormal3fv(acvt->ridge->getNormal(id2));
				glVertex3f(v2[0]-ctr[0],v2[1]-ctr[1],v2[2]-ctr[2]);
			}
		}
		for(int i = 0 ; i != acvt->ravine->edgeNumber; i++)
		{
			int id1 = acvt->ravine->edges->at(i).first.first;
			int id2 = acvt->ravine->edges->at(i).first.second;
			if (acvt->ravine->is_feature(id1) && acvt->ravine->is_feature(id2))
			{
				float *v1 = acvt->ravine->getVertex(id1);
				float *v2 = acvt->ravine->getVertex(id2);

				//glNormal3fv(acvt->ravine->getNormal(id1));
				glVertex3f(v1[0]-ctr[0],v1[1]-ctr[1],v1[2]-ctr[2]);
				//glNormal3fv(acvt->ravine->getNormal(id2));
				glVertex3f(v2[0]-ctr[0],v2[1]-ctr[1],v2[2]-ctr[2]);
			}
		}
		glEnd();
		glPopMatrix();
#endif
	}

	// opengl setting
	glPushMatrix();
	glPolygonMode (GL_FRONT_AND_BACK, GL_POINT);
	glPointSize(settings->pointSize);
	glMaterialfv(GL_FRONT, GL_AMBIENT, settings->point_ambient);
	glMaterialfv(GL_FRONT,GL_DIFFUSE, settings->point_diffuse);
	glMaterialf(GL_FRONT, GL_SHININESS, settings->point_shinness);
	glBegin(GL_POINTS);
	for(uint i = 0 ; i != vertex_number ; i++)
	{		
		float *v = acvt->mesh->getVertex(i);
		//vet = getVertexAt(i);
		glNormal3fv(acvt->getOriginalNormal(i));
		glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
	}
	glEnd();
	glPopMatrix();
}

void GLWidget::displayResultPoints()
{
	float *ctr = acvt->mesh->getCenter();

	glPushMatrix();
	glPolygonMode (GL_FRONT_AND_BACK, GL_POINT);
	glPointSize(settings->pointSize);
	glMaterialfv(GL_FRONT, GL_AMBIENT, settings->point_ambient);
	glMaterialfv(GL_FRONT,GL_DIFFUSE, settings->point_diffuse);
	glMaterialf(GL_FRONT, GL_SHININESS, settings->point_shinness);
	glBegin(GL_POINTS);
	for(uint i = 0 ; i != acvt->vertex_number; i++)
	{	
		VertexQuadrics *vtxQuad = acvt->vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			VectorType v = acvt->TargetVertex(i);
			glNormal3fv(acvt->getCurrentNormal(i));
			glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
		}
	}
	glEnd();
	glPopMatrix();
}



void GLWidget::displayClusters()
{
	float *ctr = acvt->mesh->getCenter();

	glPolygonMode (GL_FRONT_AND_BACK, GL_POINT);
	glPointSize(settings->pointSize);
	glMaterialf(GL_FRONT, GL_SHININESS, 20);
	glBegin(GL_POINTS);
	for(uint i = 0 ; i != acvt->vertex_number; i++)
	{
		VertexQuadrics *vtxQuad = acvt->vtxQuadricsList->at(i);
		if (!vtxQuad->cluster->empty())
		{
			for (set<VertexID>::iterator iter = vtxQuad->cluster->begin();
				iter != vtxQuad->cluster->end(); ++iter)
			{
				float *v = acvt->mesh->getVertex(*iter);
				//glMaterialfv(GL_FRONT, GL_AMBIENT, colors[vtxQuad->get_cluster_color_no()]);
				//glMaterialfv(GL_FRONT,GL_DIFFUSE, colors[vtxQuad->get_cluster_color_no()]);
				//glMaterialfv(GL_FRONT, GL_AMBIENT, vtxQuad->clusterColor);
				glMaterialfv(GL_FRONT, GL_AMBIENT, vtxQuad->clusterColor);
				glMaterialfv(GL_FRONT,GL_DIFFUSE, vtxQuad->clusterColor);
				glNormal3fv(acvt->getOriginalNormal(*iter));
				glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
			}	
		}	
	}
	glEnd();	

	//glPushMatrix();
	//glPolygonMode (GL_FRONT_AND_BACK, GL_POINT);
	//glPointSize(settings->face_shinnes);
	//glMaterialfv(GL_FRONT, GL_AMBIENT, settings->feature_ambient);
	//glMaterialfv(GL_FRONT,GL_DIFFUSE, settings->feature_diffuse);
	//glMaterialf(GL_FRONT, GL_SHININESS, settings->feature_shinness);
	//glBegin(GL_POINTS);
	//for(uint i = 0 ; i != acvt->vertex_number; i++)
	//{	
	//	VertexQuadrics *vtxQuad = acvt->vtxQuadricsList->at(i);
	//	if (!vtxQuad->cluster->empty())
	//	{
	//		VectorType v = acvt->TargetVertex(i);
	//		glNormal3fv(acvt->getCurrentNormal(i));
	//		glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
	//	}
	//}
	//glEnd();
	//glPopMatrix();
}

void GLWidget::displayResultEdges_QEM()
{
	//uint face_number = acvt->face_number;
	//VectorType ctr = acvt->ctr;
	float *ctr = acvt->mesh->getCenter();

	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	glMaterialfv(GL_FRONT, GL_AMBIENT, settings->line_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, settings->line_diffuse);
	glMaterialf(GL_FRONT, GL_SHININESS, settings->line_shinness);
	glLineWidth(settings->lineSize);
	glBegin(GL_LINES);
	for(int i = 0 ; i != acvt->faceList->size() ; i++)
	{
		FaceInfo *info_f = acvt->faceList->at(i);
		if (info_f->is_exist())
		{
			uint *t = info_f->v;
			for(int j = 0 ; j < 3 ; j++)
			{
				VectorType v1 = acvt->TargetVertex(t[j]);
				VectorType v2 = acvt->TargetVertex(t[(j+1)%3]);

				glNormal3fv(acvt->getCurrentNormal(t[j]));
				glVertex3f(v1[0]-ctr[0],v1[1]-ctr[1],v1[2]-ctr[2]);
				glNormal3fv(acvt->getCurrentNormal(t[(j+1)%3]));
				glVertex3f(v2[0]-ctr[0],v2[1]-ctr[1],v2[2]-ctr[2]);
			}
		}
	}

	glEnd();
}

void GLWidget::displayResultEdges_Clusters()
{
	//uint face_number = acvt->face_number;
	//VectorType ctr = acvt->ctr;
	float *ctr = acvt->mesh->getCenter();
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	glMaterialfv(GL_FRONT, GL_AMBIENT, settings->line_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, settings->line_diffuse);
	glMaterialf(GL_FRONT, GL_SHININESS, settings->line_shinness);
	glLineWidth(settings->lineSize);
	glBegin(GL_LINES);
	for (map<pair<uint,uint>, uint>::iterator iter = acvt->edgePairs_new.begin(); 
		iter != acvt->edgePairs_new.end(); ++iter)
	{
		// Create this edge
		VertexID a = (*iter).first.first, b = (*iter).first.second;
		VectorType v1 = acvt->TargetVertex(a);
		VectorType v2 = acvt->TargetVertex(b);

		glNormal3fv(acvt->getCurrentNormal(a));
		glVertex3f(v1[0]-ctr[0],v1[1]-ctr[1],v1[2]-ctr[2]);
		glNormal3fv(acvt->getCurrentNormal(b));
		glVertex3f(v2[0]-ctr[0],v2[1]-ctr[1],v2[2]-ctr[2]);
	}
	glEnd();
}

void GLWidget::displayEdges()
{
	uint face_number = acvt->face_number;
	float *ctr = acvt->mesh->getCenter();

	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	glMaterialfv(GL_FRONT, GL_AMBIENT, settings->line_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, settings->line_diffuse);
	glMaterialf(GL_FRONT, GL_SHININESS, settings->line_shinness);
	glLineWidth(settings->lineSize);
	glBegin(GL_LINES);
	for(int i = 0 ; i != face_number ; i++)
	{
		int *t = acvt->mesh->getTriangle(i);
		for(int j = 0 ; j < 3 ; j++)
		{
			float *v1 = acvt->mesh->getVertex(t[j]);
			float *v2 = acvt->mesh->getVertex(t[(j+1)%3]);

			glNormal3fv(acvt->getOriginalNormal(t[j]));
			glVertex3f(v1[0]-ctr[0],v1[1]-ctr[1],v1[2]-ctr[2]);
			glNormal3fv(acvt->getOriginalNormal(t[(j+1)%3]));
			glVertex3f(v2[0]-ctr[0],v2[1]-ctr[1],v2[2]-ctr[2]);
		}
	}
	glEnd();
//#ifdef BORDERS
//	glMaterialfv(GL_FRONT, GL_AMBIENT, settings->feature_ambient);
//	glMaterialfv(GL_FRONT, GL_DIFFUSE, settings->feature_diffuse);
//	glMaterialf(GL_FRONT, GL_SHININESS, settings->feature_shinness);
//	glLineWidth(settings->lineSize);
//	glBegin(GL_LINES);
//	for (map<pair<uint,uint>, uint>::iterator iter = acvt->edgePairs.begin(); 
//		iter != acvt->edgePairs.end(); ++iter)
//	{
//		if ((*iter).second == 1)
//		{
//			VertexID a = (*iter).first.first, b = (*iter).first.second;
//			float *v1 = acvt->mesh->getVertex(a);
//			float *v2 = acvt->mesh->getVertex(b);
//
//			glNormal3fv(acvt->getOriginalNormal(a));
//			glVertex3f(v1[0]-ctr[0],v1[1]-ctr[1],v1[2]-ctr[2]);
//			glNormal3fv(acvt->getOriginalNormal(b));
//			glVertex3f(v2[0]-ctr[0],v2[1]-ctr[1],v2[2]-ctr[2]);
//		}
//	}
//	glEnd();
//#endif

	if (show_crest_line_flag)
	{
#ifdef CRESTLINE
		glPushMatrix();
		glPolygonMode (GL_FRONT, GL_LINE);
		glMaterialfv(GL_FRONT, GL_AMBIENT, settings->feature_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, settings->feature_diffuse);
		glMaterialf(GL_FRONT, GL_SHININESS, settings->feature_shinness);
		glLineWidth(settings->feature_line_size);
		glBegin(GL_LINES);
		for(int i = 0 ; i != acvt->ridge->edgeNumber; i++)
		{
			int id1 = acvt->ridge->edges->at(i).first.first;
			int id2 = acvt->ridge->edges->at(i).first.second;
			if (acvt->ridge->is_feature(id1) && acvt->ridge->is_feature(id2))
			{
				float *v1 = acvt->ridge->getVertex(id1);
				float *v2 = acvt->ridge->getVertex(id2);

				//glNormal3fv(acvt->ridge->getNormal(id1));
				glVertex3f(v1[0]-ctr[0],v1[1]-ctr[1],v1[2]-ctr[2]);
				//glNormal3fv(acvt->ridge->getNormal(id2));
				glVertex3f(v2[0]-ctr[0],v2[1]-ctr[1],v2[2]-ctr[2]);
			}
		}
		for(int i = 0 ; i != acvt->ravine->edgeNumber; i++)
		{
			int id1 = acvt->ravine->edges->at(i).first.first;
			int id2 = acvt->ravine->edges->at(i).first.second;
			if (acvt->ravine->is_feature(id1) && acvt->ravine->is_feature(id2))
			{
				float *v1 = acvt->ravine->getVertex(id1);
				float *v2 = acvt->ravine->getVertex(id2);

				//glNormal3fv(acvt->ravine->getNormal(id1));
				glVertex3f(v1[0]-ctr[0],v1[1]-ctr[1],v1[2]-ctr[2]);
				//glNormal3fv(acvt->ravine->getNormal(id2));
				glVertex3f(v2[0]-ctr[0],v2[1]-ctr[1],v2[2]-ctr[2]);
			}
		}
		glEnd();
		glPopMatrix();
#endif
	}
}

void GLWidget::displayHiddenLinesRemoval()
{
	uint face_number = acvt->face_number;
	float *ctr = acvt->mesh->getCenter();

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glFrontFace(GL_CCW);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

	glMaterialfv(GL_FRONT, GL_AMBIENT, settings->line_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, settings->line_diffuse);
	glMaterialf(GL_FRONT, GL_SHININESS, settings->line_shinness);
	glBegin(GL_TRIANGLES);
	for(int i = 0 ; i != face_number ; i++)
	{
		int *t = acvt->mesh->getTriangle(i);
		for(int j = 0 ; j < 3 ; j++)
		{
			float *v = acvt->mesh->getVertex(t[j]);
			glNormal3fv(acvt->getOriginalNormal(t[j]));
			glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
		}
	}
	glEnd();
	glDisable(GL_CULL_FACE);
}

void GLWidget::displayTriangles(GLenum disType)
{
	uint face_number = acvt->face_number;
	float *ctr = acvt->mesh->getCenter();

	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, settings->face_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, settings->face_diffuse);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, settings->face_shinnes);
	glLineWidth(settings->lineSize);
	glShadeModel(disType);
	glBegin(GL_TRIANGLES);
	for(int i = 0 ; i != face_number ; i++)
	{
		int *t = acvt->mesh->getTriangle(i);
		for(int j = 0 ; j < 3 ; j++)
		{
			float *v = acvt->mesh->getVertex(t[j]);
			glNormal3fv(acvt->getOriginalNormal(t[j]));
			glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
		}
	}
	glEnd();
}

void GLWidget::displayResultTriangles_QEM(GLenum disType)
{
	float *ctr = acvt->mesh->getCenter();

	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, settings->face_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, settings->face_diffuse);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, settings->face_shinnes);
	glLineWidth(settings->lineSize);
	glShadeModel(disType);
	glBegin(GL_TRIANGLES);
	for(uint i = 0 ; i != acvt->faceList->size() ; i++)
	{
		FaceInfo *info_f = acvt->faceList->at(i);
		if (info_f->is_exist())
		{
			for (uint j = 0; j != 3; ++j)
			{
				VectorType v = acvt->TargetVertex(info_f->v[j]);
				glNormal3fv(acvt->getCurrentNormal(info_f->v[j]));
				glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
			}
		}
	}
	glEnd();
}

void GLWidget::displayResultTriangles_Clusters(GLenum disType)
{
	float *ctr = acvt->mesh->getCenter();
	float colors[4];
	colors[0] = colors[1] = colors[2] = 0.7;
	colors[3] = 1.0;

	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, settings->face_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colors);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, settings->face_shinnes);
	glLineWidth(settings->lineSize);
	glShadeModel(disType);
	glBegin(GL_TRIANGLES);
	for(uint i = 0 ; i != acvt->faceList_new->size() ; i++)
	{
		FaceInfo *info_f = acvt->faceList_new->at(i);
		if (info_f->is_exist())
		{
			for (uint j = 0; j != 3; ++j)
			{
				VectorType v = acvt->TargetVertex(info_f->v[j]);
				glNormal3fv(acvt->getCurrentNormal(info_f->v[j]));
				glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
			}
		}
	}
	glEnd();
}

void GLWidget::displayResultTriangles_Clusters_ColorInverseNormalFace(GLenum disType)
{
	float *ctr = acvt->mesh->getCenter();
	float colors[4];
	colors[0] = colors[1] = colors[2] = 0.7;
	colors[3] = 1.0;
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glLineWidth(settings->lineSize);
	glShadeModel(disType);
	glBegin(GL_TRIANGLES);
	for(uint i = 0 ; i != acvt->faceList_new->size() ; i++)
	{
		FaceInfo *info_f = acvt->faceList_new->at(i);
		if (info_f->is_exist())
		{
			for (uint j = 0; j != 3; ++j)
			{
				VectorType v = acvt->TargetVertex(info_f->v[j]);
				glMaterialfv(GL_FRONT, GL_AMBIENT, settings->face_ambient);
				glMaterialfv(GL_FRONT, GL_DIFFUSE, colors);
				glMaterialf(GL_FRONT, GL_SHININESS, settings->face_shinnes);
				glNormal3fv(acvt->getCurrentNormal(info_f->v[j]));
				glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
			}
			/*if (i == 0)
			{
			for (uint j = 0; j != 3; ++j)
			{
			VectorType v = acvt->TargetVertex(info_f->v[j]);
			glMaterialfv(GL_FRONT, GL_AMBIENT, settings->face_ambient);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, colors);
			glMaterialf(GL_FRONT, GL_SHININESS, settings->face_shinnes);
			glNormal3fv(acvt->getCurrentNormal(info_f->v[j]));
			glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
			}
			}
			else
			{
			FaceInfo *info_f_last = acvt->faceList_new->at(i-1);
			if (info_f_last->no != info_f->no)
			{
			for (uint j = 0; j != 3; ++j)
			{
			VectorType v = acvt->TargetVertex(info_f->v[j]);
			glMaterialfv(GL_FRONT, GL_AMBIENT, settings->feature_ambient);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, settings->feature_diffuse);
			glMaterialf(GL_FRONT, GL_SHININESS, settings->feature_shinness);
			glNormal3fv(acvt->getCurrentNormal(info_f->v[j]));
			glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
			}
			}
			else
			{
			for (uint j = 0; j != 3; ++j)
			{
			VectorType v = acvt->TargetVertex(info_f->v[j]);
			glMaterialfv(GL_FRONT, GL_AMBIENT, settings->face_ambient);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, colors);
			glMaterialf(GL_FRONT, GL_SHININESS, settings->face_shinnes);
			glNormal3fv(acvt->getCurrentNormal(info_f->v[j]));
			glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
			}
			}
			}*/
		}
	}
	glEnd();
}

void GLWidget::displayResultClustersInTriangles(GLenum disType)
{
	float *ctr = acvt->mesh->getCenter();

	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glLineWidth(settings->lineSize);
	glShadeModel(disType);
	glBegin(GL_TRIANGLES);
	for(uint i = 0 ; i != acvt->face_number ; i++)
	{
		int *t = acvt->mesh->getTriangle(i);
		for (uint j = 0; j != 3; ++j)
		{
			float* v = acvt->mesh->getVertex(t[j]);
			VertexID cid = acvt->getClusterBelong(t[j]);
			VertexQuadrics *vtxQuad = acvt->vtxQuadricsList->at(cid);
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, vtxQuad->clusterColor);
			glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE, vtxQuad->clusterColor);
			glNormal3fv(acvt->getOriginalNormal(t[j]));
			glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);

			/*	float *n = acvt->getOriginalNormal(t[j]);
			float nn[3];
			nn[0] = -n[0]; nn[1] = -n[1]; nn[2] = -n[2]; 
			glMaterialfv(GL_FRONT, GL_AMBIENT, vtxQuad->clusterColor);
			glMaterialfv(GL_FRONT,GL_DIFFUSE, vtxQuad->clusterColor);
			glNormal3fv(nn);
			glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);*/
		}
	}
	glEnd();
}

void GLWidget::displayHalfEdges()
{
	float *ctr = acvt->mesh->getCenter();

	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	float dif[4];
	dif[0] = dif[1] = dif[2] = 0.0; dif[3] = 1.0; 
	glMaterialfv(GL_FRONT, GL_AMBIENT, settings->line_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	glMaterialf(GL_FRONT, GL_SHININESS, settings->line_shinness);
	glLineWidth(5.0);
	glBegin(GL_LINES);
	for (uint i = 0; i != acvt->halfEdges->size(); ++i)
	{
		vector<HalfEdge*>* vecTemp = acvt->halfEdges->at(i);
		if (!vecTemp->empty())
		{
			for (uint j = 0; j != vecTemp->size(); ++j)
			{
				// Create this edge
				HalfEdge *halfedgenow = vecTemp->at(j);
				if (!halfedgenow->inside_face_flag)
				{
					VertexID a = halfedgenow->v1, b = halfedgenow->v2;
					VectorType v1 = acvt->TargetVertex(a);
					VectorType v2 = acvt->TargetVertex(b);

					glNormal3fv(acvt->getCurrentNormal(a));
					glVertex3f(v1[0]-ctr[0],v1[1]-ctr[1],v1[2]-ctr[2]);
					glNormal3fv(acvt->getCurrentNormal(b));
					glVertex3f(v2[0]-ctr[0],v2[1]-ctr[1],v2[2]-ctr[2]);
				}
			}
		}
	}
	glEnd();
}

void GLWidget::displayTextures()
{
	uint face_number = acvt->face_number;
	float *ctr = acvt->mesh->getCenter();
	float maxVal = 0, minVal = 0;
	float materialClr[4];
	if (acvt->mesh->extraTypeNumber > 0)
	{
		acvt->mesh->getVtxMaxMinExtraData(0, maxVal, minVal);
	}
	

	//glHint (GL_POLYGON_SMOOTH, GL_NICEST);
	//glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	//glLineWidth(settings->lineSize);
	//glShadeModel(GL_SMOOTH);
	//glBegin(GL_TRIANGLES);
	//for(int i = 0 ; i != face_number ; i++)
	//{
	//	int *t = acvt->mesh->getTriangle(i);
	//	for(int j = 0 ; j < 3 ; j++)
	//	{
	//		float *v = acvt->mesh->getVertex(t[j]);
	//		float extraVal = acvt->mesh->getVtxExtraData(t[j],0);
	//		computeRGBcolorOfWeight(minVal, maxVal, extraVal, materialClr);
	//		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialClr);
	//		glNormal3fv(acvt->getOriginalNormal(t[j]));
	//		glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
	//	}
	//}
	//glEnd();

	glPolygonMode (GL_FRONT_AND_BACK, GL_POINT);
	glPointSize(settings->pointSize);
	glMaterialf(GL_FRONT, GL_SHININESS, 20);
	glBegin(GL_POINTS);
	for(uint i = 0 ; i != acvt->vertex_number; i++)
	{
		float *v = acvt->mesh->getVertex(i);
		float extraVal = 0;
		if (acvt->mesh->extraTypeNumber > 0)
		{
			extraVal = acvt->mesh->getVtxExtraData(i,0);
		}
		computeRGBcolorOfWeight(minVal, maxVal, extraVal, materialClr);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialClr);
		glNormal3fv(acvt->getOriginalNormal(i));
		glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
	}
	glEnd();	
}

void GLWidget::displayVertexWeightDistributions()
{
	uint face_number = acvt->face_number;
	float *ctr = acvt->mesh->getCenter();
	float materialClr[4];

	float maxVal = *max_element(acvt->vtxQuadricsWeights->begin(), acvt->vtxQuadricsWeights->end());
	float minVal = *min_element(acvt->vtxQuadricsWeights->begin(), acvt->vtxQuadricsWeights->end());

	glPolygonMode (GL_FRONT_AND_BACK, GL_POINT);
	glPointSize(settings->pointSize);
	glMaterialf(GL_FRONT, GL_SHININESS, 20);
	glBegin(GL_POINTS);
	for(uint i = 0 ; i != acvt->vertex_number; i++)
	{
		float *v = acvt->mesh->getVertex(i);
		float extraVal = acvt->vtxQuadricsWeights->at(i);
		computeRGBcolorOfWeight(minVal, maxVal, extraVal, materialClr);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialClr);
		glNormal3fv(acvt->getOriginalNormal(i));
		glVertex3f(v[0]-ctr[0], v[1]-ctr[1], v[2]-ctr[2]);
	}
	glEnd();	
}