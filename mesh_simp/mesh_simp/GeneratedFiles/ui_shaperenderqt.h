/********************************************************************************
** Form generated from reading UI file 'shaperenderqt.ui'
**
** Created by: Qt User Interface Compiler version 5.3.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SHAPERENDERQT_H
#define UI_SHAPERENDERQT_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QScrollBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ShapeRenderQtClass
{
public:
    QAction *action_Open;
    QAction *action_Exit;
    QAction *action_Save;
    QAction *action_New;
    QAction *action_QEM;
    QAction *action_QEM_Optimization;
    QAction *action_Points;
    QAction *action_Edges;
    QAction *action_Hidden_Lines;
    QAction *action_Flat;
    QAction *action_Smooth;
    QAction *action_Textures;
    QWidget *centralWidget;
    QMenuBar *menuBar;
    QMenu *menu_File;
    QMenu *menu_View;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;
    QToolBar *viewToolBar;
    QDockWidget *dockWidget_QEM;
    QWidget *dockWidgetContents;
    QVBoxLayout *verticalLayout;
    QGroupBox *groupBox_CrestLine;
    QVBoxLayout *verticalLayout_2;
    QCheckBox *checkBox_ShowCrestLine;
    QLabel *label_Ridgeness;
    QScrollBar *horizontalScrollBar_Ridgeness;
    QLabel *label_Sphericalness;
    QScrollBar *horizontalScrollBar_Sphericalness;
    QLabel *label_Cyclideness;
    QScrollBar *horizontalScrollBar_Cyclideness;
    QCheckBox *checkBox_UseCrestLine;
    QLineEdit *lineEdit_eta;
    QCheckBox *checkBox_UseVertexQuality;
    QLineEdit *lineEdit_gamma;
    QLabel *label_VertexQualityThreshold;
    QLineEdit *lineEdit_VertexQualityThreshold;
    QPushButton *pushButton_ComputeWeights;
    QGroupBox *groupBox_QEM;
    QVBoxLayout *verticalLayout_3;
    QLabel *label_Vertex_Number;
    QLineEdit *lineEdit_FinalVertexNumber;
    QLabel *label_ReductionRatio;
    QLineEdit *lineEdit_ReductionRatio;
    QCheckBox *checkBox_ShowClusters;
    QPushButton *pushButton_AllSteps;
    QGroupBox *groupBox_SeparateSteps;
    QVBoxLayout *verticalLayout_4;
    QPushButton *pushButton_QEM;
    QPushButton *pushButton_QEM_Optimization;
    QPushButton *pushButton_SplitClusters;
    QPushButton *pushButton_ReQEM;
    QPushButton *pushButton_SplitNonManifold;
    QPushButton *pushButton_CreateConnectivity;
    QPushButton *pushButton_ReverseNormals;

    void setupUi(QMainWindow *ShapeRenderQtClass)
    {
        if (ShapeRenderQtClass->objectName().isEmpty())
            ShapeRenderQtClass->setObjectName(QStringLiteral("ShapeRenderQtClass"));
        ShapeRenderQtClass->resize(1100, 928);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(ShapeRenderQtClass->sizePolicy().hasHeightForWidth());
        ShapeRenderQtClass->setSizePolicy(sizePolicy);
        QIcon icon;
        icon.addFile(QStringLiteral(":/ShapeRenderQt/ref/square.png"), QSize(), QIcon::Normal, QIcon::Off);
        ShapeRenderQtClass->setWindowIcon(icon);
        action_Open = new QAction(ShapeRenderQtClass);
        action_Open->setObjectName(QStringLiteral("action_Open"));
        QIcon icon1;
        icon1.addFile(QStringLiteral(":/ShapeRenderQt/ref/open.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_Open->setIcon(icon1);
        action_Exit = new QAction(ShapeRenderQtClass);
        action_Exit->setObjectName(QStringLiteral("action_Exit"));
        QIcon icon2;
        icon2.addFile(QStringLiteral(":/ShapeRenderQt/ref/exit.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_Exit->setIcon(icon2);
        action_Save = new QAction(ShapeRenderQtClass);
        action_Save->setObjectName(QStringLiteral("action_Save"));
        QIcon icon3;
        icon3.addFile(QStringLiteral(":/ShapeRenderQt/ref/save.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_Save->setIcon(icon3);
        action_New = new QAction(ShapeRenderQtClass);
        action_New->setObjectName(QStringLiteral("action_New"));
        QIcon icon4;
        icon4.addFile(QStringLiteral(":/ShapeRenderQt/ref/new.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_New->setIcon(icon4);
        action_QEM = new QAction(ShapeRenderQtClass);
        action_QEM->setObjectName(QStringLiteral("action_QEM"));
        action_QEM_Optimization = new QAction(ShapeRenderQtClass);
        action_QEM_Optimization->setObjectName(QStringLiteral("action_QEM_Optimization"));
        action_Points = new QAction(ShapeRenderQtClass);
        action_Points->setObjectName(QStringLiteral("action_Points"));
        QIcon icon5;
        icon5.addFile(QStringLiteral(":/ShapeRenderQt/ref/points.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_Points->setIcon(icon5);
        action_Edges = new QAction(ShapeRenderQtClass);
        action_Edges->setObjectName(QStringLiteral("action_Edges"));
        QIcon icon6;
        icon6.addFile(QStringLiteral(":/ShapeRenderQt/ref/wire.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_Edges->setIcon(icon6);
        action_Hidden_Lines = new QAction(ShapeRenderQtClass);
        action_Hidden_Lines->setObjectName(QStringLiteral("action_Hidden_Lines"));
        QIcon icon7;
        icon7.addFile(QStringLiteral(":/ShapeRenderQt/ref/backlines.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_Hidden_Lines->setIcon(icon7);
        action_Flat = new QAction(ShapeRenderQtClass);
        action_Flat->setObjectName(QStringLiteral("action_Flat"));
        QIcon icon8;
        icon8.addFile(QStringLiteral(":/ShapeRenderQt/ref/flat.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_Flat->setIcon(icon8);
        action_Smooth = new QAction(ShapeRenderQtClass);
        action_Smooth->setObjectName(QStringLiteral("action_Smooth"));
        QIcon icon9;
        icon9.addFile(QStringLiteral(":/ShapeRenderQt/ref/smooth.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_Smooth->setIcon(icon9);
        action_Textures = new QAction(ShapeRenderQtClass);
        action_Textures->setObjectName(QStringLiteral("action_Textures"));
        QIcon icon10;
        icon10.addFile(QStringLiteral(":/ShapeRenderQt/ref/textures.png"), QSize(), QIcon::Normal, QIcon::Off);
        action_Textures->setIcon(icon10);
        centralWidget = new QWidget(ShapeRenderQtClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        ShapeRenderQtClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(ShapeRenderQtClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1100, 23));
        QFont font;
        font.setPointSize(10);
        menuBar->setFont(font);
        menu_File = new QMenu(menuBar);
        menu_File->setObjectName(QStringLiteral("menu_File"));
        menu_View = new QMenu(menuBar);
        menu_View->setObjectName(QStringLiteral("menu_View"));
        ShapeRenderQtClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(ShapeRenderQtClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        mainToolBar->setIconSize(QSize(36, 36));
        ShapeRenderQtClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(ShapeRenderQtClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        ShapeRenderQtClass->setStatusBar(statusBar);
        viewToolBar = new QToolBar(ShapeRenderQtClass);
        viewToolBar->setObjectName(QStringLiteral("viewToolBar"));
        viewToolBar->setIconSize(QSize(36, 36));
        ShapeRenderQtClass->addToolBar(Qt::TopToolBarArea, viewToolBar);
        dockWidget_QEM = new QDockWidget(ShapeRenderQtClass);
        dockWidget_QEM->setObjectName(QStringLiteral("dockWidget_QEM"));
        dockWidget_QEM->setEnabled(true);
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(dockWidget_QEM->sizePolicy().hasHeightForWidth());
        dockWidget_QEM->setSizePolicy(sizePolicy1);
        QFont font1;
        font1.setPointSize(9);
        dockWidget_QEM->setFont(font1);
        dockWidget_QEM->setAcceptDrops(false);
        dockWidget_QEM->setFeatures(QDockWidget::AllDockWidgetFeatures);
        dockWidgetContents = new QWidget();
        dockWidgetContents->setObjectName(QStringLiteral("dockWidgetContents"));
        verticalLayout = new QVBoxLayout(dockWidgetContents);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        groupBox_CrestLine = new QGroupBox(dockWidgetContents);
        groupBox_CrestLine->setObjectName(QStringLiteral("groupBox_CrestLine"));
        verticalLayout_2 = new QVBoxLayout(groupBox_CrestLine);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        checkBox_ShowCrestLine = new QCheckBox(groupBox_CrestLine);
        checkBox_ShowCrestLine->setObjectName(QStringLiteral("checkBox_ShowCrestLine"));

        verticalLayout_2->addWidget(checkBox_ShowCrestLine);

        label_Ridgeness = new QLabel(groupBox_CrestLine);
        label_Ridgeness->setObjectName(QStringLiteral("label_Ridgeness"));

        verticalLayout_2->addWidget(label_Ridgeness);

        horizontalScrollBar_Ridgeness = new QScrollBar(groupBox_CrestLine);
        horizontalScrollBar_Ridgeness->setObjectName(QStringLiteral("horizontalScrollBar_Ridgeness"));
        horizontalScrollBar_Ridgeness->setEnabled(false);
        horizontalScrollBar_Ridgeness->setMaximum(100);
        horizontalScrollBar_Ridgeness->setSingleStep(1);
        horizontalScrollBar_Ridgeness->setPageStep(10);
        horizontalScrollBar_Ridgeness->setValue(0);
        horizontalScrollBar_Ridgeness->setOrientation(Qt::Horizontal);

        verticalLayout_2->addWidget(horizontalScrollBar_Ridgeness);

        label_Sphericalness = new QLabel(groupBox_CrestLine);
        label_Sphericalness->setObjectName(QStringLiteral("label_Sphericalness"));

        verticalLayout_2->addWidget(label_Sphericalness);

        horizontalScrollBar_Sphericalness = new QScrollBar(groupBox_CrestLine);
        horizontalScrollBar_Sphericalness->setObjectName(QStringLiteral("horizontalScrollBar_Sphericalness"));
        horizontalScrollBar_Sphericalness->setEnabled(false);
        horizontalScrollBar_Sphericalness->setMaximum(100);
        horizontalScrollBar_Sphericalness->setOrientation(Qt::Horizontal);

        verticalLayout_2->addWidget(horizontalScrollBar_Sphericalness);

        label_Cyclideness = new QLabel(groupBox_CrestLine);
        label_Cyclideness->setObjectName(QStringLiteral("label_Cyclideness"));

        verticalLayout_2->addWidget(label_Cyclideness);

        horizontalScrollBar_Cyclideness = new QScrollBar(groupBox_CrestLine);
        horizontalScrollBar_Cyclideness->setObjectName(QStringLiteral("horizontalScrollBar_Cyclideness"));
        horizontalScrollBar_Cyclideness->setEnabled(false);
        horizontalScrollBar_Cyclideness->setMaximum(100);
        horizontalScrollBar_Cyclideness->setOrientation(Qt::Horizontal);

        verticalLayout_2->addWidget(horizontalScrollBar_Cyclideness);

        checkBox_UseCrestLine = new QCheckBox(groupBox_CrestLine);
        checkBox_UseCrestLine->setObjectName(QStringLiteral("checkBox_UseCrestLine"));

        verticalLayout_2->addWidget(checkBox_UseCrestLine);

        lineEdit_eta = new QLineEdit(groupBox_CrestLine);
        lineEdit_eta->setObjectName(QStringLiteral("lineEdit_eta"));

        verticalLayout_2->addWidget(lineEdit_eta);

        checkBox_UseVertexQuality = new QCheckBox(groupBox_CrestLine);
        checkBox_UseVertexQuality->setObjectName(QStringLiteral("checkBox_UseVertexQuality"));

        verticalLayout_2->addWidget(checkBox_UseVertexQuality);

        lineEdit_gamma = new QLineEdit(groupBox_CrestLine);
        lineEdit_gamma->setObjectName(QStringLiteral("lineEdit_gamma"));

        verticalLayout_2->addWidget(lineEdit_gamma);

        label_VertexQualityThreshold = new QLabel(groupBox_CrestLine);
        label_VertexQualityThreshold->setObjectName(QStringLiteral("label_VertexQualityThreshold"));

        verticalLayout_2->addWidget(label_VertexQualityThreshold);

        lineEdit_VertexQualityThreshold = new QLineEdit(groupBox_CrestLine);
        lineEdit_VertexQualityThreshold->setObjectName(QStringLiteral("lineEdit_VertexQualityThreshold"));

        verticalLayout_2->addWidget(lineEdit_VertexQualityThreshold);

        pushButton_ComputeWeights = new QPushButton(groupBox_CrestLine);
        pushButton_ComputeWeights->setObjectName(QStringLiteral("pushButton_ComputeWeights"));

        verticalLayout_2->addWidget(pushButton_ComputeWeights);


        verticalLayout->addWidget(groupBox_CrestLine);

        groupBox_QEM = new QGroupBox(dockWidgetContents);
        groupBox_QEM->setObjectName(QStringLiteral("groupBox_QEM"));
        verticalLayout_3 = new QVBoxLayout(groupBox_QEM);
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setContentsMargins(11, 11, 11, 11);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        label_Vertex_Number = new QLabel(groupBox_QEM);
        label_Vertex_Number->setObjectName(QStringLiteral("label_Vertex_Number"));

        verticalLayout_3->addWidget(label_Vertex_Number);

        lineEdit_FinalVertexNumber = new QLineEdit(groupBox_QEM);
        lineEdit_FinalVertexNumber->setObjectName(QStringLiteral("lineEdit_FinalVertexNumber"));

        verticalLayout_3->addWidget(lineEdit_FinalVertexNumber);

        label_ReductionRatio = new QLabel(groupBox_QEM);
        label_ReductionRatio->setObjectName(QStringLiteral("label_ReductionRatio"));

        verticalLayout_3->addWidget(label_ReductionRatio);

        lineEdit_ReductionRatio = new QLineEdit(groupBox_QEM);
        lineEdit_ReductionRatio->setObjectName(QStringLiteral("lineEdit_ReductionRatio"));

        verticalLayout_3->addWidget(lineEdit_ReductionRatio);

        checkBox_ShowClusters = new QCheckBox(groupBox_QEM);
        checkBox_ShowClusters->setObjectName(QStringLiteral("checkBox_ShowClusters"));
        checkBox_ShowClusters->setChecked(true);

        verticalLayout_3->addWidget(checkBox_ShowClusters);

        pushButton_AllSteps = new QPushButton(groupBox_QEM);
        pushButton_AllSteps->setObjectName(QStringLiteral("pushButton_AllSteps"));

        verticalLayout_3->addWidget(pushButton_AllSteps);

        groupBox_SeparateSteps = new QGroupBox(groupBox_QEM);
        groupBox_SeparateSteps->setObjectName(QStringLiteral("groupBox_SeparateSteps"));
        verticalLayout_4 = new QVBoxLayout(groupBox_SeparateSteps);
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setContentsMargins(11, 11, 11, 11);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        pushButton_QEM = new QPushButton(groupBox_SeparateSteps);
        pushButton_QEM->setObjectName(QStringLiteral("pushButton_QEM"));
        pushButton_QEM->setFlat(false);

        verticalLayout_4->addWidget(pushButton_QEM);

        pushButton_QEM_Optimization = new QPushButton(groupBox_SeparateSteps);
        pushButton_QEM_Optimization->setObjectName(QStringLiteral("pushButton_QEM_Optimization"));

        verticalLayout_4->addWidget(pushButton_QEM_Optimization);

        pushButton_SplitClusters = new QPushButton(groupBox_SeparateSteps);
        pushButton_SplitClusters->setObjectName(QStringLiteral("pushButton_SplitClusters"));

        verticalLayout_4->addWidget(pushButton_SplitClusters);

        pushButton_ReQEM = new QPushButton(groupBox_SeparateSteps);
        pushButton_ReQEM->setObjectName(QStringLiteral("pushButton_ReQEM"));

        verticalLayout_4->addWidget(pushButton_ReQEM);

        pushButton_SplitNonManifold = new QPushButton(groupBox_SeparateSteps);
        pushButton_SplitNonManifold->setObjectName(QStringLiteral("pushButton_SplitNonManifold"));

        verticalLayout_4->addWidget(pushButton_SplitNonManifold);

        pushButton_CreateConnectivity = new QPushButton(groupBox_SeparateSteps);
        pushButton_CreateConnectivity->setObjectName(QStringLiteral("pushButton_CreateConnectivity"));

        verticalLayout_4->addWidget(pushButton_CreateConnectivity);

        pushButton_ReverseNormals = new QPushButton(groupBox_SeparateSteps);
        pushButton_ReverseNormals->setObjectName(QStringLiteral("pushButton_ReverseNormals"));

        verticalLayout_4->addWidget(pushButton_ReverseNormals);


        verticalLayout_3->addWidget(groupBox_SeparateSteps);


        verticalLayout->addWidget(groupBox_QEM);

        dockWidget_QEM->setWidget(dockWidgetContents);
        ShapeRenderQtClass->addDockWidget(static_cast<Qt::DockWidgetArea>(2), dockWidget_QEM);

        menuBar->addAction(menu_File->menuAction());
        menuBar->addAction(menu_View->menuAction());
        menu_File->addAction(action_New);
        menu_File->addAction(action_Open);
        menu_File->addAction(action_Save);
        menu_File->addSeparator();
        menu_File->addAction(action_Exit);
        menu_View->addAction(action_Points);
        menu_View->addAction(action_Edges);
        menu_View->addAction(action_Hidden_Lines);
        menu_View->addAction(action_Flat);
        menu_View->addAction(action_Smooth);
        menu_View->addAction(action_Textures);
        mainToolBar->addAction(action_New);
        mainToolBar->addAction(action_Open);
        mainToolBar->addAction(action_Save);
        viewToolBar->addAction(action_Points);
        viewToolBar->addAction(action_Edges);
        viewToolBar->addAction(action_Hidden_Lines);
        viewToolBar->addAction(action_Flat);
        viewToolBar->addAction(action_Smooth);
        viewToolBar->addAction(action_Textures);

        retranslateUi(ShapeRenderQtClass);

        QMetaObject::connectSlotsByName(ShapeRenderQtClass);
    } // setupUi

    void retranslateUi(QMainWindow *ShapeRenderQtClass)
    {
        ShapeRenderQtClass->setWindowTitle(QApplication::translate("ShapeRenderQtClass", "Mesh Simplification", 0));
        action_Open->setText(QApplication::translate("ShapeRenderQtClass", "&Open", 0));
        action_Open->setShortcut(QApplication::translate("ShapeRenderQtClass", "Ctrl+O", 0));
        action_Exit->setText(QApplication::translate("ShapeRenderQtClass", "&Exit", 0));
        action_Exit->setShortcut(QApplication::translate("ShapeRenderQtClass", "Ctrl+E", 0));
        action_Save->setText(QApplication::translate("ShapeRenderQtClass", "&Save", 0));
        action_Save->setShortcut(QApplication::translate("ShapeRenderQtClass", "Ctrl+S", 0));
        action_New->setText(QApplication::translate("ShapeRenderQtClass", "&New", 0));
        action_New->setShortcut(QApplication::translate("ShapeRenderQtClass", "Ctrl+N", 0));
        action_QEM->setText(QApplication::translate("ShapeRenderQtClass", "QEM", 0));
        action_QEM_Optimization->setText(QApplication::translate("ShapeRenderQtClass", "QEM Optimization", 0));
        action_Points->setText(QApplication::translate("ShapeRenderQtClass", "&Points", 0));
        action_Edges->setText(QApplication::translate("ShapeRenderQtClass", "&Edges", 0));
        action_Hidden_Lines->setText(QApplication::translate("ShapeRenderQtClass", "&Hidden Lines", 0));
        action_Flat->setText(QApplication::translate("ShapeRenderQtClass", "&Flat", 0));
        action_Smooth->setText(QApplication::translate("ShapeRenderQtClass", "&Smooth", 0));
        action_Textures->setText(QApplication::translate("ShapeRenderQtClass", "&Textures", 0));
        menu_File->setTitle(QApplication::translate("ShapeRenderQtClass", "&File", 0));
        menu_View->setTitle(QApplication::translate("ShapeRenderQtClass", "&View", 0));
        mainToolBar->setWindowTitle(QApplication::translate("ShapeRenderQtClass", "File", 0));
        viewToolBar->setWindowTitle(QApplication::translate("ShapeRenderQtClass", "View", 0));
        dockWidget_QEM->setWindowTitle(QApplication::translate("ShapeRenderQtClass", "Mesh Simplification Panel", 0));
        groupBox_CrestLine->setTitle(QApplication::translate("ShapeRenderQtClass", "Crest Lines and Vertex Quality Confidence", 0));
        checkBox_ShowCrestLine->setText(QApplication::translate("ShapeRenderQtClass", "Show Crest Lines", 0));
        label_Ridgeness->setText(QApplication::translate("ShapeRenderQtClass", "Ridgeness", 0));
        label_Sphericalness->setText(QApplication::translate("ShapeRenderQtClass", "Sphericalness", 0));
        label_Cyclideness->setText(QApplication::translate("ShapeRenderQtClass", "Cyclideness", 0));
        checkBox_UseCrestLine->setText(QApplication::translate("ShapeRenderQtClass", "Use Crest Line Weight: eta = ", 0));
        lineEdit_eta->setText(QApplication::translate("ShapeRenderQtClass", "3", 0));
        checkBox_UseVertexQuality->setText(QApplication::translate("ShapeRenderQtClass", "Use Vertex Quality Weight: gamma = ", 0));
        lineEdit_gamma->setText(QApplication::translate("ShapeRenderQtClass", "3", 0));
        label_VertexQualityThreshold->setText(QApplication::translate("ShapeRenderQtClass", "Vertex Quality Threshold", 0));
        lineEdit_VertexQualityThreshold->setText(QApplication::translate("ShapeRenderQtClass", "0.05", 0));
        pushButton_ComputeWeights->setText(QApplication::translate("ShapeRenderQtClass", "Compute Vertex Weights", 0));
        groupBox_QEM->setTitle(QApplication::translate("ShapeRenderQtClass", "Mesh Simplification", 0));
        label_Vertex_Number->setText(QApplication::translate("ShapeRenderQtClass", "Final Vertex Number", 0));
        lineEdit_FinalVertexNumber->setText(QApplication::translate("ShapeRenderQtClass", "1000", 0));
        label_ReductionRatio->setText(QApplication::translate("ShapeRenderQtClass", "Reduction Ratio (0-1)", 0));
        lineEdit_ReductionRatio->setText(QApplication::translate("ShapeRenderQtClass", "0.01", 0));
        checkBox_ShowClusters->setText(QApplication::translate("ShapeRenderQtClass", "Show Clusters", 0));
        pushButton_AllSteps->setText(QApplication::translate("ShapeRenderQtClass", "QEM-based Simplification", 0));
        groupBox_SeparateSteps->setTitle(QApplication::translate("ShapeRenderQtClass", "Separate Steps", 0));
        pushButton_QEM->setText(QApplication::translate("ShapeRenderQtClass", "QSlim", 0));
        pushButton_QEM_Optimization->setText(QApplication::translate("ShapeRenderQtClass", "Swap", 0));
        pushButton_SplitClusters->setText(QApplication::translate("ShapeRenderQtClass", "Split Clusters", 0));
        pushButton_ReQEM->setText(QApplication::translate("ShapeRenderQtClass", "Re-QEM After Split", 0));
        pushButton_SplitNonManifold->setText(QApplication::translate("ShapeRenderQtClass", "Split Non-Manifold Clusters", 0));
        pushButton_CreateConnectivity->setText(QApplication::translate("ShapeRenderQtClass", "Create New Connectivity", 0));
        pushButton_ReverseNormals->setText(QApplication::translate("ShapeRenderQtClass", "Reverse Normals", 0));
    } // retranslateUi

};

namespace Ui {
    class ShapeRenderQtClass: public Ui_ShapeRenderQtClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SHAPERENDERQT_H
