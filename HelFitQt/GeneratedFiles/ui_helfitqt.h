/********************************************************************************
** Form generated from reading UI file 'helfitqt.ui'
**
** Created by: Qt User Interface Compiler version 5.6.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_HELFITQT_H
#define UI_HELFITQT_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_HelFitQtClass
{
public:
    QAction *menuFileExit;
    QAction *menuLoadModel;
    QAction *menuDebyeCurrModel;
    QAction *lblcurveweight;
    QAction *menuWeight0;
    QAction *menuWeight1;
    QAction *menuWeight2;
    QAction *menuRefit;
    QAction *menuGenerate_random_seed;
    QAction *actionChange_fitting_params;
    QWidget *centralWidget;
    QTabWidget *tabWidget;
    QWidget *tabDataView;
    QCustomPlot *wdgtDataPlot;
    QCheckBox *chkbLogLogPlot;
    QWidget *tab3Dview;
    QSlider *hsliderPhiRot;
    QLabel *label_12;
    QLabel *label_13;
    QCustomPlot *wdgtXZplot;
    QCustomPlot *wdgtYZplot;
    QCustomPlot *wdgtXYplot;
    QLabel *label_14;
    QLabel *label_15;
    QWidget *tabTextoutput;
    QTextEdit *textOutput;
    QGroupBox *groupBox_1;
    QTextEdit *txtFilePath;
    QLabel *label;
    QPushButton *btnLoadFile;
    QSpinBox *spbDataMin;
    QSpinBox *spbDataMax;
    QLabel *label_2;
    QLabel *label_3;
    QLabel *label_4;
    QLabel *label_5;
    QLabel *label_6;
    QSpinBox *spbFitMin;
    QSpinBox *spbFitMax;
    QLabel *label_7;
    QLabel *label_8;
    QLabel *label_9;
    QLineEdit *lineQmin;
    QLineEdit *lineQmax;
    QLabel *label_10;
    QLabel *label_11;
    QGroupBox *groupBox_2;
    QSpinBox *spbCalcNrPoints;
    QLabel *label_16;
    QLineEdit *lineCalcQmin;
    QLabel *label_17;
    QLabel *label_18;
    QLineEdit *lineCalcQmax;
    QLabel *label_19;
    QLabel *label_20;
    QLabel *label_21;
    QSpinBox *spbCalcNrStacks;
    QLabel *label_22;
    QSpinBox *spbNrCores;
    QLineEdit *lineStackSpacing;
    QLabel *label_23;
    QLabel *label_24;
    QPushButton *btnFit;
    QMenuBar *menuBar;
    QMenu *menu_File;
    QMenu *menuRunDebyeCalc;
    QMenu *menuOptions;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *HelFitQtClass)
    {
        if (HelFitQtClass->objectName().isEmpty())
            HelFitQtClass->setObjectName(QStringLiteral("HelFitQtClass"));
        HelFitQtClass->setEnabled(true);
        HelFitQtClass->resize(1024, 639);
        HelFitQtClass->setAnimated(true);
        menuFileExit = new QAction(HelFitQtClass);
        menuFileExit->setObjectName(QStringLiteral("menuFileExit"));
        menuLoadModel = new QAction(HelFitQtClass);
        menuLoadModel->setObjectName(QStringLiteral("menuLoadModel"));
        menuDebyeCurrModel = new QAction(HelFitQtClass);
        menuDebyeCurrModel->setObjectName(QStringLiteral("menuDebyeCurrModel"));
        menuDebyeCurrModel->setEnabled(false);
        lblcurveweight = new QAction(HelFitQtClass);
        lblcurveweight->setObjectName(QStringLiteral("lblcurveweight"));
        lblcurveweight->setEnabled(false);
        menuWeight0 = new QAction(HelFitQtClass);
        menuWeight0->setObjectName(QStringLiteral("menuWeight0"));
        menuWeight0->setCheckable(true);
        menuWeight1 = new QAction(HelFitQtClass);
        menuWeight1->setObjectName(QStringLiteral("menuWeight1"));
        menuWeight1->setCheckable(true);
        menuWeight1->setChecked(true);
        menuWeight2 = new QAction(HelFitQtClass);
        menuWeight2->setObjectName(QStringLiteral("menuWeight2"));
        menuWeight2->setCheckable(true);
        menuRefit = new QAction(HelFitQtClass);
        menuRefit->setObjectName(QStringLiteral("menuRefit"));
        menuRefit->setEnabled(false);
        menuGenerate_random_seed = new QAction(HelFitQtClass);
        menuGenerate_random_seed->setObjectName(QStringLiteral("menuGenerate_random_seed"));
        actionChange_fitting_params = new QAction(HelFitQtClass);
        actionChange_fitting_params->setObjectName(QStringLiteral("actionChange_fitting_params"));
        centralWidget = new QWidget(HelFitQtClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        tabWidget = new QTabWidget(centralWidget);
        tabWidget->setObjectName(QStringLiteral("tabWidget"));
        tabWidget->setEnabled(true);
        tabWidget->setGeometry(QRect(4, 0, 761, 531));
        tabWidget->setLayoutDirection(Qt::LeftToRight);
        tabWidget->setAutoFillBackground(false);
        tabWidget->setTabPosition(QTabWidget::North);
        tabWidget->setTabShape(QTabWidget::Rounded);
        tabWidget->setTabsClosable(false);
        tabDataView = new QWidget();
        tabDataView->setObjectName(QStringLiteral("tabDataView"));
        wdgtDataPlot = new QCustomPlot(tabDataView);
        wdgtDataPlot->setObjectName(QStringLiteral("wdgtDataPlot"));
        wdgtDataPlot->setGeometry(QRect(9, 9, 731, 491));
        chkbLogLogPlot = new QCheckBox(wdgtDataPlot);
        chkbLogLogPlot->setObjectName(QStringLiteral("chkbLogLogPlot"));
        chkbLogLogPlot->setEnabled(false);
        chkbLogLogPlot->setGeometry(QRect(660, 0, 81, 20));
        chkbLogLogPlot->setChecked(false);
        tabWidget->addTab(tabDataView, QString());
        tab3Dview = new QWidget();
        tab3Dview->setObjectName(QStringLiteral("tab3Dview"));
        hsliderPhiRot = new QSlider(tab3Dview);
        hsliderPhiRot->setObjectName(QStringLiteral("hsliderPhiRot"));
        hsliderPhiRot->setEnabled(true);
        hsliderPhiRot->setGeometry(QRect(70, 430, 591, 22));
        hsliderPhiRot->setMaximum(1000);
        hsliderPhiRot->setValue(500);
        hsliderPhiRot->setOrientation(Qt::Horizontal);
        hsliderPhiRot->setInvertedAppearance(false);
        hsliderPhiRot->setInvertedControls(false);
        hsliderPhiRot->setTickPosition(QSlider::TicksBelow);
        hsliderPhiRot->setTickInterval(250);
        label_12 = new QLabel(tab3Dview);
        label_12->setObjectName(QStringLiteral("label_12"));
        label_12->setGeometry(QRect(310, 400, 121, 16));
        QFont font;
        font.setBold(true);
        font.setWeight(75);
        label_12->setFont(font);
        label_13 = new QLabel(tab3Dview);
        label_13->setObjectName(QStringLiteral("label_13"));
        label_13->setGeometry(QRect(70, 460, 611, 16));
        wdgtXZplot = new QCustomPlot(tab3Dview);
        wdgtXZplot->setObjectName(QStringLiteral("wdgtXZplot"));
        wdgtXZplot->setGeometry(QRect(30, 70, 200, 250));
        wdgtYZplot = new QCustomPlot(tab3Dview);
        wdgtYZplot->setObjectName(QStringLiteral("wdgtYZplot"));
        wdgtYZplot->setGeometry(QRect(280, 70, 200, 250));
        wdgtXYplot = new QCustomPlot(tab3Dview);
        wdgtXYplot->setObjectName(QStringLiteral("wdgtXYplot"));
        wdgtXYplot->setGeometry(QRect(520, 90, 200, 200));
        label_14 = new QLabel(tab3Dview);
        label_14->setObjectName(QStringLiteral("label_14"));
        label_14->setGeometry(QRect(110, 50, 571, 20));
        label_14->setFont(font);
        label_15 = new QLabel(tab3Dview);
        label_15->setObjectName(QStringLiteral("label_15"));
        label_15->setGeometry(QRect(190, 20, 511, 20));
        QFont font1;
        font1.setPointSize(10);
        font1.setBold(true);
        font1.setWeight(75);
        label_15->setFont(font1);
        tabWidget->addTab(tab3Dview, QString());
        tabTextoutput = new QWidget();
        tabTextoutput->setObjectName(QStringLiteral("tabTextoutput"));
        textOutput = new QTextEdit(tabTextoutput);
        textOutput->setObjectName(QStringLiteral("textOutput"));
        textOutput->setGeometry(QRect(10, 10, 731, 481));
        textOutput->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
        textOutput->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        textOutput->setReadOnly(true);
        tabWidget->addTab(tabTextoutput, QString());
        groupBox_1 = new QGroupBox(centralWidget);
        groupBox_1->setObjectName(QStringLiteral("groupBox_1"));
        groupBox_1->setEnabled(true);
        groupBox_1->setGeometry(QRect(780, 20, 231, 241));
        groupBox_1->setAutoFillBackground(true);
        groupBox_1->setStyleSheet(QStringLiteral(""));
        groupBox_1->setFlat(false);
        groupBox_1->setCheckable(false);
        txtFilePath = new QTextEdit(groupBox_1);
        txtFilePath->setObjectName(QStringLiteral("txtFilePath"));
        txtFilePath->setEnabled(true);
        txtFilePath->setGeometry(QRect(30, 50, 171, 31));
        label = new QLabel(groupBox_1);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(30, 30, 131, 16));
        btnLoadFile = new QPushButton(groupBox_1);
        btnLoadFile->setObjectName(QStringLiteral("btnLoadFile"));
        btnLoadFile->setGeometry(QRect(160, 27, 51, 21));
        spbDataMin = new QSpinBox(groupBox_1);
        spbDataMin->setObjectName(QStringLiteral("spbDataMin"));
        spbDataMin->setEnabled(false);
        spbDataMin->setGeometry(QRect(60, 110, 42, 22));
        spbDataMax = new QSpinBox(groupBox_1);
        spbDataMax->setObjectName(QStringLiteral("spbDataMax"));
        spbDataMax->setEnabled(false);
        spbDataMax->setGeometry(QRect(160, 110, 42, 22));
        spbDataMax->setMaximum(99999);
        spbDataMax->setValue(999);
        label_2 = new QLabel(groupBox_1);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(30, 90, 171, 16));
        label_3 = new QLabel(groupBox_1);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(30, 110, 31, 16));
        label_4 = new QLabel(groupBox_1);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(120, 110, 31, 16));
        label_5 = new QLabel(groupBox_1);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setGeometry(QRect(30, 160, 31, 16));
        label_6 = new QLabel(groupBox_1);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setGeometry(QRect(30, 140, 171, 16));
        spbFitMin = new QSpinBox(groupBox_1);
        spbFitMin->setObjectName(QStringLiteral("spbFitMin"));
        spbFitMin->setEnabled(false);
        spbFitMin->setGeometry(QRect(60, 160, 42, 22));
        spbFitMax = new QSpinBox(groupBox_1);
        spbFitMax->setObjectName(QStringLiteral("spbFitMax"));
        spbFitMax->setEnabled(false);
        spbFitMax->setGeometry(QRect(160, 160, 42, 22));
        spbFitMax->setMaximum(99999);
        spbFitMax->setValue(999);
        label_7 = new QLabel(groupBox_1);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setGeometry(QRect(120, 160, 31, 16));
        label_8 = new QLabel(groupBox_1);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setGeometry(QRect(40, 190, 50, 16));
        label_8->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_9 = new QLabel(groupBox_1);
        label_9->setObjectName(QStringLiteral("label_9"));
        label_9->setGeometry(QRect(40, 210, 50, 20));
        label_9->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        lineQmin = new QLineEdit(groupBox_1);
        lineQmin->setObjectName(QStringLiteral("lineQmin"));
        lineQmin->setEnabled(false);
        lineQmin->setGeometry(QRect(90, 190, 61, 22));
        lineQmax = new QLineEdit(groupBox_1);
        lineQmax->setObjectName(QStringLiteral("lineQmax"));
        lineQmax->setEnabled(false);
        lineQmax->setGeometry(QRect(90, 210, 61, 22));
        label_10 = new QLabel(groupBox_1);
        label_10->setObjectName(QStringLiteral("label_10"));
        label_10->setGeometry(QRect(160, 191, 51, 16));
        label_11 = new QLabel(groupBox_1);
        label_11->setObjectName(QStringLiteral("label_11"));
        label_11->setGeometry(QRect(160, 212, 51, 16));
        groupBox_2 = new QGroupBox(centralWidget);
        groupBox_2->setObjectName(QStringLiteral("groupBox_2"));
        groupBox_2->setGeometry(QRect(780, 270, 231, 181));
        groupBox_2->setAutoFillBackground(true);
        spbCalcNrPoints = new QSpinBox(groupBox_2);
        spbCalcNrPoints->setObjectName(QStringLiteral("spbCalcNrPoints"));
        spbCalcNrPoints->setEnabled(true);
        spbCalcNrPoints->setGeometry(QRect(95, 60, 60, 22));
        spbCalcNrPoints->setMinimum(16);
        spbCalcNrPoints->setMaximum(512);
        spbCalcNrPoints->setValue(256);
        label_16 = new QLabel(groupBox_2);
        label_16->setObjectName(QStringLiteral("label_16"));
        label_16->setGeometry(QRect(165, 42, 51, 16));
        lineCalcQmin = new QLineEdit(groupBox_2);
        lineCalcQmin->setObjectName(QStringLiteral("lineCalcQmin"));
        lineCalcQmin->setEnabled(true);
        lineCalcQmin->setGeometry(QRect(95, 20, 60, 22));
        label_17 = new QLabel(groupBox_2);
        label_17->setObjectName(QStringLiteral("label_17"));
        label_17->setGeometry(QRect(40, 20, 50, 16));
        label_17->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_18 = new QLabel(groupBox_2);
        label_18->setObjectName(QStringLiteral("label_18"));
        label_18->setGeometry(QRect(165, 21, 51, 16));
        lineCalcQmax = new QLineEdit(groupBox_2);
        lineCalcQmax->setObjectName(QStringLiteral("lineCalcQmax"));
        lineCalcQmax->setEnabled(true);
        lineCalcQmax->setGeometry(QRect(95, 40, 60, 22));
        label_19 = new QLabel(groupBox_2);
        label_19->setObjectName(QStringLiteral("label_19"));
        label_19->setGeometry(QRect(40, 40, 50, 20));
        label_19->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_20 = new QLabel(groupBox_2);
        label_20->setObjectName(QStringLiteral("label_20"));
        label_20->setGeometry(QRect(40, 60, 50, 20));
        label_20->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_21 = new QLabel(groupBox_2);
        label_21->setObjectName(QStringLiteral("label_21"));
        label_21->setGeometry(QRect(20, 130, 70, 20));
        label_21->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spbCalcNrStacks = new QSpinBox(groupBox_2);
        spbCalcNrStacks->setObjectName(QStringLiteral("spbCalcNrStacks"));
        spbCalcNrStacks->setEnabled(true);
        spbCalcNrStacks->setGeometry(QRect(95, 130, 60, 22));
        spbCalcNrStacks->setMinimum(0);
        spbCalcNrStacks->setMaximum(50);
        spbCalcNrStacks->setValue(15);
        label_22 = new QLabel(groupBox_2);
        label_22->setObjectName(QStringLiteral("label_22"));
        label_22->setGeometry(QRect(20, 150, 70, 20));
        label_22->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spbNrCores = new QSpinBox(groupBox_2);
        spbNrCores->setObjectName(QStringLiteral("spbNrCores"));
        spbNrCores->setEnabled(true);
        spbNrCores->setGeometry(QRect(95, 150, 60, 22));
        spbNrCores->setMinimum(1);
        spbNrCores->setMaximum(500);
        spbNrCores->setValue(500);
        lineStackSpacing = new QLineEdit(groupBox_2);
        lineStackSpacing->setObjectName(QStringLiteral("lineStackSpacing"));
        lineStackSpacing->setEnabled(true);
        lineStackSpacing->setGeometry(QRect(95, 110, 60, 22));
        label_23 = new QLabel(groupBox_2);
        label_23->setObjectName(QStringLiteral("label_23"));
        label_23->setGeometry(QRect(20, 110, 70, 20));
        label_23->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_24 = new QLabel(groupBox_2);
        label_24->setObjectName(QStringLiteral("label_24"));
        label_24->setGeometry(QRect(165, 110, 51, 16));
        btnFit = new QPushButton(centralWidget);
        btnFit->setObjectName(QStringLiteral("btnFit"));
        btnFit->setGeometry(QRect(830, 480, 151, 71));
        HelFitQtClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(HelFitQtClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setEnabled(true);
        menuBar->setGeometry(QRect(0, 0, 1024, 26));
        menu_File = new QMenu(menuBar);
        menu_File->setObjectName(QStringLiteral("menu_File"));
        menuRunDebyeCalc = new QMenu(menuBar);
        menuRunDebyeCalc->setObjectName(QStringLiteral("menuRunDebyeCalc"));
        menuRunDebyeCalc->setEnabled(true);
        menuOptions = new QMenu(menuBar);
        menuOptions->setObjectName(QStringLiteral("menuOptions"));
        HelFitQtClass->setMenuBar(menuBar);
        statusBar = new QStatusBar(HelFitQtClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        HelFitQtClass->setStatusBar(statusBar);
#ifndef QT_NO_SHORTCUT
        label_3->setBuddy(spbDataMin);
        label_4->setBuddy(spbDataMax);
        label_5->setBuddy(spbFitMin);
        label_7->setBuddy(spbFitMax);
        label_8->setBuddy(lineQmin);
        label_9->setBuddy(lineQmax);
        label_10->setBuddy(lineQmin);
        label_11->setBuddy(lineQmax);
        label_16->setBuddy(lineCalcQmax);
        label_17->setBuddy(lineCalcQmin);
        label_18->setBuddy(lineCalcQmin);
        label_19->setBuddy(lineCalcQmax);
        label_20->setBuddy(spbCalcNrPoints);
        label_21->setBuddy(spbCalcNrStacks);
        label_22->setBuddy(spbNrCores);
        label_23->setBuddy(lineStackSpacing);
        label_24->setBuddy(lineStackSpacing);
#endif // QT_NO_SHORTCUT

        menuBar->addAction(menu_File->menuAction());
        menuBar->addAction(menuRunDebyeCalc->menuAction());
        menuBar->addAction(menuOptions->menuAction());
        menu_File->addAction(menuFileExit);
        menuRunDebyeCalc->addAction(menuLoadModel);
        menuRunDebyeCalc->addAction(menuGenerate_random_seed);
        menuRunDebyeCalc->addSeparator();
        menuRunDebyeCalc->addAction(menuDebyeCurrModel);
        menuRunDebyeCalc->addAction(menuRefit);
        menuOptions->addAction(lblcurveweight);
        menuOptions->addAction(menuWeight0);
        menuOptions->addAction(menuWeight1);
        menuOptions->addAction(menuWeight2);
        menuOptions->addSeparator();
        menuOptions->addAction(actionChange_fitting_params);

        retranslateUi(HelFitQtClass);
        QObject::connect(menuFileExit, SIGNAL(triggered()), HelFitQtClass, SLOT(close()));
        QObject::connect(chkbLogLogPlot, SIGNAL(clicked()), HelFitQtClass, SLOT(replotData()));
        QObject::connect(menuLoadModel, SIGNAL(triggered()), HelFitQtClass, SLOT(loadModelFromData()));
        QObject::connect(spbDataMin, SIGNAL(valueChanged(int)), HelFitQtClass, SLOT(changedMinDataRange()));
        QObject::connect(spbDataMax, SIGNAL(valueChanged(int)), HelFitQtClass, SLOT(changedMaxDataRange()));
        QObject::connect(spbFitMin, SIGNAL(valueChanged(int)), HelFitQtClass, SLOT(changedMinFitRange()));
        QObject::connect(spbFitMax, SIGNAL(valueChanged(int)), HelFitQtClass, SLOT(changedMaxFitRange()));
        QObject::connect(hsliderPhiRot, SIGNAL(valueChanged(int)), HelFitQtClass, SLOT(changedModelRotationPhi(int)));
        QObject::connect(menuDebyeCurrModel, SIGNAL(triggered()), HelFitQtClass, SLOT(calcDebyeStackCurrent()));
        QObject::connect(menuWeight0, SIGNAL(triggered()), HelFitQtClass, SLOT(actionWeight0()));
        QObject::connect(menuWeight1, SIGNAL(triggered()), HelFitQtClass, SLOT(actionWeight1()));
        QObject::connect(menuWeight2, SIGNAL(triggered()), HelFitQtClass, SLOT(actionWeight2()));
        QObject::connect(menuRefit, SIGNAL(triggered()), HelFitQtClass, SLOT(refitDebyeStackCurrent()));
        QObject::connect(menuGenerate_random_seed, SIGNAL(triggered()), HelFitQtClass, SLOT(actionGenerateRandModel()));

        tabWidget->setCurrentIndex(2);


        QMetaObject::connectSlotsByName(HelFitQtClass);
    } // setupUi

    void retranslateUi(QMainWindow *HelFitQtClass)
    {
        HelFitQtClass->setWindowTitle(QApplication::translate("HelFitQtClass", "HelFitQt", 0));
        menuFileExit->setText(QApplication::translate("HelFitQtClass", "Exit", 0));
        menuLoadModel->setText(QApplication::translate("HelFitQtClass", "Load existing model...", 0));
        menuDebyeCurrModel->setText(QApplication::translate("HelFitQtClass", "Calc. current model!", 0));
        lblcurveweight->setText(QApplication::translate("HelFitQtClass", "Curve weight...", 0));
        menuWeight0->setText(QApplication::translate("HelFitQtClass", "... I(q) * q", 0));
        menuWeight1->setText(QApplication::translate("HelFitQtClass", "... I(q) * q^2", 0));
        menuWeight2->setText(QApplication::translate("HelFitQtClass", "... ln(I(q))", 0));
        menuRefit->setText(QApplication::translate("HelFitQtClass", "Refit current modeldata", 0));
        menuGenerate_random_seed->setText(QApplication::translate("HelFitQtClass", "Generate random model...", 0));
        actionChange_fitting_params->setText(QApplication::translate("HelFitQtClass", "Change fitting params...", 0));
        chkbLogLogPlot->setText(QApplication::translate("HelFitQtClass", "log-log?", 0));
        tabWidget->setTabText(tabWidget->indexOf(tabDataView), QApplication::translate("HelFitQtClass", "Scattering Data", 0));
        label_12->setText(QApplication::translate("HelFitQtClass", "Model Rotation:", 0));
        label_13->setText(QApplication::translate("HelFitQtClass", "0\302\260                                90\302\260                              180\302\260                              270\302\260                             360\302\260", 0));
        label_14->setText(QApplication::translate("HelFitQtClass", "XZ                                                            YZ                                                        XY", 0));
        label_15->setText(QApplication::translate("HelFitQtClass", "Side Views:                                                         Top View:", 0));
        tabWidget->setTabText(tabWidget->indexOf(tab3Dview), QApplication::translate("HelFitQtClass", "3D View", 0));
        tabWidget->setTabText(tabWidget->indexOf(tabTextoutput), QApplication::translate("HelFitQtClass", "output", 0));
        groupBox_1->setTitle(QApplication::translate("HelFitQtClass", "Scattering Data Toolbar", 0));
        label->setText(QApplication::translate("HelFitQtClass", "Input Scattering File:", 0));
        btnLoadFile->setText(QApplication::translate("HelFitQtClass", "Load!", 0));
        label_2->setText(QApplication::translate("HelFitQtClass", "Plotting Range (datapoints):", 0));
        label_3->setText(QApplication::translate("HelFitQtClass", "Min:", 0));
        label_4->setText(QApplication::translate("HelFitQtClass", "Max:", 0));
        label_5->setText(QApplication::translate("HelFitQtClass", "Min:", 0));
        label_6->setText(QApplication::translate("HelFitQtClass", "Fitting Range (datapoints):", 0));
        label_7->setText(QApplication::translate("HelFitQtClass", "Max:", 0));
        label_8->setText(QApplication::translate("HelFitQtClass", "qMin:", 0));
        label_9->setText(QApplication::translate("HelFitQtClass", "qMax:", 0));
        label_10->setText(QApplication::translate("HelFitQtClass", "(1/nm)", 0));
        label_11->setText(QApplication::translate("HelFitQtClass", "(1/nm)", 0));
        groupBox_2->setTitle(QApplication::translate("HelFitQtClass", "DebyeStack parameters:", 0));
        label_16->setText(QApplication::translate("HelFitQtClass", "(1/nm)", 0));
        lineCalcQmin->setText(QApplication::translate("HelFitQtClass", "0", 0));
        label_17->setText(QApplication::translate("HelFitQtClass", "qMin:", 0));
        label_18->setText(QApplication::translate("HelFitQtClass", "(1/nm)", 0));
        lineCalcQmax->setText(QApplication::translate("HelFitQtClass", "5.0", 0));
        label_19->setText(QApplication::translate("HelFitQtClass", "qMax:", 0));
        label_20->setText(QApplication::translate("HelFitQtClass", "points:", 0));
        label_21->setText(QApplication::translate("HelFitQtClass", "No. stacks:", 0));
        label_22->setText(QApplication::translate("HelFitQtClass", "CPU cores:", 0));
        lineStackSpacing->setText(QApplication::translate("HelFitQtClass", "25.0", 0));
        label_23->setText(QApplication::translate("HelFitQtClass", "stack dist.:", 0));
        label_24->setText(QApplication::translate("HelFitQtClass", "(nm)", 0));
        btnFit->setText(QApplication::translate("HelFitQtClass", "FIT!", 0));
        menu_File->setTitle(QApplication::translate("HelFitQtClass", "File", 0));
        menuRunDebyeCalc->setTitle(QApplication::translate("HelFitQtClass", "DebyeStack current model", 0));
        menuOptions->setTitle(QApplication::translate("HelFitQtClass", "Options", 0));
    } // retranslateUi

};

namespace Ui {
    class HelFitQtClass: public Ui_HelFitQtClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_HELFITQT_H
