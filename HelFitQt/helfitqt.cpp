#include "helfitqt.h"
#include "data_def.h"
#include "data_utils.h"
#include "thread_utils.h"

#include "ui_helfitqt.h"
#include "qfiledialog.h"

#include "qfileinfo.h"
#include "qmessagebox.h"
#include <QMouseEvent>
#include "qevent.h"
#include "qobject.h"
#include <QDoubleValidator>
#include "qstring.h"
#include <algorithm>

#include <boost/scoped_ptr.hpp>
#include <boost/chrono.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>

//Initialize Sinc Lookup
std::vector<double> HelFitQt::sinc_lookup(100001, 1.0f);

HelFitQt::HelFitQt(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
}

HelFitQt::~HelFitQt()
{

}

void HelFitQt::init()
{
	initialize_graphs();
	//Initialize boxvalidator
	ui.lineStackSpacing->setValidator(new QDoubleValidator(0, 1000, 4, this));
	ui.lineCalcQmin->setValidator(new QDoubleValidator(0, 20, 4, this));
	ui.lineCalcQmax->setValidator(new QDoubleValidator(0, 20, 4, this));
	//Initialize Sinc Lookup
	HelFitQt::sinc_lookup[0] = 1;
	for (std::size_t i = 1; i < sinc_lookup.size(); ++i)
	{
		HelFitQt::sinc_lookup[i] = boost::math::sinc_pi((double)i / 1000);
	}
	//Get number of cores for debyecalc
	core_number = boost::thread::hardware_concurrency();
	ui.spbNrCores->setMaximum(core_number);
	
	//Install event filters
	ui.spbDataMax->installEventFilter(this);
	ui.spbDataMin->installEventFilter(this);
	ui.spbFitMax->installEventFilter(this);
	ui.spbFitMin->installEventFilter(this);

	writetolog("Start of Helfitqt....");
	writetolog("------------------------");
	writetolog("Welcome!");
	writetolog("------------------------");
	writetolog(" ");
}

// --------------------------------------------------------------------------------------------------
// Event Handling - Buttons/Menuitems
// --------------------------------------------------------------------------------------------------
void HelFitQt::on_btnLoadFile_clicked()
{
	if (data_loaded) clear_data();
	scatterfilepath = saxs::filedialog("Select data file:", "Data (*.qI);;All (*.*)");
	if (scatterfilepath != "none")
	{
		//scatterfilepath = QFileDialog::getOpenFileName(this, tr("Select File"), "/path/to/file/", tr("Data (*.qI);;All (*.*)"),, QFileDialog::DontUseNativeDialog);
		QFileInfo fi(scatterfilepath);
		scatterfilename = fi.fileName();
		ui.txtFilePath->setText(scatterfilename);

		//Import data from file inte scatteringdata object
		saxs::import_scatteringdata(scatterfilepath, imported_data);
		//convert object data into globals x_expdata y_expdata
		saxs::convertScatteringToVector(imported_data, x_expdata, y_expdata, e_expdata);
		//rescale data so no values are <1
		saxs::rescaleScatteringData(imported_data, y_expdata, e_expdata);

		//determine datarange globlas
		min_datarange = 0;
		min_fitrange = 0;
		num_datapoints = x_expdata.size();
		max_datarange = x_expdata.size() - 1;
		max_fitrange = x_expdata.size() - 1;
		num_calcpoints = max_fitrange + 1 - min_fitrange;

		///Altering userinterface for plotting
		ui.spbCalcNrPoints->setValue(num_calcpoints);
		ui.spbDataMin->setEnabled(true);
		ui.spbDataMax->setEnabled(true);
		ui.spbDataMin->setMinimum(1);
		ui.spbDataMin->setMaximum(num_datapoints);
		ui.spbDataMax->setMinimum(1);
		ui.spbDataMax->setMaximum(num_datapoints);
		ui.spbDataMin->setValue(1);
		ui.spbFitMin->setEnabled(true);
		ui.spbFitMax->setEnabled(true);
		ui.spbFitMin->setMinimum(1);
		ui.spbFitMin->setMaximum(num_datapoints);
		ui.spbFitMax->setMinimum(1);
		ui.spbFitMax->setMaximum(num_datapoints);
		ui.spbFitMin->setValue(1);
		ui.lineQmin->setText(QString::number(x_expdata[min_datarange]));
		ui.lineQmax->setText(QString::number(x_expdata[max_datarange]));
		ui.chkbLogLogPlot->setEnabled(true);
		//Altering userinterface for calc
		ui.lineCalcQmin->setEnabled(false);
		ui.lineCalcQmax->setEnabled(false);
		ui.spbCalcNrPoints->setEnabled(false);
		ui.lineCalcQmax->setText(QString::number(x_expdata[max_datarange]));
		ui.lineCalcQmin->setText(QString::number(x_expdata[min_datarange]));
		ui.spbCalcNrPoints->setValue(max_datarange - min_datarange + 1);

		//set global dataload
		data_loaded = true;

		//plot data
		plot_data();
		ui.tabWidget->setCurrentIndex(0);

		//Write info in ouput field
		writetolog(" ");
		std::string str = " ";
		writetolog(("Succesfully loaded data from: " + scatterfilepath.toStdString()));
		str = str+ boost::lexical_cast<std::string>(num_datapoints);
		str = str+ " points\t\t";
		str = str + "Qmin: " + ui.lineQmin->text().toStdString() + "\t\tQmax: " + ui.lineQmax->text().toStdString();
		writetolog(str);
	}
	else
	{
		ui.txtFilePath->setText("none selected");
	}
}

void HelFitQt::loadModelFromData()
{
	//Load data from file
	QString modelfilepath;
	modelfilepath = saxs::filedialog("Select model file:", "XYZ (*.xyz);;All (*.*)");
	if (modelfilepath != "none")
	{
		if (model_loaded) clear_model();
		saxs::import_xyz_model(modelfilepath, imported_model);
		//Copy original data into plotting qvectors
		saxs::convertCoordinateToVector(imported_model, model_x, model_y, model_z);
		std::vector<double> model_boundaries(3);
		model_boundaries = saxs::convertCoordinateToCylinder(imported_model, plot_model_cyl);
		plot_model(model_boundaries);
		phi_moodel_rotation = saxs::pi;
		ui.hsliderPhiRot->setValue(500);
		model_loaded = true;
		ui.menuDebyeCurrModel->setEnabled(true);
		ui.tabWidget->setCurrentIndex(1);

		//Copying data into globalfittingobject
		globalfittingobject->m_model.clear();
		globalfittingobject->m_model.resize(imported_model.size());
		globalfittingobject->m_model = imported_model;

		//Write info in ouput field
		writetolog(" ");
		std::string str = " ";
		writetolog(("Succesfully loaded model from: " + modelfilepath.toStdString()));
		str = str + "The model includes ";
		str = str + boost::lexical_cast<std::string>(imported_model.size());
		str = str + " DAs";
		writetolog(str);
	}
}

void HelFitQt::on_btnFit_clicked()
{
	//Check if model is loaded
	if (!model_loaded && !model_generated)
	{
		saxs::error_message("No model loaded!");
		return;
	}

	if (!data_loaded)
	{
		saxs::error_message("No data loaded!");
		return;
	}

	//Select target filename
	QString filters("Chi-files (*.chi);;PDB-files (*.pdb);;All files (*.*)");
	QString defaultFilter("Chi-files (*.chi)");
	/* Static method approach */
	QString filename = QFileDialog::getSaveFileName(0, "Save files to..", QDir::currentPath(),
		filters, &defaultFilter);
	if (filename == NULL) return;
	savefilename = filename.toStdString().substr(0, filename.size() - 4);

	//Make sure current ui values are in globals
	num_calcpoints = int(ui.spbCalcNrPoints->value());
	calc_qmin = double(ui.lineCalcQmin->text().toDouble());
	calc_qmax = double(ui.lineCalcQmax->text().toDouble());

	//get calculation parameters from gui
	//build fittingobject
	movedatatofittingobject();

	//Write info in ouput field
	writetolog(" ");
	writetolog("******* Debye Fit *******");
	std::string str = " ";
	str = "Calculating model data between q = " + ui.lineCalcQmin->text().toStdString() + " - " + ui.lineCalcQmax->text().toStdString();
	str = str + " using " + boost::lexical_cast<std::string>(num_calcpoints) + " points.";
	writetolog(str);
	writetolog("Parameters:");
	writetolog("Stacking distance\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_stack_spacing) + " nm");
	writetolog("Number of stacks\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_num_stacks));
	writetolog("Number of cores\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_num_cores));

	writetolog("START!");
	qApp->processEvents();
	//Get threads start time
	start = boost::chrono::system_clock::now();

	//TODO!!!!!!!!!!!!!!!!!!!!!!!!
	globalfittingobject->m_diameter = 25;
	//TODO!!!!!!!!!!!!!!!!!!!!!!!!

	//Disabling interface
	ui.menuBar->setEnabled(false);
	ui.groupBox_1->setEnabled(false);
	ui.groupBox_2->setEnabled(false);
	ui.btnFit->setEnabled(false);

	//Send to signaling function
	startDebyeFitCurrentModel();
}

void HelFitQt::calcDebyeStackCurrent()
{
	//Check if model is loaded
	if (!model_loaded && !model_generated)
	{
		saxs::error_message("No model loaded!");
		return;
	}

	//Make sure current ui values are in globals
	num_calcpoints = int(ui.spbCalcNrPoints->value());
	calc_qmin = double(ui.lineCalcQmin->text().toDouble());
	calc_qmax = double(ui.lineCalcQmax->text().toDouble());

	//get calculation parameters from gui
	//build fittingobject
	movedatatofittingobject();

	//Write info in ouput field
	writetolog(" ");
	writetolog("******* Debye Stack Calculation *******");
	std::string str = " ";
	str = "Calculating model data between q = " + ui.lineCalcQmin->text().toStdString() + " - " + ui.lineCalcQmax->text().toStdString();
	str = str +  " using " + boost::lexical_cast<std::string>(num_calcpoints) +" points.";
	writetolog(str);
	writetolog("Parameters:");
	writetolog("Stacking distance\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_stack_spacing)+" nm");
	writetolog("Number of stacks\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_num_stacks));
	writetolog("Number of cores\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_num_cores));

	writetolog("START!");
	qApp->processEvents();
	//Get threads start time
	auto start = boost::chrono::system_clock::now();
		
	saxs::calcDebyeStackCurrentModel(globalfittingobject);

	//Get threads end time
	auto end = boost::chrono::system_clock::now();

	//Calculate elapsed time
	boost::uint64_t elapsed_seconds =boost::chrono::duration_cast<boost::chrono::seconds>(end - start).count();
	writetolog("FINISHED!");
	writetolog("Duration\t\t = \t " + boost::lexical_cast<std::string>(elapsed_seconds) + " seconds.");

	//Clear previous calcdata
	//Old Datastructure!!!!
	imported_model_data.clear();
	x_imported_model_data.clear();
	y_imported_model_data.clear();

	//Move data to plotable qvector and normalize to I[0]
	for (int i = 0; i < globalfittingobject->m_data_q.size(); i++)
	{
		x_imported_model_data.push_back(globalfittingobject->m_data_q[i]);
		y_imported_model_data.push_back(globalfittingobject->m_fitted_I[i]);
		imported_model_data.push_back(saxs::scatteringdata_sp(new saxs::scatteringdata(
			globalfittingobject->m_data_q[i], globalfittingobject->m_fitted_I[i])));
	}

	plotCalcModel();
	
	//Enabling refit menu item
	ui.menuRefit->setEnabled(true);
}

void HelFitQt::refitDebyeStackCurrent()
{
	globalfittingobject->m_weighing = curvefit_weight;
	linreg_results = refit_modeltoexp(curvefit_weight);
	
	//Clear previous calcdata
	//Old Datastructure!!!!
	imported_model_data.clear();
	x_imported_model_data.clear();
	y_imported_model_data.clear();

	//Move data to plotable qvector and normalize to I[0]
	for (int i = 0; i < globalfittingobject->m_data_q.size(); i++)
	{
		x_imported_model_data.push_back(globalfittingobject->m_data_q[i]);
		y_imported_model_data.push_back(globalfittingobject->m_fitted_I[i]);
		imported_model_data.push_back(saxs::scatteringdata_sp(new saxs::scatteringdata(
			globalfittingobject->m_data_q[i], globalfittingobject->m_fitted_I[i])));
	}
	
	plotCalcModel();
}

void HelFitQt::actionWeight0()
{
	curvefit_weight = 0;
	ui.menuWeight0->setChecked(true);
	ui.menuWeight1->setChecked(false);
	ui.menuWeight2->setChecked(false);
}

void HelFitQt::actionWeight1()
{
	curvefit_weight = 1;
	ui.menuWeight0->setChecked(false);
	ui.menuWeight1->setChecked(true);
	ui.menuWeight2->setChecked(false);
}

void HelFitQt::actionWeight2()
{
	curvefit_weight = 2;
	ui.menuWeight0->setChecked(false);
	ui.menuWeight1->setChecked(false);
	ui.menuWeight2->setChecked(true);
}

void HelFitQt::actionGenerateRandModel()
{
	dialog_modelvars *mydialog_modelvars;
	mydialog_modelvars = new dialog_modelvars;
	mydialog_modelvars->exec();
	bool dialogOk = mydialog_modelvars->result();
	
	//Check if parameters are accepted
	if (!dialogOk)
	{
		saxs::error_message("No parameters specified..."); 
		//Destroy object
		delete mydialog_modelvars;
		return;
	}
	//Transfer values into HelFitQt class
	int nr_atoms = mydialog_modelvars->m_nr_atoms;
	double height = mydialog_modelvars->m_stackspacing;
	ui.lineStackSpacing->setText(QString::number(height));
	double diameter = mydialog_modelvars->m_diameter;
	
	//Reinit structures
	generated_model.clear();
	plot_model_cyl.clear();
	model_x.clear();
	model_y.clear();
	model_z.clear();

	//Generate ranodm cylindrical coords:
	saxs::gen_cyl_randomseed(plot_model_cyl, height, diameter, nr_atoms);
	//initialize plotting vectors and load z coords
	saxs::convertCoordinateToVector(plot_model_cyl, model_x, model_y, model_z);
	//Now convert variables into cartesion vectors
	saxs::convertCylinderToQVectors(plot_model_cyl, model_x, model_y, 0);
	//Now move everything into generated_model object
	saxs::convertQvectorsToCoordinate(generated_model, model_x, model_y, model_z);
	
	//Plotting
	std::vector<double> model_boundaries(3);
	model_boundaries = saxs::convertCoordinateToCylinder(generated_model, plot_model_cyl);
	plot_model(model_boundaries);
	phi_moodel_rotation = saxs::pi;
	ui.hsliderPhiRot->setValue(500);
	model_generated = true;
	ui.menuDebyeCurrModel->setEnabled(true);
	ui.tabWidget->setCurrentIndex(1);

	//Copying data into globalfittingobject
	globalfittingobject->m_model.clear();
	globalfittingobject->m_model.resize(generated_model.size());
	globalfittingobject->m_model = generated_model;

	//Write info in ouput field
	writetolog(" ");
	std::string str = " ";
	writetolog("Succesfully generated random model!");
	str = str + "The model includes ";
	str = str + boost::lexical_cast<std::string>(generated_model.size());
	str = str + " DAs";
	writetolog(str);

}

// --------------------------------------------------------------------------------------------------
// Event Handling - Sliders/Spinboxes/etc...
// --------------------------------------------------------------------------------------------------
void HelFitQt::mousePressEvent(QMouseEvent *e)
{
	if (e->button() == Qt::RightButton)
	{
		saxs::error_message("Test Right click!");
	}
}

bool HelFitQt::eventFilter(QObject *obj, QEvent *ev)
{
	//Filtering out keyboard input and opening popup
	if (obj == ui.spbDataMax || obj == ui.spbDataMin || obj == ui.spbFitMin || obj == ui.spbFitMax)
	{
		if (ev->type() == QEvent::KeyPress || ev->type() == QEvent::MouseButtonDblClick) {
			//QMouseEvent * mouseEvent = static_cast <QMouseEvent *> (ev);
			QSpinBox * spb = qobject_cast<QSpinBox * > (obj);
			int point;
			int min = spb->minimum();
			int max = spb->maximum();
			int curr = spb->value();
			QString text = " ";
			if (obj == ui.spbDataMax) text = "maximum plotting range:";
			if (obj == ui.spbDataMin) text = "minimum plotting range:";
			if (obj == ui.spbFitMax) text = "maximum fitting range:";
			if (obj == ui.spbFitMin) text = "minimum fitting range:";
			point = saxs::input_message(text, curr, min, max);
			spb->setValue(point);
			//saxs::error_message("Test Double click!" + QString::number(min));
			return true;
		}
		else return false;
	}
	else {
		// standard event processing
		return HelFitQt::eventFilter(obj, ev);
	}
}

void HelFitQt::changedMinDataRange()
{
	if (int(ui.spbDataMin->value() - 1)<max_datarange)
	{
		min_datarange = int(ui.spbDataMin->value()-1);
		ui.spbFitMin->setMinimum(min_datarange + 1);
	}
	else
	{
		ui.spbDataMin->setValue(min_datarange + 1);
	}
	ui.wdgtDataPlot->xAxis->setRange(x_expdata[min_datarange], x_expdata[max_datarange]);
	ui.wdgtDataPlot->replot();
}

void HelFitQt::changedMaxDataRange()
{
	if (int(ui.spbDataMax->value() - 1) > min_datarange)
	{
		max_datarange = int(ui.spbDataMax->value() - 1);
		ui.spbFitMax->setMaximum(max_datarange + 1);
	}
	else
	{
		ui.spbDataMax->setValue(max_datarange + 1);
	}
	ui.wdgtDataPlot->xAxis->setRange(x_expdata[min_datarange], x_expdata[max_datarange]);
	ui.wdgtDataPlot->replot();
}

void HelFitQt::changedMinFitRange()
{
	if (int(ui.spbFitMin->value() - 1) < max_fitrange)
	{
		min_fitrange = int(ui.spbFitMin->value() - 1);
		num_calcpoints = max_fitrange  + 1 - min_fitrange;
		replotData();
		ui.lineQmin->setText(QString::number(x_expdata[min_fitrange]));
		ui.lineCalcQmin->setText(QString::number(x_expdata[min_fitrange]));
		ui.spbCalcNrPoints->setValue(num_calcpoints);
	}
	else
	{
		ui.spbFitMin->setValue(min_fitrange + 1);
	}
}

void HelFitQt::changedMaxFitRange()
{
	if (int(ui.spbFitMax ->value()-1)>min_fitrange)
	{
		max_fitrange = int(ui.spbFitMax->value() - 1);
		num_calcpoints = max_fitrange + 1 - min_fitrange;
		replotData();
		ui.lineQmax->setText(QString::number(x_expdata[max_fitrange]));
		ui.lineCalcQmax->setText(QString::number(x_expdata[max_fitrange]));
		ui.spbCalcNrPoints->setValue(num_calcpoints);
	}
	else
	{
		ui.spbFitMax->setValue(max_fitrange + 2);
	}
}

void HelFitQt::changedModelRotationPhi(int sliderpos)
{
	phi_moodel_rotation = double(sliderpos) / 500 * saxs::pi;
	if (model_loaded || model_generated) replotModel();
}

// --------------------------------------------------------------------------------------------------
// Graphing tools
// --------------------------------------------------------------------------------------------------
void HelFitQt::initialize_graphs()
{
	//data plotting graphs
	ui.wdgtDataPlot->addGraph();
	// give the axes some labels:
	ui.wdgtDataPlot->xAxis->setLabel("q (1/nm)");
	ui.wdgtDataPlot->yAxis->setLabel("I(q) (-)");
	ui.wdgtDataPlot->yAxis->setRange(0, 1);
	ui.wdgtDataPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
	ui.wdgtDataPlot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
	ui.wdgtDataPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
	ui.wdgtDataPlot->yAxis->setScaleLogBase(10);
	ui.wdgtDataPlot->replot();

	//model plotting graphs
	ui.wdgtXZplot->addGraph();
	ui.wdgtXZplot->xAxis->setRange(-20, 20);
	ui.wdgtXZplot->xAxis->grid()->setVisible(false);
	ui.wdgtXZplot->xAxis2->setVisible(true);
	ui.wdgtXZplot->xAxis2->setTickLabels(false);
	ui.wdgtXZplot->yAxis->setRange(0, 50);
	ui.wdgtXZplot->yAxis2->setVisible(true);
	ui.wdgtXZplot->yAxis2->setTickLabels(false);
	ui.wdgtXZplot->yAxis->grid()->setVisible(false);
	ui.wdgtXZplot->replot();

	ui.wdgtYZplot->addGraph();
	ui.wdgtYZplot->xAxis->setRange(-20, 20);
	ui.wdgtYZplot->xAxis->grid()->setVisible(false);
	ui.wdgtYZplot->xAxis2->setVisible(true);
	ui.wdgtYZplot->xAxis2->setTickLabels(false);
	ui.wdgtYZplot->yAxis->setRange(0, 50);
	ui.wdgtYZplot->yAxis2->setVisible(true);
	ui.wdgtYZplot->yAxis2->setTickLabels(false);
	ui.wdgtYZplot->yAxis->grid()->setVisible(false);
	ui.wdgtYZplot->replot();

	ui.wdgtXYplot->addGraph();
	ui.wdgtXYplot->xAxis->setRange(-20, 20);
	ui.wdgtXYplot->xAxis->grid()->setVisible(false);
	ui.wdgtXYplot->xAxis2->setVisible(true);
	ui.wdgtXYplot->xAxis2->setTickLabels(false);
	ui.wdgtXYplot->yAxis->setRange(-20, 20);
	ui.wdgtXYplot->yAxis2->setVisible(true);
	ui.wdgtXYplot->yAxis2->setTickLabels(false);
	ui.wdgtXYplot->yAxis->grid()->setVisible(false);
	ui.wdgtXYplot->replot();
}

void HelFitQt::plot_data()
{
	x_fitdata.resize(num_calcpoints);
	y_fitdata.resize(num_calcpoints);
	e_fitdata.resize(num_calcpoints);
	x_fitdata = x_expdata;
	y_fitdata = y_expdata;
	e_fitdata = e_expdata;

	//add experimental data
	ui.wdgtDataPlot->addGraph();
	ui.wdgtDataPlot->graph(0)->setPen(QPen(Qt::blue));
	ui.wdgtDataPlot->graph(0)->setData(x_expdata, y_expdata);

	ui.wdgtDataPlot->addGraph();
	QPen RedPen;
	RedPen.setColor(Qt::red);
	RedPen.setWidthF(2);
	ui.wdgtDataPlot->graph(1)->setPen(RedPen);
	ui.wdgtDataPlot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::red, Qt::red, 1));
	ui.wdgtDataPlot->graph(1)->setErrorType(QCPGraph::etValue);
	ui.wdgtDataPlot->graph(1)->setErrorPen(QPen(Qt::red));
	ui.wdgtDataPlot->graph(1)->setDataValueError(x_fitdata, y_fitdata, e_fitdata);

	// give the axes some labels:
	ui.wdgtDataPlot->xAxis->setLabel("q (1/nm)");
	ui.wdgtDataPlot->yAxis->setLabel("I(q) (-)");


	// set axes ranges, so we see all data:
	double xmin = x_expdata[min_datarange];
	double xmax = x_expdata[max_datarange];
	ui.wdgtDataPlot->xAxis->setRange(xmin, xmax);
	if (ui.chkbLogLogPlot->checkState() == Qt::Checked)
	{
		ui.wdgtDataPlot->xAxis->setScaleLogBase(10);
		ui.wdgtDataPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
	}
	else
	{
		ui.wdgtDataPlot->xAxis->setScaleType(QCPAxis::stLinear);
	}

	double ymin = *std::min_element(y_expdata.constBegin() + min_datarange, y_expdata.constBegin() + max_datarange);
	double ymax = *std::max_element(y_expdata.constBegin() + min_datarange, y_expdata.constBegin() + max_datarange);
	ui.wdgtDataPlot->yAxis->setRange(ymin*0.5, 1.5*ymax);
	ui.wdgtDataPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
	ui.wdgtDataPlot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
	ui.wdgtDataPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
	ui.wdgtDataPlot->yAxis->setScaleLogBase(10);
	ui.wdgtDataPlot->replot();
}

void HelFitQt::replotData()
{
	x_fitdata.clear();
	y_fitdata.clear();
	e_fitdata.clear();
	x_fitdata.resize(num_calcpoints);
	y_fitdata.resize(num_calcpoints);
	e_fitdata.resize(num_calcpoints);

	for (int i = min_fitrange; i < max_fitrange+1; i++)
	{
		x_fitdata[i- min_fitrange] = x_expdata[i];
		y_fitdata[i- min_fitrange] = y_expdata[i];
		e_fitdata[i - min_fitrange] = e_expdata[i];
	}
	if (ui.chkbLogLogPlot->checkState() == Qt::Checked)
	{
		ui.wdgtDataPlot->xAxis->setScaleLogBase(10);
		ui.wdgtDataPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
	}
	else
	{
		ui.wdgtDataPlot->xAxis->setScaleType(QCPAxis::stLinear);
	}
	ui.wdgtDataPlot->graph(1)->setDataValueError(x_fitdata, y_fitdata, e_fitdata);
	ui.wdgtDataPlot->replot();
}

void HelFitQt::plot_model(std::vector<double>& boundaries)
{
	double rmax = boundaries[0] * 1.1;
	double zmax = boundaries[1];
	double zmin = boundaries[2];
	QPen RedPen;
	RedPen.setWidth(2);
	RedPen.setColor(Qt::green);
	ui.wdgtXZplot->addGraph();
	ui.wdgtXZplot->graph(1)->setPen(RedPen);
	ui.wdgtXZplot->graph(1)->setLineStyle(QCPGraph::lsNone);
	ui.wdgtXZplot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 1));
	ui.wdgtXZplot->graph(1)->setData(model_x, model_z);
	ui.wdgtXZplot->xAxis->setRange(-rmax, rmax);
	ui.wdgtXZplot->yAxis->setRange(zmin, zmax);
	ui.wdgtXZplot->replot();

	ui.wdgtYZplot->addGraph();
	ui.wdgtYZplot->graph(1)->setPen(RedPen);
	ui.wdgtYZplot->graph(1)->setLineStyle(QCPGraph::lsNone);
	ui.wdgtYZplot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 1));
	ui.wdgtYZplot->graph(1)->setData(model_y, model_z);
	ui.wdgtYZplot->xAxis->setRange(-rmax, rmax);
	ui.wdgtYZplot->yAxis->setRange(zmin, zmax);
	ui.wdgtYZplot->replot();

	ui.wdgtXYplot->addGraph();
	ui.wdgtXYplot->graph(1)->setPen(RedPen);
	ui.wdgtXYplot->graph(1)->setLineStyle(QCPGraph::lsNone);
	ui.wdgtXYplot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 1));
	ui.wdgtXYplot->graph(1)->setData(model_x, model_y);
	ui.wdgtXYplot->xAxis->setRange(-rmax, rmax);
	ui.wdgtXYplot->yAxis->setRange(-rmax, rmax);
	ui.wdgtXYplot->replot();
}

void HelFitQt::replotModel() 
{
	saxs::convertCylinderToQVectors(plot_model_cyl,model_x, model_y,phi_moodel_rotation);
	ui.wdgtXZplot->graph(1)->setData(model_x, model_z);
	ui.wdgtXZplot->replot();
	ui.wdgtYZplot->graph(1)->setData(model_y, model_z);
	ui.wdgtYZplot->replot();
	ui.wdgtXYplot->graph(1)->setData(model_x, model_y);
	ui.wdgtXYplot->replot();
}

void HelFitQt::plotCalcModel()
{
	if (data_loaded)
	{
		double chi = globalfittingobject->m_chi;
		std::string str;
		switch (curvefit_weight) {
			case (0) : str = "Chisquare (I(q)*q)\t = \t"; break;
			case (1) : str = "Chisquare (I(q)*q^2)\t = \t"; break;
			case (2) : str = "Chisquare (ln(I(q)))\t = \t"; break;
		}
		writetolog(str + boost::lexical_cast<std::string>(chi));
	}
	int k = 0;
	ui.wdgtDataPlot->addGraph();
	if (data_loaded) { k = 2; }
	QPen BlackPen;
	BlackPen.setColor(Qt::black);
	BlackPen.setWidthF(2);
	ui.wdgtDataPlot->graph(k)->setPen(BlackPen);
	ui.wdgtDataPlot->graph(k)->setData(x_imported_model_data, y_imported_model_data);
	if (!data_loaded)
	{
		double xmin = x_imported_model_data[0];
		double xmax = x_imported_model_data[x_imported_model_data.size()-1];
		ui.wdgtDataPlot->xAxis->setRange(xmin, xmax);
		double ymin = *std::min_element(y_imported_model_data.constBegin(), y_imported_model_data.constEnd()-1);
		double ymax = *std::max_element(y_imported_model_data.constBegin(), y_imported_model_data.constEnd()-1);
		ui.wdgtDataPlot->yAxis->setRange(0.5*ymin, 1.5*ymax);
		ui.wdgtDataPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
		ui.wdgtDataPlot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
		ui.wdgtDataPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
		ui.wdgtDataPlot->yAxis->setScaleLogBase(10);
	}

	ui.wdgtDataPlot->replot();
	ui.tabWidget->setCurrentIndex(0);
}

void HelFitQt::UpdateFittedModel(saxs::fittingobject_sp& fitobject)
{

	clear_model();
	imported_model = fitobject->m_model;
	//Copy original data into plotting qvectors
	saxs::convertCoordinateToVector(imported_model, model_x, model_y, model_z);
	std::vector<double> model_boundaries(3);
	model_boundaries = saxs::convertCoordinateToCylinder(imported_model, plot_model_cyl);
	plot_model(model_boundaries);
	phi_moodel_rotation = saxs::pi;
	ui.hsliderPhiRot->setValue(500);
	model_loaded = true;
	ui.menuDebyeCurrModel->setEnabled(true);
	ui.tabWidget->setCurrentIndex(1);
}

//************************************************************************************************************************
//Member functions of HelFitQt
//************************************************************************************************************************

//------------------------------------------------------------------------------
//	Only used to refit usinng menu - OBSOLETE!!!!
//------------------------------------------------------------------------------
std::pair<double, double> HelFitQt::refit_modeltoexp(int weighing)
{
	//Curve weighting: 0...none,1...q^2, 2...ln
	std::vector<double> exp;
	std::vector<double> model;
	std::vector<double> weight;
	std::pair<double, double> fitparams;

	double min_m = 10000000000000;
	double min_e = 10000000000000;

	//Loading data in standard vectors and finding minimum value
	for (int i = 0; i < num_calcpoints; i++)
	{
		exp.push_back(y_fitdata[i]);
		if (exp[i] < min_e) min_e = exp[i];
		model.push_back(y_imported_model_data[i]);
		if (model[i] < min_m) min_m = model[i];
	}
	weight.resize(exp.size());
	//Adding minimum value to make sure only positive values exist are used
	if (min_e < 0) min_e = min_e*1.01;
	else min_e = min_e*0.99;
	min_m = 0.99*min_m;
	for (int i = 0; i < num_calcpoints; i++)
	{
		exp[i] = exp[i] - min_e;
		model[i] = model[i] - min_m;
	}

	//Apply weighting conditions
	for (int i = 0; i < exp.size(); i++)
	{
		if (weighing == 0)
		{
			weight[i] = x_imported_model_data[i];
		}
		if (weighing == 1)
		{
			weight[i] = std::pow(x_imported_model_data[i], 2);
		}
		if (weighing == 2)
		{
			exp[i] = std::log(exp[i]);
			model[i] = std::log(model[i]);
			weight[i] = 1;
		}
	}

	//GetFitparams
	fitparams = saxs::lindatafit(exp, model, weight);

	if (weighing < 2)
	{
		for (int i = 0; i < x_imported_model_data.size(); i++)
		{
			y_imported_model_data[i] = fitparams.first + fitparams.second * y_imported_model_data[i];
		}
	}
	else
	{
		for (int i = 0; i < x_imported_model_data.size(); i++)
		{
			y_imported_model_data[i] = std::exp(fitparams.second*std::log(y_imported_model_data[i]) + fitparams.first);
		}
	}
	return fitparams;
}

//------------------------------------------------------------------------------
//	Clear temporary data for reload
//------------------------------------------------------------------------------
void HelFitQt::clear_data()
{
	imported_data.clear();
	x_expdata.clear();
	y_expdata.clear();
	x_fitdata.clear();
	y_fitdata.clear();
}

void HelFitQt::clear_model()
{
	imported_model_data.clear();
	imported_model.clear();
	plot_model_cyl.clear();
	model_x.clear();
	model_y.clear();
	model_z.clear();
	phi_moodel_rotation = 0;
}

//------------------------------------------------------------------------------
//	Write date : + text into logfield
//------------------------------------------------------------------------------
void HelFitQt::writetolog(std::string text)
{
	std::string time = saxs::now_str();
	QString qstr = QString::fromStdString(time + ": " + text);
	ui.textOutput->append(qstr);
	qApp->processEvents();
}

void HelFitQt::writetologext(std::string text)
{
	writetolog(text);
	qApp->processEvents();
}

//------------------------------------------------------------------------------
//	Create globalfittingobject
//------------------------------------------------------------------------------
void HelFitQt::movedatatofittingobject()
{

	//First move data to global object
	globalfittingobject->m_data_q.clear();
	globalfittingobject->m_data_I.clear();
	globalfittingobject->m_data_e.clear();
	globalfittingobject->m_model_I.clear();
	globalfittingobject->m_fitted_I.clear();
	if (data_loaded)
	{
		globalfittingobject->m_data_q.resize(x_fitdata.size());
		globalfittingobject->m_data_q = x_fitdata.toStdVector();
		globalfittingobject->m_data_I.resize(x_fitdata.size());
		globalfittingobject->m_data_I = y_fitdata.toStdVector();
		globalfittingobject->m_data_e.resize(x_fitdata.size());
		globalfittingobject->m_data_e = e_fitdata.toStdVector();
		globalfittingobject->m_model_I.resize(x_fitdata.size());
		globalfittingobject->m_fitted_I.resize(x_fitdata.size());
	}
	else
	{
		//Resize data intensity to 0 (used as boolean later on)
		globalfittingobject->m_data_I.resize(0);
		//Generate linsapce data
		double step = double(calc_qmax - calc_qmin) / double(num_calcpoints - 1);
		double temp = calc_qmin;
		for (int i = 0; i < num_calcpoints; i++)
		{
			globalfittingobject->m_data_q.push_back(temp);
			temp += step;
		}
		globalfittingobject->m_model_I.resize(num_calcpoints);
		globalfittingobject->m_fitted_I.resize(num_calcpoints);

	}
	std::fill(globalfittingobject->m_model_I.begin(), globalfittingobject->m_model_I.end(), 0);
	std::fill(globalfittingobject->m_fitted_I.begin(), globalfittingobject->m_fitted_I.end(), 0);

	//Now move params to global object
	globalfittingobject->m_num_stacks = int(ui.spbCalcNrStacks->value());
	globalfittingobject->m_stack_spacing = double(ui.lineStackSpacing->text().toDouble());
	globalfittingobject->m_num_cores = int(ui.spbNrCores->value());
	globalfittingobject->m_weighing = curvefit_weight;

}

//------------------------------------------------------------------------------
//	FITTING PROCEDURE!!! Node function using qthread worker
//------------------------------------------------------------------------------
void HelFitQt::startDebyeFitCurrentModel()
{
	//TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int n_runs = 100;
	double starttemp = 0.7;
	double deltatemp = 0.99;
	double endtemp = 0.0;
	//TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!

	double chi = globalfittingobject->m_chi;
	//Helper string
	std::string str;

	//Initialize thread parameters
	saxs::calcDebyeStackCurrentModel(globalfittingobject);
	switch (globalfittingobject->m_weighing) {
	case (0) : str = "Chisquare (I(q)*q)\t = \t"; break;
	case (1) : str = "Chisquare (I(q)*q^2)\t = \t"; break;
	case (2) : str = "Chisquare (ln(I(q)))\t = \t"; break;
	}
	writetologext("Initial " + str + boost::lexical_cast<std::string>(chi));

	//Preparing tracer object:
	result_tracer.clear();
	for (unsigned int i = 0; i < globalfittingobject->m_num_cores; i++)
	{
		result_tracer.push_back(saxs::model_sp(new saxs::model));
		result_tracer[i]->m_chi = globalfittingobject->m_chi;
		result_tracer[i]->m_chi_tracer.push_back(globalfittingobject->m_chi);
		result_tracer[i]->m_thread_status = false;
		for (int j = 0; j < globalfittingobject->m_model.size(); j++)
		{
			result_tracer[i]->m_coordinate.push_back(saxs::coordinate_sp(new saxs::coordinate(
				globalfittingobject->m_model[j]->m_x,
				globalfittingobject->m_model[j]->m_y,
				globalfittingobject->m_model[j]->m_z)));
		}
		for (int j = 0; j < globalfittingobject->m_model_I.size(); j++)
		{
			result_tracer[i]->m_model_I.push_back(0.0f);
			result_tracer[i]->m_fitted_I.push_back(0.0f);
		}
	}

	//Now everything is ready for calculation
	//Instantiate threadvector
	std::vector<DebyeFitThread*> thread_vector;

	qRegisterMetaType<std::string>();
	for (int i = 0; i < globalfittingobject->m_num_cores; i++)
	{
		thread_vector.push_back(new DebyeFitThread);
		thread_vector[i]->reftofittingobject = globalfittingobject;
		thread_vector[i]->result_tracer = result_tracer;
		thread_vector[i]->n_runs = n_runs;
		thread_vector[i]->starttemp = starttemp;
		thread_vector[i]->endtemp = endtemp;
		thread_vector[i]->deltatemp = deltatemp;
		thread_vector[i]->coreid = i;
		QObject::connect(thread_vector[i], SIGNAL(sig_writelogext(std::string)), this, SLOT(writetologext(std::string)));
		QObject::connect(thread_vector[i], SIGNAL(sig_fittingthread_done(int)), this, SLOT(fittingthread_done(int)));
	}

	//Start Threads
	for (int i = 0; i < globalfittingobject->m_num_cores; i++)
	{
		thread_vector[i]->start();
	}
}

void HelFitQt::fittingthread_done(int finished_thread)
{
	result_tracer[finished_thread]->m_thread_status = true;
	bool helper = true;
	for (int i = 0; i < result_tracer.size(); i++)
	{
		if (!result_tracer[i]->m_thread_status) helper = false;
	}
	if (helper) finishedDebyeFitCurrentModel();
}

void HelFitQt::finishedDebyeFitCurrentModel()
{
	//Save results to files
	std::string tempstring;
	for (unsigned int i = 0; i < globalfittingobject->m_num_cores; i++)
	{
		tempstring = savefilename + "_" + std::to_string(i);
		saxs::writemodeltopdb(tempstring, result_tracer[i]->m_coordinate);
		saxs::writedatatochi(tempstring, globalfittingobject->m_data_q, globalfittingobject->m_data_I, result_tracer[i]->m_fitted_I);
	}
	saxs::writechisquretofile(savefilename, result_tracer);

	//Check which had the best result
	double chimin = globalfittingobject->m_chi;
	int chipntr = 0;

	//Check which solution was best
	for (unsigned int i = 0; i < globalfittingobject->m_num_cores; i++)
	{
		if (result_tracer[i]->m_chi < chimin)
		{
			chimin = result_tracer[i]->m_chi;
			chipntr = i;
		}
	}

	//If there is a better solution write best one into globalfittingobject
	if (globalfittingobject->m_chi != chimin)
	{
		globalfittingobject->m_chi = result_tracer[chipntr]->m_chi;
		globalfittingobject->m_model_I = result_tracer[chipntr]->m_model_I;
		for (int j = 0; j < globalfittingobject->m_model.size(); j++)
		{
			globalfittingobject->m_model[j]->m_x = result_tracer[chipntr]->m_coordinate[j]->m_x;
			globalfittingobject->m_model[j]->m_y = result_tracer[chipntr]->m_coordinate[j]->m_y;
			globalfittingobject->m_model[j]->m_z = result_tracer[chipntr]->m_coordinate[j]->m_z;
		}
	}

	//Output to extfield
	std::stringstream ss;
	ss << "The best result has core  " << chipntr << " with chi=";
	ss << std::scientific;
	ss << chimin;
	writetolog(ss.str());


	writetolog("Results saved in " + savefilename + "_*.chi");

	//Get threads end time
	boost::chrono::time_point<boost::chrono::system_clock> end = boost::chrono::system_clock::now();

	//Calculate elapsed time
	boost::uint64_t elapsed_seconds = boost::chrono::duration_cast<boost::chrono::seconds>(end - start).count();
	writetolog("FINISHED!");
	writetolog("Duration\t\t = \t " + boost::lexical_cast<std::string>(elapsed_seconds) + " seconds.");

	//Clear previous calcdata
	//Old Datastructure!!!!
	imported_model_data.clear();
	x_imported_model_data.clear();
	y_imported_model_data.clear();

	//Move data to plotable qvector and normalize to I[0]
	for (int i = 0; i < globalfittingobject->m_data_q.size(); i++)
	{
		x_imported_model_data.push_back(globalfittingobject->m_data_q[i]);
		y_imported_model_data.push_back(globalfittingobject->m_fitted_I[i]);
		imported_model_data.push_back(saxs::scatteringdata_sp(new saxs::scatteringdata(
			globalfittingobject->m_data_q[i], globalfittingobject->m_fitted_I[i])));
	}


	//GUI changes
	plotCalcModel();
	UpdateFittedModel(globalfittingobject);
	//Enabling refit menu item
	ui.menuRefit->setEnabled(true);

	//Re-enabking interface
	ui.menuBar->setEnabled(true);
	ui.groupBox_1->setEnabled(true);
	ui.groupBox_2->setEnabled(true);
	ui.btnFit->setEnabled(true);


}