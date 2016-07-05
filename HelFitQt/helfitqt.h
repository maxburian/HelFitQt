#ifndef HELFITQT_H
#define HELFITQT_H

//#define BOOST_ALL_DYN_LINK

#include <QtWidgets/QMainWindow>
#include <QMouseEvent>
#include "qvector.h"
#include <string.h>
#include <vector>
#include <boost/thread.hpp>
#include <qapplication.h>
#include <boost/chrono.hpp>

#include "data_def.h"
#include "ui_helfitqt.h"
#include "dialog_modelvars.h"

class HelFitQt : public QMainWindow
{
	Q_OBJECT

public:
	HelFitQt(QWidget *parent = 0);
	~HelFitQt();
	void initialize_graphs();
	void init();

	//The main Fittingobject!!
	saxs::fittingobject_sp globalfittingobject{new saxs::fittingobject};
	std::vector<saxs::model_sp> result_tracer;

	//Fileloadparams
	QString scatterfilepath;
	QString scatterfilename;
	std::string savefilename;

	//RandomModelGeneration
	std::vector<saxs::coordinate_sp> generated_model;
	bool model_generated = false;

	//Scatteringdata objets
	std::vector<saxs::scatteringdata_sp> imported_data;
	QVector<double> x_expdata;
	QVector<double> y_expdata;
	QVector<double> e_expdata;
	bool data_loaded = false;
	QVector<double> x_fitdata;
	QVector<double> y_fitdata;
	QVector<double> e_fitdata;
	int min_datarange;
	int max_datarange;
	int num_datapoints;
	int min_fitrange;
	int max_fitrange;

	//Imported Model objects
	std::vector<saxs::scatteringdata_sp> imported_model_data;
	std::vector<saxs::coordinate_sp> imported_model;
	std::vector<saxs::coordinate_sp> plot_model_cyl;
	QVector<double> model_x;
	QVector<double> model_y;
	QVector<double> model_z;
	double phi_moodel_rotation;
	bool model_loaded = false;

	//Helperobjects
	static std::vector<double> sinc_lookup;
	boost::chrono::time_point<boost::chrono::system_clock> start;

	//Imported model calculation objects
	//ranges from x_imported_model_data min_calcrange to max_calcrange
	QVector<double> x_imported_model_data;
	QVector<double> y_imported_model_data;
	std::pair<double, double> linreg_results;
	double calc_qmin;
	double calc_qmax;
	int num_calcpoints;
	int num_stacks;
	int core_number;

	//Weight for curvefitting: 0=Is, 1=Iss, 2=ln(I)
	int curvefit_weight=1;

public slots:
	void writetologext(std::string text);
	void fittingthread_done(int);

private slots:
	//data load gui
	void clear_data();
	void on_btnLoadFile_clicked();
	void plot_data();
	void changedMinDataRange();
	void changedMaxDataRange();
	void changedMinFitRange();
	void changedMaxFitRange();
	void mousePressEvent(QMouseEvent *e);
	void replotData();

	//Fittingweight
	void actionWeight0();
	void actionWeight1();
	void actionWeight2();

	//model load gui
	void clear_model();
	void plot_model(std::vector<double>& boundaries);
	void loadModelFromData();
	void changedModelRotationPhi(int sliderpos);
	void replotModel();
	void plotCalcModel();
	void actionGenerateRandModel();

	//calc gui
	void calcDebyeStackCurrent();
	void refitDebyeStackCurrent();
	std::pair<double, double>refit_modeltoexp(int weighing);

	//Fitting
	void on_btnFit_clicked();
	void UpdateFittedModel(saxs::fittingobject_sp& fitobject);
	void startDebyeFitCurrentModel();
	void finishedDebyeFitCurrentModel();

	//Others
	void movedatatofittingobject();
	void writetolog(std::string text);
	

private:
	Ui::HelFitQtClass ui;
	bool eventFilter(QObject *obj, QEvent *ev);
	QFutureWatcher<void>  *m_future_watcher;
	


};
#endif // HELFITQT_H
