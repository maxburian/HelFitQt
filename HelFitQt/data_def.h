#ifndef DATA_DEF_H
#define DATA_DEF_H

#include <boost\thread.hpp>
#include <boost\shared_ptr.hpp>
#include <QtCore>
#include <qthread.h>

Q_DECLARE_METATYPE(std::string);

namespace saxs
{
	const double pi = 3.14159265359;
//------------------------------------------------------------------------------
//	Coordinate struct definition and typedef
//------------------------------------------------------------------------------
struct coordinate
{
	coordinate(double x = 0, double y = 0, double z = 0) :
		m_x(x), m_y(y), m_z(z) { }

	double m_x;
	double m_y;
	double m_z;
};

typedef boost::shared_ptr<coordinate> coordinate_sp;

struct model
{
	std::vector<coordinate_sp> m_coordinate;
	double m_chi;
	std::vector<double> m_model_I;
	std::vector<double> m_fitted_I;
	std::vector<double> m_chi_tracer;
	bool m_thread_status; //false - running / true - done
};

typedef boost::shared_ptr<model> model_sp;

//------------------------------------------------------------------------------
//	Scatteringdata struct definition and typedef
//------------------------------------------------------------------------------
struct scatteringdata
{
	scatteringdata(double q = 0, double I = 0, double e=0) :
		m_q(q), m_I(I), m_e(e) { }

	double m_q;
	double m_I;
	double m_e;
};

typedef boost::shared_ptr<scatteringdata> scatteringdata_sp;


//------------------------------------------------------------------------------
//	Fitting object struct definition and typedef
//------------------------------------------------------------------------------
struct fittingobject
{
	/*fittingobject(std::vector<saxs::coordinate_sp> model, std::vector<double> data_q, \
		std::vector<double> data_I, std::vector<double> model_I, double chi=0) : 
		m_model(model), m_data_q (data_q), m_data_I(data_I), m_model_I(model_I), m_chi(chi) {}
		*/
	std::vector<saxs::coordinate_sp> m_model;
	std::vector<double> m_data_q;
	std::vector<double> m_data_I;
	std::vector<double> m_data_e;
	std::vector<double> m_model_I;
	std::vector<double> m_fitted_I;
	double m_chi;
	std::pair<double, double> m_linreg_results;
	double m_norm;
	double m_norm_offset;
	//Parameters
	int m_num_stacks;
	double m_stack_spacing;
	double m_diameter;
	//SA Parameters
	int m_num_cores;
	int m_num_runspertemp;
	int m_num_runs;
	int m_weighing;
};
typedef boost::shared_ptr<fittingobject> fittingobject_sp;

}	//End of namespace



//QT THREAD!!!!
class DebyeFitThread : public QThread
{
	Q_OBJECT

	public:
		saxs::fittingobject_sp reftofittingobject;
		std::vector<saxs::model_sp> result_tracer;
		int coreid;
		int n_runs;
		double starttemp;
		double endtemp;
		double deltatemp;

		void run();

	signals:
		void sig_writelogext(std::string);
		void sig_fittingthread_done(int);

};

#endif /*!DATA_DEF_H*/