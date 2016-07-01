#ifndef DATA_UTILS_H
#define DATA_UTILS_H

#include <boost/make_shared.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "data_def.h"
#include "qtextstream.h"
#include "qstring.h"
#include "qvector.h"
#include "qevent.h"
#include "qinputdialog.h"
#include "qmessagebox.h"
#include "qfiledialog.h"
#include "helfitqt.h"
#include "thread_utils.h"

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <time.h>       /* time */


namespace saxs
{

	//------------------------------------------------------------------------------
	//	Error Message box
	//------------------------------------------------------------------------------
	void error_message(QString errormessage)
	{
		//Predefine Messagebox for error warnings
		QMessageBox msgBox;
		msgBox.setText("Error");
		msgBox.setStandardButtons(QMessageBox::Ok);
		msgBox.setIcon(QMessageBox::Critical);

		msgBox.setInformativeText(errormessage);
		int ret = msgBox.exec();

	}

	//------------------------------------------------------------------------------
	//	Integer inputbox
	//------------------------------------------------------------------------------
	int input_message(QString text,int curr, int min, int max)
	{
		//Predefine Messagebox for error warnings
		bool ok;
		int input = -1;
		input = QInputDialog::getInt(NULL, "Specify..",text, curr, min, max, 1, &ok);

		if (ok && input > 0) return input;
		else return curr;
	}

	//------------------------------------------------------------------------------
	//	Extract values from scatteringdataobject into globals 
	//------------------------------------------------------------------------------
	QString filedialog(QString dialogTitle, QString FileFilter)
	{
		QFileDialog dialog;
		dialog.setNameFilter(FileFilter);
		dialog.setFileMode(QFileDialog::ExistingFile);
		dialog.setWindowTitle(dialogTitle);
		dialog.setViewMode(QFileDialog::Detail);
		dialog.setOption(QFileDialog::DontUseNativeDialog, true);
		QStringList fileNames;
		if (dialog.exec())
		{
			fileNames = dialog.selectedFiles();
		}
		else
		{
			fileNames.append("none");
		}
		return fileNames[0];
	}

	//------------------------------------------------------------------------------
	//	Import scattering data from file
	//------------------------------------------------------------------------------
	void import_scatteringdata(const QString file_name, std::vector<scatteringdata_sp>& data)
	{
		//transform qstring into string
		std::string file_name_string = file_name.toStdString();

		//Predefine Messagebox for error warnings
		QMessageBox msgBox;
		msgBox.setText("Error");
		msgBox.setStandardButtons(QMessageBox::Ok);
		msgBox.setIcon(QMessageBox::Critical);

		//Open file 
		std::ifstream infile(file_name_string.c_str(), std::ifstream::in);

		//Check if correctly open
		if (!infile.is_open())
		{
			std::stringstream error_stream;
			msgBox.setInformativeText("Could not open file " + file_name + "!");
			int ret = msgBox.exec();
			return;
		}
		//Temporary variable used during import
		double q = 0;
		double I = 0;

		//Temporary line used during import
		std::string line;

		//Read and parse file line by line
		int counter = 0;
		while (getline(infile, line))
		{
			counter += 1;

			std::istringstream iss(line);
			iss >> q >> I;

			//Check for errors
			if (iss.fail())
			{
				std::stringstream error_stream;
				msgBox.setInformativeText("Error in " + file_name + " at line " + QString::number(counter));
				int ret = msgBox.exec();
				return;
			}

			//Store coordinates
			data.push_back(scatteringdata_sp(new scatteringdata(q, I)));
		}

		infile.close();
	}

	//------------------------------------------------------------------------------
	//	Import 3d model data from file
	//------------------------------------------------------------------------------
	void import_xyz_model(const QString file_name, std::vector<coordinate_sp>& data)
	{
		//transform qstring into string
		std::string file_name_string = file_name.toStdString();

		//Predefine Messagebox for error warnings
		QMessageBox msgBox;
		msgBox.setText("Error");
		msgBox.setStandardButtons(QMessageBox::Ok);
		msgBox.setIcon(QMessageBox::Critical);

		//Open file 
		std::ifstream infile(file_name_string.c_str(), std::ifstream::in);

		//Check if correctly open
		if (!infile.is_open())
		{
			std::stringstream error_stream;
			msgBox.setInformativeText("Could not open file " + file_name + "!");
			int ret = msgBox.exec();
			return;
		}
		//Temporary variable used during import
		double x = 0;
		double y = 0;
		double z = 0;

		//Temporary line used during import
		std::string line;

		//Read and parse file line by line
		int counter = 0;
		while (getline(infile, line))
		{
			counter += 1;

			std::istringstream iss(line);
			iss >> x >> y >> z;

			//Check for errors
			if (iss.fail())
			{
				std::stringstream error_stream;
				msgBox.setInformativeText("Error in " + file_name + " at line " + QString::number(counter));
				int ret = msgBox.exec();
				return;
			}

			//Store coordinates
			data.push_back(coordinate_sp(new coordinate(x, y, z)));
		}

		infile.close();
	}

	//Resorting points by z coordinate

	//Sort helper function
	std::vector<size_t> ordered(std::vector<double> const& values) 
	{
		std::vector<size_t> indices(values.size());
		std::iota(begin(indices), end(indices), static_cast<size_t>(0));

		std::sort(
			begin(indices), end(indices),
			[&](size_t a, size_t b) { return values[a] < values[b]; }
		);
		return indices;
	}

	void resortbyzcoord(fittingobject_sp& reftofittingobject)
	{
		std::vector<size_t> pointers;
		std::vector<double> z_helper,y_helper,x_helper;
		pointers.resize(reftofittingobject->m_model.size());
		z_helper.resize(reftofittingobject->m_model.size());
		y_helper.resize(reftofittingobject->m_model.size());
		x_helper.resize(reftofittingobject->m_model.size());
		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
		{
			z_helper[i] = reftofittingobject->m_model[i]->m_z;
			y_helper[i] = reftofittingobject->m_model[i]->m_y;
			x_helper[i] = reftofittingobject->m_model[i]->m_x;
		}
		pointers = ordered(z_helper);
		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
		{
			reftofittingobject->m_model[i]->m_x = x_helper[pointers[i]];
			reftofittingobject->m_model[i]->m_y = y_helper[pointers[i]];
			reftofittingobject->m_model[i]->m_z = z_helper[pointers[i]];
		}

	}

	//------------------------------------------------------------------------------
	//	Extract values from scatteringdataobject into globals 
	//------------------------------------------------------------------------------
	void convertScatteringToVector(std::vector<saxs::scatteringdata_sp> scatterobject, QVector<double>& x, QVector<double>& y)
	{
		int n = scatterobject.size();
		int i = 0;
		x.resize(n);
		y.resize(n);

		for (int i = 0; i<n; ++i)
		{
			x[i] = scatterobject[i]->m_q;
			y[i] = scatterobject[i]->m_I;
		}
	}

	//------------------------------------------------------------------------------
	//	Extract values from Coordinateobject into globals 
	//------------------------------------------------------------------------------
	void convertCoordinateToVector(std::vector<saxs::coordinate_sp> coordinate, QVector<double>& x, QVector<double>& y, QVector<double>& z)
	{
		int n = coordinate.size();
		int i = 0;
		x.resize(n);
		y.resize(n);
		z.resize(n);

		for (int i = 0; i<n; ++i)
		{
			x[i] = coordinate[i]->m_x;
			y[i] = coordinate[i]->m_y;
			z[i] = coordinate[i]->m_z;
		}
	}
	
	//------------------------------------------------------------------------------
	//	Convert QVector values into Coordinateobject
	//------------------------------------------------------------------------------
	void convertQvectorsToCoordinate(std::vector<saxs::coordinate_sp>& loc_coordinate, QVector<double>& x, QVector<double>& y, QVector<double>& z)
	{
		int n = x.size();
		int i = 0;
		loc_coordinate.clear();
		for (int i = 0; i<n; ++i)
		{
			loc_coordinate.push_back(coordinate_sp(new coordinate(x[i], y[i], z[i])));
		}
	}

	//------------------------------------------------------------------------------
	//	Convert values from Coordinateobject into cylindrical cartesian system 
	//------------------------------------------------------------------------------
	std::vector<double> convertCoordinateToCylinder(const std::vector<saxs::coordinate_sp>& cartesian, std::vector<saxs::coordinate_sp>& cylinder)
	{
		cylinder.clear();
		int n = cartesian.size();
		int i = 0;
		std::vector<double> boundaries(3);
		double r, phi, z;
		double rmax=0;
		double zmax=0;
		double zmin=100;
		for (int i = 0; i<n; ++i)
		{
			r = std::pow(std::pow(cartesian[i]->m_x, 2)+ std::pow(cartesian[i]->m_y, 2), 0.5);
			phi= std::atan2(cartesian[i]->m_y, cartesian[i]->m_x);
			z= cartesian[i]->m_z;
			if (rmax < r) {rmax = r;}
			if (zmax < z) { zmax = z; }
			if (zmin > z) { zmin = z; }
			cylinder.push_back(coordinate_sp(new coordinate(r, phi, z)));
		}
		boundaries[0] = rmax;
		boundaries[1] = zmax;
		boundaries[2] = zmin;
		return boundaries;
	}

	//------------------------------------------------------------------------------
	//	Convert values from cylindrical object into plotting QVectors 
	//------------------------------------------------------------------------------
	void convertCylinderToQVectors(const std::vector<saxs::coordinate_sp>& cylindrical,
													QVector<double>& x, QVector<double>& y,	const double phi)
	{
		int n = cylindrical.size();
		int i = 0;
		for (int i = 0; i<n; ++i)
		{
			x[i] = cylindrical[i]->m_x * std::cos(cylindrical[i]->m_y + phi);
			y[i] = cylindrical[i]->m_x * std::sin(cylindrical[i]->m_y + phi);
		}
	}
	
	//------------------------------------------------------------------------------
	//	Transfer cartesian into cylindrical coordinates
	//------------------------------------------------------------------------------
	void dupl_coord_by_ref(std::vector<saxs::coordinate_sp>& cart, std::vector<saxs::coordinate_sp>& cyl)
	{
		int n = cart.size();
		int i = 0;
	}

	//------------------------------------------------------------------------------
	//	Returns current time as string
	//------------------------------------------------------------------------------
	std::string now_str()
	{
		// Get current time from the clock, using microseconds resolution
		const boost::posix_time::ptime now =
			boost::posix_time::microsec_clock::local_time();

		// Get the time offset in current day
		const boost::posix_time::time_duration td = now.time_of_day();

		//
		// Extract hours, minutes, seconds and milliseconds.
		//
		// Since there is no direct accessor ".milliseconds()",
		// milliseconds are computed _by difference_ between total milliseconds
		// (for which there is an accessor), and the hours/minutes/seconds
		// values previously fetched.
		//
		const long hours = td.hours();
		const long minutes = td.minutes();
		const long seconds = td.seconds();
		const long milliseconds = td.total_milliseconds() -
			((hours * 3600 + minutes * 60 + seconds) * 1000);

		//
		// Format like this:
		//
		//      hh:mm:ss.SSS
		//
		// e.g. 02:15:40:321
		//
		//      ^          ^
		//      |          |
		//      123456789*12
		//      ---------10-     --> 12 chars + \0 --> 13 chars should suffice
		//  
		// 
		char buf[40];
		sprintf(buf, "%02ld:%02ld:%02ld.%03ld",
			hours, minutes, seconds, milliseconds);

		return buf;
	}

	//------------------------------------------------------------------------------
	//	Base function to save model as *.pdb file and scattering data as *.chi
	//------------------------------------------------------------------------------
	//Helper functions
	std::string prd(double x, int decDigits, const int width)
	{
		std::stringstream ss;
		ss << std::right;
		ss.fill(' ');        // fill space around displayed #
		ss.width(width);     // set  width around displayed #
		ss.precision(decDigits); // set # places after decimal
		ss << std::fixed;
		ss << x;
		return ss.str();
	}
	std::string pri(int x, const int width)
	{
		std::stringstream ss;
		ss << std::right;
		ss.fill(' ');        // fill space around displayed #
		ss.width(width);     // set  width around displayed #
		ss.precision(0);	// set # places after decimal
		ss << std::fixed;
		ss << x;
		return ss.str();
	}

	void writemodeltopdb(std::string file_name, std::vector<coordinate_sp>& m_model)
	{
		//add fileextension
		file_name += ".pdb";

		//Open file 
		std::ofstream outfile(file_name.c_str(), std::ifstream::out);

		//Check if correctly open
		if (!outfile.is_open())
		{
			std::stringstream error_stream;
			error_message(QString::fromStdString("Cannot open " + file_name));
			return;
		}

		//Store results to file
		for (int i = 0; i<m_model.size(); ++i)
		{
			outfile << "ATOM  " << pri((i + 1), 5) << "  CA  ASP" << pri((1 + int(i / 10)), 5) << "    "
				<< prd(m_model[i]->m_x, 3, 8)
				<< prd(m_model[i]->m_y, 3, 8)
				<< prd(m_model[i]->m_z, 3, 8)
				<< prd(1.0, 2, 6)
				<< prd(20.0, 2, 6)
				<< "           C"
				<< std::endl;
		}
		//Close file
		outfile.close();
	}
	void writedatatochi(std::string file_name, const std::vector<double>& m_data_q, const std::vector<double>& m_data_I,
		const std::vector<double>& m_fitted_I)
	{
		//add fileextension
		file_name += ".chi";

		//Open file 
		std::ofstream outfile(file_name.c_str(), std::ifstream::out);

		//Check if correctly open
		if (!outfile.is_open())
		{
			std::stringstream error_stream;
			error_message(QString::fromStdString("Cannot open " + file_name));
			return;
		}

		//Set precision and notation
		outfile.precision(16);
		outfile << std::scientific;

		//Store results to file
		for (int i = 0; i<m_data_q.size(); ++i)
		{
			outfile << " " << m_data_q[i]
				<< " " << m_data_I[i]
				<< " " << m_fitted_I[i]
				<< std::endl;
		}
		//Close file
		outfile.close();
	}
	void writechisquretofile(std::string file_name, std::vector<saxs::model_sp>& result_tracer)
	{
		//add fileextension
		file_name += ".log";

		//Open file 
		std::ofstream outfile(file_name.c_str(), std::ifstream::out);

		//Check if correctly open
		if (!outfile.is_open())
		{
			std::stringstream error_stream;
			error_message(QString::fromStdString("Cannot open " + file_name));
			return;
		}

		//Set precision and notation
		outfile.precision(8);
		outfile << std::scientific;

		//Store results to file
		
		for (int j = 0; j < result_tracer[0]->m_chi_tracer.size(); j++)
		{
			for (int i = 0; i<result_tracer.size(); ++i)
			{
				outfile << "\t" << result_tracer[i]->m_chi_tracer[j];
			}
			outfile  << std::endl;
		}
		//Close file
		outfile.close();
	}



	//************************************************
	// MAIN FITTING FUNCTIONS!!!!!!!
	//************************************************
	//------------------------------------------------------------------------------
	//	Takes ftiingobject and fits model to exp by leastsquare
	//------------------------------------------------------------------------------
	void calcDebyeStackCurrentModel(fittingobject_sp& reftofittingobject)
	{
		//resortbyzcoord(reftofittingobject);
		//Now everything is ready for calculation
		//Thread container
		boost::scoped_ptr<boost::thread_group> threadGroup_sp(new boost::thread_group);

		//Mutex to protect IO operations
		boost::mutex io_mutex;

		//Initialize thread parameters
		int thread_step = static_cast<int>(reftofittingobject->m_data_q.size() / reftofittingobject->m_num_cores);
		int index_delta = reftofittingobject->m_data_q.size() - thread_step*reftofittingobject->m_num_cores;
		//Taking into account index offset due to core number
		//Putting indexboundaries into boundaries vector of size cores+1
		std::vector<int> boundaries;
		boundaries.push_back(0);
		int offset = 1;
		for (int i = 0; i < reftofittingobject->m_num_cores; i++)
		{
			if (i < index_delta) offset = 1;
			else offset = 0;
			boundaries.push_back(boundaries[i] + thread_step + offset);
		}

		//Instantiate and start threads
		for (unsigned int i = 0; i<reftofittingobject->m_num_cores; i++)
		{
			threadGroup_sp->add_thread(new boost::thread(calculate_debye,boost::ref(io_mutex), boundaries[i], boundaries[i + 1],boost::ref(reftofittingobject)));
		}

		//Wait until the end of their jobs
		if (threadGroup_sp)
			threadGroup_sp->join_all();

		//normfittingdata
		if (reftofittingobject->m_data_I.size() == 0)
		{
			double norm = *std::max_element(reftofittingobject->m_model_I.begin(), reftofittingobject->m_model_I.end());
			double norm_offset = *std::min_element(reftofittingobject->m_model_I.begin(), reftofittingobject->m_model_I.end());
			for (int i = 0; i < reftofittingobject->m_model_I.size(); i++)
			{
				reftofittingobject->m_fitted_I[i] = (reftofittingobject->m_model_I[i] - 0.99*norm_offset) / (norm - 0.99*norm_offset);
			}
			reftofittingobject->m_chi = NAN;
			reftofittingobject->m_linreg_results.first = NAN;
			reftofittingobject->m_linreg_results.second = NAN;
		}
		else
		{
			fit_modeltoexp(reftofittingobject);
		}
		//Make linear regression
	}

	//------------------------------------------------------------------------------
	//	FITTING PROCEDURE!!! Node function
	//------------------------------------------------------------------------------
	void calcDebyeFitCurrentModel(HelFitQt& helfitqt,fittingobject_sp& reftofittingobject)
	{
		int n_runs = 1;
		double starttemp = 0.7;
		double deltatemp = 0.99;
		int n_temps = 1;
		double endtemp = 0.0;
		double temp;
		double chi = reftofittingobject->m_chi;
		std::string str;

		//Now everything is ready for calculation
		//Thread container
		boost::scoped_ptr<boost::thread_group> threadGroup_sp(new boost::thread_group);

		//Mutex to protect IO operations
		boost::mutex io_mutex;

		//Initialize thread parameters
		/*int thread_step = static_cast<int>(reftofittingobject->m_model.size() / reftofittingobject->m_num_cores);
		int index_delta = reftofittingobject->m_model.size() - thread_step*reftofittingobject->m_num_cores;
		//Taking into account index offset due to core number
		//Putting indexboundaries into boundaries vector of size cores+1
		
		std::vector<int> boundaries;
		boundaries.push_back(0);
		int offset = 1;
		for (int i = 0; i < reftofittingobject->m_num_cores; i++)
		{
			if (i < index_delta) offset = 1;
			else offset = 0;
			boundaries.push_back(boundaries[i] + thread_step + offset);
		}
		*/
		calcDebyeStackCurrentModel(reftofittingobject);
		chi = reftofittingobject->m_chi;
		switch (reftofittingobject->m_weighing) {
			case (0) : str = "Chisquare (I(q)*q)\t = \t"; break;
			case (1) : str = "Chisquare (I(q)*q^2)\t = \t"; break;
			case (2) : str = "Chisquare (ln(I(q)))\t = \t"; break;
		}
		helfitqt.writetologext(str + boost::lexical_cast<std::string>(chi));

		//Preparing tracer object:
		std::vector<model_sp> result_tracer;
		for (unsigned int i = 0; i < reftofittingobject->m_num_cores; i++)
		{
			result_tracer.push_back(model_sp(new model));
			result_tracer[i]->m_chi = reftofittingobject->m_chi;
			for (int j = 0; j < reftofittingobject->m_model.size(); j++)
			{
				result_tracer[i]->m_coordinate.push_back(coordinate_sp(new coordinate(
					reftofittingobject->m_model[j]->m_x,
					reftofittingobject->m_model[j]->m_y,
					reftofittingobject->m_model[j]->m_z)));
			}
			for (int j = 0; j < reftofittingobject->m_model_I.size();j++)
			{
				result_tracer[i]->m_model_I.push_back(0.0f);
			}
		}

		for (int k = 0; k < n_temps; k++)
		{
			temp = endtemp+ (starttemp - endtemp)*std::pow(deltatemp ,k);
			//Instantiate and start threads
			for (unsigned int i = 0; i<reftofittingobject->m_num_cores; i++)
			{
				threadGroup_sp->add_thread(new boost::thread(debyefit_thread,
					boost::ref(io_mutex), boost::ref(reftofittingobject),boost::ref(result_tracer),i,n_runs, temp));
			}
			//Wait until the end of their jobs
			if (threadGroup_sp)
				threadGroup_sp->join_all();
			
			//Check which had the best result
			double chimin = reftofittingobject->m_chi;	
			int chipntr = 0;


			//Check which solution was best
			for (unsigned int i = 0; i < reftofittingobject->m_num_cores; i++)
			{
				if (result_tracer[i]->m_chi < chimin)
				{
					chimin = result_tracer[i]->m_chi;
					chipntr = i;
				}
			}

			double deltachi = double(1) - (chimin / reftofittingobject->m_chi);

			//If there is a better solution write best one into globalfittingobject
			if (reftofittingobject->m_chi != chimin)
			{
				reftofittingobject->m_chi = result_tracer[chipntr]->m_chi;
				reftofittingobject->m_model_I = result_tracer[chipntr]->m_model_I;
				for (int j = 0; j < reftofittingobject->m_model.size(); j++)
				{
					reftofittingobject->m_model[j]->m_x = result_tracer[chipntr]->m_coordinate[j]->m_x;
					reftofittingobject->m_model[j]->m_y = result_tracer[chipntr]->m_coordinate[j]->m_y;
					reftofittingobject->m_model[j]->m_z = result_tracer[chipntr]->m_coordinate[j]->m_z;
				}
			}
			//Redetermine linear regression parameters
			fit_modeltoexp(reftofittingobject);

			//Reinit result tracer
			for (int i = 0; i < reftofittingobject->m_num_cores; i++)
			{
				result_tracer[i]->m_chi = reftofittingobject->m_chi;
				result_tracer[i]->m_model_I = reftofittingobject->m_model_I;
				for (int j = 0; j < reftofittingobject->m_model.size(); j++)
				{
					result_tracer[i]->m_coordinate[j]->m_x = reftofittingobject->m_model[j]->m_x;
					result_tracer[i]->m_coordinate[j]->m_y = reftofittingobject->m_model[j]->m_y;
					result_tracer[i]->m_coordinate[j]->m_z = reftofittingobject->m_model[j]->m_z;
				}
			}

			//calcDebyeStackCurrentModel(reftofittingobject);+
			str = "Step " + boost::lexical_cast<std::string>(k + 1) + " of " + boost::lexical_cast<std::string>(n_temps) + ": ";
			str = str + "\t temp = " + boost::lexical_cast<std::string>(temp) + "\n";
			chi = reftofittingobject->m_chi;
			switch (reftofittingobject->m_weighing) {
				case (0) : str = str + "Chisquare (I(q)*q)\t = "; break;
				case (1) : str = str + "Chisquare (I(q)*q^2)\t = "; break;
				case (2) : str = str + "Chisquare (ln(I(q)))\t = "; break;
			}
			str = str + boost::lexical_cast<std::string>(chi);
			str = str + "\t rel. change = " + boost::lexical_cast<std::string>(deltachi);
			helfitqt.writetologext(str);

		}
		//Make linear regression
		//fit_modeltoexp(reftofittingobject);
	}

	//------------------------------------------------------------------------------
	//	FITTING PROCEDURE!!! Node function2
	//------------------------------------------------------------------------------
	void calcDebyeFitCurrentModel_2(HelFitQt& helfitqt, fittingobject_sp& reftofittingobject, std::string cut_filename)
	{
		//TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!
		int n_runs = 4;
		double starttemp = 0.7;
		double deltatemp = 0.99;
		double endtemp = 0.0;
		//TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!

		double chi = reftofittingobject->m_chi;
		//Helper string
		std::string str;

		//Now everything is ready for calculation
		//Thread container
		boost::scoped_ptr<boost::thread_group> threadGroup_sp(new boost::thread_group);

		//Mutex to protect IO operations
		boost::mutex io_mutex;

		//Initialize thread parameters
		/*int thread_step = static_cast<int>(reftofittingobject->m_model.size() / reftofittingobject->m_num_cores);
		int index_delta = reftofittingobject->m_model.size() - thread_step*reftofittingobject->m_num_cores;
		//Taking into account index offset due to core number
		//Putting indexboundaries into boundaries vector of size cores+1

		std::vector<int> boundaries;
		boundaries.push_back(0);
		int offset = 1;
		for (int i = 0; i < reftofittingobject->m_num_cores; i++)
		{
		if (i < index_delta) offset = 1;
		else offset = 0;
		boundaries.push_back(boundaries[i] + thread_step + offset);
		}
		*/
		calcDebyeStackCurrentModel(reftofittingobject);
		chi = reftofittingobject->m_chi;
		switch (reftofittingobject->m_weighing) {
		case (0) : str = "Chisquare (I(q)*q)\t = \t"; break;
		case (1) : str = "Chisquare (I(q)*q^2)\t = \t"; break;
		case (2) : str = "Chisquare (ln(I(q)))\t = \t"; break;
		}
		helfitqt.writetologext("Initial " +str + boost::lexical_cast<std::string>(chi));

		//Preparing tracer object:
		std::vector<model_sp> result_tracer;
		for (unsigned int i = 0; i < reftofittingobject->m_num_cores; i++)
		{
			result_tracer.push_back(model_sp(new model));
			result_tracer[i]->m_chi = reftofittingobject->m_chi;
			result_tracer[i]->m_chi_tracer.push_back(reftofittingobject->m_chi);
			for (int j = 0; j < reftofittingobject->m_model.size(); j++)
			{
				result_tracer[i]->m_coordinate.push_back(coordinate_sp(new coordinate(
					reftofittingobject->m_model[j]->m_x,
					reftofittingobject->m_model[j]->m_y,
					reftofittingobject->m_model[j]->m_z)));
			}
			for (int j = 0; j < reftofittingobject->m_model_I.size(); j++)
			{
				result_tracer[i]->m_model_I.push_back(0.0f);
				result_tracer[i]->m_fitted_I.push_back(0.0f);
			}
		}

		//Instantiate and start threads
		for (unsigned int i = 0; i<reftofittingobject->m_num_cores; i++)
		{
			threadGroup_sp->add_thread(new boost::thread(debyefit_thread_2,
				boost::ref(helfitqt), boost::ref(reftofittingobject), boost::ref(result_tracer), i, n_runs, starttemp, endtemp,deltatemp));
		}
		//Wait until the end of their jobs
		if (threadGroup_sp)
			threadGroup_sp->join_all();

		//Save results to files
		std::string tempstring;
		for (unsigned int i = 0; i < reftofittingobject->m_num_cores; i++)
		{
			tempstring = cut_filename + "_" + std::to_string(i);
			writemodeltopdb(tempstring, result_tracer[i]->m_coordinate);
			writedatatochi(tempstring, reftofittingobject->m_data_q, reftofittingobject->m_data_I, result_tracer[i]->m_fitted_I);
		}

		//Check which had the best result
		double chimin = reftofittingobject->m_chi;
		int chipntr = 0;
		
		//Check which solution was best
		for (unsigned int i = 0; i < reftofittingobject->m_num_cores; i++)
		{
			if (result_tracer[i]->m_chi < chimin)
			{
				chimin = result_tracer[i]->m_chi;
				chipntr = i;
			}
		}

		//If there is a better solution write best one into globalfittingobject
		if (reftofittingobject->m_chi != chimin)
		{
			reftofittingobject->m_chi = result_tracer[chipntr]->m_chi;
			reftofittingobject->m_model_I = result_tracer[chipntr]->m_model_I;
			for (int j = 0; j < reftofittingobject->m_model.size(); j++)
			{
				reftofittingobject->m_model[j]->m_x = result_tracer[chipntr]->m_coordinate[j]->m_x;
				reftofittingobject->m_model[j]->m_y = result_tracer[chipntr]->m_coordinate[j]->m_y;
				reftofittingobject->m_model[j]->m_z = result_tracer[chipntr]->m_coordinate[j]->m_z;
			}
		}

		//Output to extfield
		std::stringstream ss;
		ss << "The best result has core  " << chipntr << " with chi=";
		ss << std::scientific;
		ss << chimin;
		str = ss.str();
		helfitqt.writetologext(str);

		//Make linear regression
		//fit_modeltoexp(reftofittingobject);
	}

	//------------------------------------------------------------------------------
	//	Generate random cylindrical coordinates into coordinate structure
	//------------------------------------------------------------------------------
	void gen_cyl_randomseed(std::vector<coordinate_sp>& cylindrical_coord, double height, double diameter, const int nr_atoms)
	{
		//pseudo coords: m_x = r, m_y = phi, m_z = z

		typedef boost::mt19937 RNGType;
		RNGType rng(std::time(0));
		boost::uniform_real<> zero_to_one(0, 1);
		boost::variate_generator<RNGType, boost::uniform_real<>> dice(rng, zero_to_one);

		for (int i = 0; i < nr_atoms; i++)
			cylindrical_coord.push_back(coordinate_sp(new coordinate(
				double(std::pow(dice(),0.5) * diameter/2),
				double(dice() * 2 * 3.14159265359),
				double(dice() * height)
				)));
	}

}//End of SAXS namespace




/*
//------------------------------------------------------------------------------
//	Import coordinates from file
//------------------------------------------------------------------------------
void import_coordinates(const std::string& file_name, std::vector<coordinate_sp>& coordinates)
{
	//Create a logger source for import funtion
	//boost::log::sources::logger lg;

	//BOOST_LOG(lg) << "Importing coordinates from: " << file_name;

	//Open file 
	std::ifstream infile(file_name.c_str(), std::ifstream::in);

	//Check if correctly open
	if (!infile.is_open())
	{
		std::stringstream error_stream;
		error_stream << "Cannot open \'" << file_name << "\'" << std::endl;
		throw std::runtime_error(error_stream.str());
	}

	//Temporary variable used during import
	double x = 0;
	double y = 0;
	double z = 0;

	//Temporary line used during import
	std::string line;

	//Read and parse file line by line
	int counter = 0;
	while (getline(infile, line))
	{
		counter += 1;

		std::istringstream iss(line);
		iss >> x >> y >> z;

		//Check for errors
		if (iss.fail())
		{
			std::stringstream error_stream;
			error_stream << "Error reading file " << file_name
				<< " at line " << counter << std::endl;
			throw std::runtime_error(error_stream.str());
		}

		//Store coordinates
		coordinates.push_back( coordinate_sp( new coordinate(x, y, z) ) );
	}

	//Close file
	infile.close();
}

//------------------------------------------------------------------------------
//	Store results to file
//------------------------------------------------------------------------------
void store_results(const std::string& file_name, std::vector<scatteringdata_sp>& results)
{
	//Create a logger source for import funtion
	//boost::log::sources::logger lg;

	//BOOST_LOG(lg) << "Storing results to: " << file_name;

	//Open file 
	std::ofstream outfile(file_name.c_str(), std::ifstream::out);

	//Check if correctly open
	if (!outfile.is_open())
	{
		std::stringstream error_stream;
		error_stream << "Cannot open \'" << file_name << "\'" << std::endl;
		throw std::runtime_error(error_stream.str());
	}

	//Set precision and notation
	outfile.precision(6);
	outfile << std::scientific;

	//Store results to file
	std::vector<scatteringdata_sp>::iterator it;
	for (it = results.begin(); it != results.end(); ++it)
	{
		outfile << (*it)->m_q << "\t" << (*it)->m_sq << std::endl;
	}

	//Close file
	outfile.close();
}

}	//End of namespace
*/
#endif /*!DATA_UTILS_H*/