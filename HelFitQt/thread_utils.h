#ifndef THREAD_UTILS_H
#define THREAD_UTILS_H

#include "data_def.h"
#include "data_utils.h"
#include <QtCore>
#include <QApplication>

#include <boost\shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/thread.hpp>
#include <boost/math/special_functions/sinc.hpp>

// Used in randomization
#include <iostream>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>

#include <vector>
#include <random>
#include <stdint.h>

namespace saxs
{
	//helper function
	uint32_t randseed = 7;  // 100% random seed value

	double xor128(void) {
		randseed ^= randseed << 13;
		randseed ^= randseed >> 17;
		randseed ^= randseed << 5;
		return double(randseed) / 2 / double(2147483647);
	}

	//------------------------------------------------------------------------------
	//	Calculate debye function
	//------------------------------------------------------------------------------
	void calculate_debye(boost::mutex& io_mutex, int start_index, int end_index,fittingobject_sp& fittingobject)
	{
		//Thread local temporary vector
		std::vector<double> local_results(fittingobject->m_data_q.size(), 0.0f);

		//Keep size of lookup table
		std::size_t sinc_lookup_size = HelFitQt::sinc_lookup.size();

		//Temporary variables used in Debye function
		double r = 0;
		double qr = 0;

		//loca coordinate vector
		std::vector<coordinate_sp> coordinates = fittingobject->m_model;

		//Compute Debye function [Max Code]
		//Loop around every atom in BuildingBlock
		for (std::size_t i = 0; i < (coordinates.size()); ++i)
		{
			//BuildingBlockCodnition
			if (i < (coordinates.size() - 1))
			{
				//loop over j points
				for (std::size_t j = i + 1; j < (coordinates.size()); ++j)
				{
					r = std::pow(
						std::pow(coordinates[i]->m_x - coordinates[j]->m_x, 2) +
						std::pow(coordinates[i]->m_y - coordinates[j]->m_y, 2) +
						std::pow(coordinates[i]->m_z - coordinates[j]->m_z, 2),
						0.5);//loop over q-points

					for (int k = start_index; k < end_index; k++)
					{
						qr = fittingobject->m_data_q[k] * r * 1000 + 0.5;
						//Lookup Condition
						if (qr < sinc_lookup_size)
						{
							local_results[k] += HelFitQt::sinc_lookup[(int)qr] * 2.0*double(fittingobject->m_num_stacks);
						}
						else
						{
							local_results[k] += boost::math::sinc_pi(double(qr) / 1000)*2.0*double(fittingobject->m_num_stacks);
						}
					}//k index loop
				}//j index loop
			}//if Building block Condition

			 //Loop over stacks
			for (int l = 1; l < fittingobject->m_num_stacks; ++l)
			{
				//loop over j - radial pairs
				for (std::size_t j = 0; j < (coordinates.size()); ++j)
				{
					r = std::pow(
						std::pow(coordinates[i]->m_x - coordinates[j]->m_x, 2) +
						std::pow(coordinates[i]->m_y - coordinates[j]->m_y, 2) +
						std::pow(coordinates[i]->m_z - (coordinates[j]->m_z + double(l*fittingobject->m_stack_spacing)), 2),
						0.5);
					//loop over q-points
					for (int k = start_index; k < end_index; k++)
					{
						qr = fittingobject->m_data_q[k] * r * 1000 + 0.5;
						//Lookup Condition
						if (qr < sinc_lookup_size)
						{
							local_results[k] += 2 * HelFitQt::sinc_lookup[(int)qr] * double(fittingobject->m_num_stacks - l);
						}
						else
						{
							local_results[k] += 2 * boost::math::sinc_pi(double(qr) / 1000)* double(fittingobject->m_num_stacks - l);
						}
					}// k index
				}//j index loop
			}//l index loop
			 //Notify user to current status


		}// i index loop

		 //Lock mutex
		boost::mutex::scoped_lock lock(io_mutex);

		//Copy data from local temporary vector to results
		for (unsigned int k = start_index; k < end_index; k++)
		{
			fittingobject->m_model_I[k] += local_results[k] + double(fittingobject->m_num_stacks*fittingobject->m_model.size());
		}

		//Unlock mutex
		lock.unlock();
	}

	//*************************************************************************************
	// Debye fitting functions
	//*************************************************************************************
	//------------------------------------------------------------------------------
	//	linear fit of two vectors without public vars
	//------------------------------------------------------------------------------
	double chisquare(std::vector<double>& exp, std::vector<double>& model, std::vector<double>& x, int& weighing)
	{
		//Chisquare between data1(exp) and data2(model)
		double chi = 0;
		double weight = 0;
		for (int i = 0; i < exp.size(); i++)
		{
			chi += std::pow((exp[i] - model[i])/0.01, 2)/(exp[i]);
			weight += 1;
		}
		return chi / weight;
	}

	std::pair<double, double> lindatafit(std::vector<double>& data, std::vector<double>& fit, std::vector<double>& weight)
	{
		double f2sum = 0;
		double fsum = 0;
		double dfsum = 0;
		double dsum = 0;
		double weightsum = 0;

		for (int i = 0; i < data.size(); i++)
		{
			f2sum = f2sum + fit[i] * fit[i]* weight[i];
			fsum = fsum + fit[i] * weight[i];
			dfsum = dfsum + data[i] * fit[i] * weight[i];
			dsum = dsum + data[i] * weight[i];
			weightsum = weightsum + weight[i];
		}
		double k = (weightsum*dfsum - dsum*fsum) / (weightsum*f2sum - fsum*fsum);
		double d = (dfsum - f2sum*k) / fsum;

		std::pair<double, double> fitresults(d, k);
		return fitresults;
	}

	void fit_modeltoexp(fittingobject_sp& reftofittingobject)
	{
		//Curve weighting: 0...none,1...q^2, 2...ln
		std::vector<double> exp;
		std::vector<double> model;
		std::vector<double> weight;
		std::pair<double, double> fitparams;

		reftofittingobject->m_norm = *std::max_element(reftofittingobject->m_model_I.begin(), reftofittingobject->m_model_I.end());
		reftofittingobject->m_norm_offset = *std::min_element(reftofittingobject->m_model_I.begin(), reftofittingobject->m_model_I.end());
		for (int i = 0; i < reftofittingobject->m_model_I.size(); i++)
		{
			reftofittingobject->m_fitted_I[i] = (reftofittingobject->m_model_I[i] - 
				0.99*reftofittingobject->m_norm_offset) / (reftofittingobject->m_norm - 0.99*reftofittingobject->m_norm_offset);
		}

		//Loading data in standard vectors and finding minimum value
		for (int i = 0; i < reftofittingobject->m_data_q.size(); i++)
		{
			exp.push_back(reftofittingobject->m_data_I[i]);
			model.push_back(reftofittingobject->m_fitted_I[i]);
		}
		weight.resize(model.size());
		//Apply weighting conditions
		for (int i = 0; i < reftofittingobject->m_data_q.size(); i++)
		{
			if (reftofittingobject->m_weighing ==0 )
			{
				weight[i] = reftofittingobject->m_data_q[i]/std::pow(exp[i],0.5);
			}
			if (reftofittingobject->m_weighing == 1)
			{
				weight[i] = std::pow(reftofittingobject->m_data_q[i], 2)/std::pow(exp[i], 0.5);
			}
			if (reftofittingobject->m_weighing == 2)
			{
				exp[i] = std::log(exp[i]);
				model[i] = std::log(model[i]);
				weight[i] =1;
			}
		}

		//GetFitparams
		fitparams = lindatafit(exp, model,weight);
		reftofittingobject->m_linreg_results = fitparams;
		if (reftofittingobject->m_weighing < 2)
		{
			for (int i = 0; i < reftofittingobject->m_data_q.size(); i++)
			{
				reftofittingobject->m_fitted_I[i] = fitparams.first + fitparams.second * reftofittingobject->m_fitted_I[i];
			}
		}
		else
		{
			for (int i = 0; i < reftofittingobject->m_data_q.size(); i++)
			{
				reftofittingobject->m_fitted_I[i] = std::exp(fitparams.second*std::log(reftofittingobject->m_fitted_I[i]) + fitparams.first);
			}
		}

		reftofittingobject->m_chi = chisquare(reftofittingobject->m_data_I, reftofittingobject->m_fitted_I,
			reftofittingobject->m_data_q, reftofittingobject->m_weighing);
	}

	void recalc_fitted_i(fittingobject_sp& localfittingobject_sp)
	{
		std::pair<double, double> fitparams = localfittingobject_sp->m_linreg_results;
		if (localfittingobject_sp->m_weighing < 2)
		{
			for (int i = 0; i < localfittingobject_sp->m_data_q.size(); i++)
			{
				localfittingobject_sp->m_fitted_I[i] = (localfittingobject_sp->m_model_I[i] -
					0.99*localfittingobject_sp->m_norm_offset) / (localfittingobject_sp->m_norm - 0.99*localfittingobject_sp->m_norm_offset);
				localfittingobject_sp->m_fitted_I[i] = localfittingobject_sp->m_linreg_results.first + localfittingobject_sp->m_linreg_results.second * localfittingobject_sp->m_fitted_I[i];
			}
		}
		else
		{
			for (int i = 0; i < localfittingobject_sp->m_data_q.size(); i++)
			{
				localfittingobject_sp->m_fitted_I[i] = (localfittingobject_sp->m_model_I[i] -
					0.99*localfittingobject_sp->m_norm_offset) / (localfittingobject_sp->m_norm - 0.99*localfittingobject_sp->m_norm_offset);
				localfittingobject_sp->m_fitted_I[i] = std::exp(localfittingobject_sp->m_linreg_results.second*std::log(localfittingobject_sp->m_fitted_I[i]) + localfittingobject_sp->m_linreg_results.first);
			}
		}
		localfittingobject_sp->m_chi = chisquare(localfittingobject_sp->m_data_I, localfittingobject_sp->m_fitted_I,
			localfittingobject_sp->m_data_q, localfittingobject_sp->m_weighing);
	}

	//Recalculates the scattering pattern if changecood is moved by movement
	//Result will be put into references new_model_i
	void thread_recalc_change(fittingobject_sp& reftofittingobject, int& changedcoord, coordinate_sp& movement)
	{
		//tempvars
		double r_old = 0;
		double qr_old = 0;
		double sinc_qr_old = 0;
		double r_new = 0;
		double qr_new = 0;
		double sinc_qr_new = 0;
		//tempvars for stack
		double r_old_st = 0;
		double qr_old_st = 0;
		double sinc_qr_old_st = 0;
		double r_new_st = 0;
		double qr_new_st = 0;
		double sinc_qr_new_st = 0;
		int sinc_lookup_size = HelFitQt::sinc_lookup.size();

		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
		{
			if (i != changedcoord)
			{
				//Building block contribution
				r_old = std::pow(
					std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
					std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
					std::pow(reftofittingobject->m_model[i]->m_z - reftofittingobject->m_model[changedcoord]->m_z, 2),
					0.5);
				r_new = std::pow(
					std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
					std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
					std::pow(reftofittingobject->m_model[i]->m_z - (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
					0.5);
				//loop over q-points
				for (int q = 0; q < reftofittingobject->m_data_q.size(); q++)
				{
					//Lookup Condition
					qr_old = reftofittingobject->m_data_q[q] * r_old * 1000 + 0.5;
					if (qr_old < sinc_lookup_size) sinc_qr_old = HelFitQt::sinc_lookup[(int)qr_old];
					else sinc_qr_old = boost::math::sinc_pi(double(qr_old) / 1000);

					qr_new = reftofittingobject->m_data_q[q] * r_new * 1000 + 0.5;
					if (qr_new < sinc_lookup_size) sinc_qr_new = HelFitQt::sinc_lookup[(int)qr_new];
					else sinc_qr_new = boost::math::sinc_pi(double(qr_new) / 1000);
					//Lookup Condition
					reftofittingobject->m_model_I[q] += 2 * (sinc_qr_new - sinc_qr_old)*double(reftofittingobject->m_num_stacks);
				}//enfor q
				 //Stack contribution
				if (reftofittingobject->m_num_stacks>1)
				{
					for (int j = 1; j < reftofittingobject->m_num_stacks; j++)
					{
						r_old = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
							std::pow(reftofittingobject->m_model[i]->m_z + double(j) * reftofittingobject->m_stack_spacing 
																				- reftofittingobject->m_model[changedcoord]->m_z, 2),
							0.5);
						r_new = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
							std::pow(reftofittingobject->m_model[i]->m_z + double(j) * reftofittingobject->m_stack_spacing 
																			- (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
							0.5);
						r_old_st = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
							std::pow(reftofittingobject->m_model[i]->m_z - double(j) * reftofittingobject->m_stack_spacing 
																				- reftofittingobject->m_model[changedcoord]->m_z, 2),
							0.5);
						r_new_st = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
							std::pow(reftofittingobject->m_model[i]->m_z - double(j) * reftofittingobject->m_stack_spacing 
																				- (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
							0.5);
						//loop over q-points
						for (int q = 0; q < reftofittingobject->m_data_q.size(); q++)
						{
							//Lookup Condition
							qr_old = reftofittingobject->m_data_q[q] * r_old * 1000 + 0.5;
							if (qr_old < sinc_lookup_size) sinc_qr_old = HelFitQt::sinc_lookup[(int)qr_old];
							else sinc_qr_old = boost::math::sinc_pi(double(qr_old) / 1000);

							qr_new = reftofittingobject->m_data_q[q] * r_new * 1000 + 0.5;
							if (qr_new < sinc_lookup_size) sinc_qr_new = HelFitQt::sinc_lookup[(int)qr_new];
							else sinc_qr_new = boost::math::sinc_pi(double(qr_new) / 1000);

							qr_old_st = reftofittingobject->m_data_q[q] * r_old_st * 1000 + 0.5;
							if (qr_old_st < sinc_lookup_size) sinc_qr_old_st = HelFitQt::sinc_lookup[(int)qr_old_st];
							else sinc_qr_old_st = boost::math::sinc_pi(double(qr_old_st) / 1000);

							qr_new_st = reftofittingobject->m_data_q[q] * r_new_st * 1000 + 0.5;
							if (qr_new_st < sinc_lookup_size) sinc_qr_new_st = HelFitQt::sinc_lookup[(int)qr_new_st];
							else sinc_qr_new_st = boost::math::sinc_pi(double(qr_new_st) / 1000);
							//Lookup Condition

							reftofittingobject->m_model_I[q] += 2 * (sinc_qr_new - sinc_qr_old)*double(reftofittingobject->m_num_stacks - j);

							reftofittingobject->m_model_I[q] += 2 * (sinc_qr_new_st - sinc_qr_old_st)*double(reftofittingobject->m_num_stacks - j);
						}//enfor q
					}//endfor stacks
				}//endif stack condition
			}//endif i!=changedcoord
		}//endfor i

	}

	//************************************************
	// MAIN FITTING THREAD!!!!!!!
	//************************************************
	void debyefit_thread(boost::mutex& io_mutex, fittingobject_sp& reftofittingobject, std::vector<model_sp>& result_tracer, int coreid, int n_runs, double temp)
	{
		//Make local copy of fitting object so changes don't affect global object
		fittingobject localfittingobject = *reftofittingobject;
		boost::shared_ptr<fittingobject> localfittingobject_sp = boost::make_shared<fittingobject>(localfittingobject);
		localfittingobject_sp->m_model.clear();
		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
			localfittingobject_sp->m_model.push_back(coordinate_sp(new coordinate(reftofittingobject->m_model[i]->m_x,
																			reftofittingobject->m_model[i]->m_y,
																			reftofittingobject->m_model[i]->m_z)));

		std::vector<double> old_model_i;
		old_model_i.resize(localfittingobject_sp->m_model_I.size());

		coordinate move = (0, 0, 0);
		boost::shared_ptr<coordinate> movement = boost::make_shared<coordinate>(move);
		double old_chi;

		//heating up randomization
		randseed = (unsigned)time(NULL);
		double rand_scalar_r = temp * localfittingobject_sp->m_diameter / 2;
		double rand_scalar_h = temp * localfittingobject_sp->m_stack_spacing / 2;



		for (int i = 0; i < n_runs; i++)
		{
			for (int j = 0; j < reftofittingobject->m_model.size(); j++)
			{
				old_model_i = localfittingobject_sp->m_model_I;
				old_chi = localfittingobject_sp->m_chi;
				movement->m_x = rand_scalar_r*(0.5 - xor128());
				movement->m_y = rand_scalar_r*(0.5 - xor128());
				movement->m_z = rand_scalar_h*(0.5 - xor128());
				thread_recalc_change(localfittingobject_sp, j, movement);
				recalc_fitted_i(localfittingobject_sp);
				//fit_modeltoexp(localfittingobject_sp);
				if (localfittingobject_sp->m_chi < old_chi*(double(1)))
				{
					localfittingobject_sp->m_model[j]->m_x += movement->m_x;
					localfittingobject_sp->m_model[j]->m_y += movement->m_y;
					localfittingobject_sp->m_model[j]->m_z += movement->m_z;
				}
				else
				{
					localfittingobject_sp->m_model_I = old_model_i;
					localfittingobject_sp->m_chi = old_chi;
				}
			}//endfor j coords

		}//endfor i runs

		 //Lock mutex
		//boost::mutex::scoped_lock lock(io_mutex);

		//Copy data from local temporary vector to results
		for (unsigned int k = 0; k < reftofittingobject->m_model.size(); k++)
		{
			result_tracer[coreid]->m_coordinate[k]->m_x = localfittingobject_sp->m_model[k]->m_x;
			result_tracer[coreid]->m_coordinate[k]->m_y = localfittingobject_sp->m_model[k]->m_y;
			result_tracer[coreid]->m_coordinate[k]->m_z = localfittingobject_sp->m_model[k]->m_z;
		}
		result_tracer[coreid]->m_chi = localfittingobject_sp->m_chi;
		result_tracer[coreid]->m_model_I = localfittingobject_sp->m_model_I;
		//Unlock mutex
		//lock.unlock();
	}

	void debyefit_thread_2(HelFitQt& helfitqt, fittingobject_sp& reftofittingobject, std::vector<model_sp>& result_tracer,
																		int coreid, int n_runs, double starttemp, double endtemp, double deltatemp)
	{
		//Make local copy of fitting object so changes don't affect global object
		fittingobject localfittingobject = *reftofittingobject;
		boost::shared_ptr<fittingobject> localfittingobject_sp = boost::make_shared<fittingobject>(localfittingobject);

		//Deepcopy model
		localfittingobject_sp->m_model.clear();
		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
			localfittingobject_sp->m_model.push_back(coordinate_sp(new coordinate(reftofittingobject->m_model[i]->m_x,
				reftofittingobject->m_model[i]->m_y,
				reftofittingobject->m_model[i]->m_z)));
		
		//Deppcopy data
		int datasize = localfittingobject_sp->m_model_I.size();
		localfittingobject_sp->m_model_I.clear();
		localfittingobject_sp->m_fitted_I.clear();
		localfittingobject_sp->m_data_I.clear();
		for (int i = 0; i < datasize; i++)
		{
			localfittingobject_sp->m_model_I.push_back(reftofittingobject->m_model_I[i]);
			localfittingobject_sp->m_fitted_I.push_back(reftofittingobject->m_fitted_I[i]);
			localfittingobject_sp->m_data_I.push_back(reftofittingobject->m_data_I[i]);
		}

		//Define helper vector for temporary storage of results
		std::vector<double> old_model_i;
		old_model_i.resize(localfittingobject_sp->m_model_I.size());

		//Initialize futher variables
		coordinate move = (0, 0, 0);
		boost::shared_ptr<coordinate> movement = boost::make_shared<coordinate>(move);
		double old_chi;
		double temp;
		double temp_z_coord;

		for (int i = 0; i < n_runs; i++)
		{
		temp = endtemp + (starttemp - endtemp)*std::pow(deltatemp, i);
		//heating up randomization
		randseed = (unsigned)time(NULL);
		double rand_scalar_r = temp * localfittingobject_sp->m_diameter / 2;
		double rand_scalar_h = temp * localfittingobject_sp->m_stack_spacing / 2;

			for (int j = 0; j < reftofittingobject->m_model.size(); j++)
			{
				old_model_i = localfittingobject_sp->m_model_I;
				old_chi = localfittingobject_sp->m_chi;
				movement->m_x = rand_scalar_r*(0.5 - xor128());
				movement->m_y = rand_scalar_r*(0.5 - xor128());
				movement->m_z = rand_scalar_h*(0.5 - xor128());
				thread_recalc_change(localfittingobject_sp, j, movement);
				recalc_fitted_i(localfittingobject_sp);
				//fit_modeltoexp(localfittingobject_sp);
				if (localfittingobject_sp->m_chi < old_chi*(double(1)))
				{
					localfittingobject_sp->m_model[j]->m_x += movement->m_x;
					localfittingobject_sp->m_model[j]->m_y += movement->m_y;
					temp_z_coord = localfittingobject_sp->m_model[j]->m_z + movement->m_z;
					if (temp_z_coord < 0)  movement->m_z += localfittingobject_sp->m_stack_spacing;
					if (temp_z_coord > localfittingobject_sp->m_stack_spacing)  movement->m_z -= localfittingobject_sp->m_stack_spacing;
					localfittingobject_sp->m_model[j]->m_z += movement->m_z;
				}
				else
				{
					localfittingobject_sp->m_model_I = old_model_i;
					localfittingobject_sp->m_chi = old_chi;
				}
			}//endfor j coords
			fit_modeltoexp(localfittingobject_sp);

			result_tracer[coreid]->m_chi_tracer.push_back(localfittingobject_sp->m_chi);
			std::string str;
			std::stringstream ss;
			ss << "Core " << coreid << " - step " << i << ": ";
			ss << std::scientific;
			ss << localfittingobject_sp->m_chi;
			str = ss.str();
			//helfitqt.writetologext(str);
		}//endfor i runs

		 //Lock mutex
		 //boost::mutex::scoped_lock lock(io_mutex);

		 //Copy data from local temporary vector to results
		for (unsigned int k = 0; k < reftofittingobject->m_model.size(); k++)
		{
			result_tracer[coreid]->m_coordinate[k]->m_x = localfittingobject_sp->m_model[k]->m_x;
			result_tracer[coreid]->m_coordinate[k]->m_y = localfittingobject_sp->m_model[k]->m_y;
			result_tracer[coreid]->m_coordinate[k]->m_z = localfittingobject_sp->m_model[k]->m_z;
		}
		result_tracer[coreid]->m_chi = localfittingobject_sp->m_chi;
		result_tracer[coreid]->m_model_I = localfittingobject_sp->m_model_I;
		result_tracer[coreid]->m_fitted_I = localfittingobject_sp->m_fitted_I;
		//Unlock mutex
		//lock.unlock();
	}

}	//End of namespace

	//Debyefit as part of DebyeFitThread-class
	void DebyeFitThread::run()
	{
		//Make local copy of fitting object so changes don't affect global object
		saxs::fittingobject localfittingobject = *reftofittingobject;
		boost::shared_ptr<saxs::fittingobject> localfittingobject_sp = boost::make_shared<saxs::fittingobject>(localfittingobject);

		//Deepcopy model
		localfittingobject_sp->m_model.clear();
		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
			localfittingobject_sp->m_model.push_back(saxs::coordinate_sp(new saxs::coordinate(reftofittingobject->m_model[i]->m_x,
				reftofittingobject->m_model[i]->m_y,
				reftofittingobject->m_model[i]->m_z)));

		//Deppcopy data
		int datasize = localfittingobject_sp->m_model_I.size();
		localfittingobject_sp->m_model_I.clear();
		localfittingobject_sp->m_fitted_I.clear();
		localfittingobject_sp->m_data_I.clear();
		for (int i = 0; i < datasize; i++)
		{
			localfittingobject_sp->m_model_I.push_back(reftofittingobject->m_model_I[i]);
			localfittingobject_sp->m_fitted_I.push_back(reftofittingobject->m_fitted_I[i]);
			localfittingobject_sp->m_data_I.push_back(reftofittingobject->m_data_I[i]);
		}

		//Define helper vector for temporary storage of results
		std::vector<double> old_model_i;
		old_model_i.resize(localfittingobject_sp->m_model_I.size());

		//Initialize futher variables
		saxs::coordinate move = (0, 0, 0);
		boost::shared_ptr<saxs::coordinate> movement = boost::make_shared<saxs::coordinate>(move);
		double old_chi;
		double temp;
		double temp_z_coord;
		double temp_r_coord;

		for (int i = 0; i < n_runs; i++)
		{
			temp = endtemp + (starttemp - endtemp)*std::pow(deltatemp, i);
			//heating up randomization
			saxs::randseed = (unsigned)time(NULL);
			double rand_scalar_r = temp * localfittingobject_sp->m_diameter / 2;
			double rand_scalar_h = temp * localfittingobject_sp->m_stack_spacing / 2;

			for (int j = 0; j < reftofittingobject->m_model.size(); j++)
			{
				old_model_i = localfittingobject_sp->m_model_I;
				old_chi = localfittingobject_sp->m_chi;
				movement->m_x = rand_scalar_r*(0.5 - saxs::xor128());
				movement->m_y = rand_scalar_r*(0.5 - saxs::xor128());
				movement->m_z = rand_scalar_h*(0.5 - saxs::xor128());
				/*temp_r_coord =std::pow((movement->m_x + localfittingobject_sp->m_model[j]->m_x), 2) +
									std::pow((movement->m_y + localfittingobject_sp->m_model[j]->m_y), 2);
				if (temp_r_coord > std::pow(localfittingobject_sp->m_diameter, 2))
				{
					movement->m_x = (0.5 - saxs::xor128())*2*localfittingobject_sp->m_diameter - localfittingobject_sp->m_model[j]->m_x;
					movement->m_y = (0.5 - saxs::xor128())*2*localfittingobject_sp->m_diameter - localfittingobject_sp->m_model[j]->m_y;
				}*/
				saxs::thread_recalc_change(localfittingobject_sp, j, movement);
				saxs::recalc_fitted_i(localfittingobject_sp);
				//fit_modeltoexp(localfittingobject_sp);
				if (localfittingobject_sp->m_chi < old_chi*(double(1)))
				{
					localfittingobject_sp->m_model[j]->m_x += movement->m_x;
					localfittingobject_sp->m_model[j]->m_y += movement->m_y;
					temp_z_coord = localfittingobject_sp->m_model[j]->m_z + movement->m_z;
					if (temp_z_coord < 0)  movement->m_z += localfittingobject_sp->m_stack_spacing;
					if (temp_z_coord > localfittingobject_sp->m_stack_spacing)  movement->m_z -= localfittingobject_sp->m_stack_spacing;
					localfittingobject_sp->m_model[j]->m_z += movement->m_z;
				}
				else
				{
					localfittingobject_sp->m_model_I = old_model_i;
					localfittingobject_sp->m_chi = old_chi;
				}
			}//endfor j coords
			fit_modeltoexp(localfittingobject_sp);

			result_tracer[coreid]->m_chi_tracer.push_back(localfittingobject_sp->m_chi);
			if (coreid==0)
			{
				std::string str;
				std::stringstream ss;
				ss << "Core " << coreid << " - step " << i << ": ";
				ss << std::scientific;
				ss << localfittingobject_sp->m_chi;
				str = ss.str();
				emit sig_writelogext(str);
			}
		}//endfor i runs

		 //Copy data from local temporary vector to results
		for (unsigned int k = 0; k < reftofittingobject->m_model.size(); k++)
		{
			result_tracer[coreid]->m_coordinate[k]->m_x = localfittingobject_sp->m_model[k]->m_x;
			result_tracer[coreid]->m_coordinate[k]->m_y = localfittingobject_sp->m_model[k]->m_y;
			result_tracer[coreid]->m_coordinate[k]->m_z = localfittingobject_sp->m_model[k]->m_z;
		}
		result_tracer[coreid]->m_chi = localfittingobject_sp->m_chi;
		result_tracer[coreid]->m_model_I = localfittingobject_sp->m_model_I;
		result_tracer[coreid]->m_fitted_I = localfittingobject_sp->m_fitted_I;

		//Send finish signal
		emit sig_fittingthread_done(coreid);
	}

#endif /*!THREAD_UTILS_H*/

/*Old Procedures
void calculate_debye_old(boost::mutex& io_mutex, int start_index, int end_index, int n_stacks, double stack_spacing,
const std::vector<double>& sinc_lookup, const std::vector<coordinate_sp>& coordinates,
std::vector<scatteringdata_sp>& results)
{
//Thread local temporary vector
std::vector<double> local_results(results.size(), 0.0f);

//Keep size of lookup table
std::size_t sinc_lookup_size = sinc_lookup.size();

//Temporary variables used in Debye function
double r = 0;
double qr = 0;

//Compute Debye function [Max Code]
//Loop around every atom in BuildingBlock
for (std::size_t i = 0; i < (coordinates.size()); ++i)
{
//BuildingBlockCodnition
if (i < (coordinates.size() - 1))
{
//loop over j points
for (std::size_t j = i + 1; j < (coordinates.size()); ++j)
{
r = std::pow(
std::pow(coordinates[i]->m_x - coordinates[j]->m_x, 2) +
std::pow(coordinates[i]->m_y - coordinates[j]->m_y, 2) +
std::pow(coordinates[i]->m_z - coordinates[j]->m_z, 2),
0.5);//loop over q-points

for (int k = start_index; k < end_index; k++)
{
qr = results[k]->m_q * r * 1000 + 0.5;
//Lookup Condition
if (qr < sinc_lookup_size)
{
local_results[k] += sinc_lookup[(int)qr]*2.0*double(n_stacks);
}
else
{
local_results[k] += boost::math::sinc_pi(double(qr) / 1000)*2.0*double(n_stacks);
}
}//k index loop
}//j index loop
}//if Building block Condition

//Loop over stacks
for (int l = 1; l < n_stacks; ++l)
{
//loop over j - radial pairs
for (std::size_t j = 0; j < (coordinates.size()); ++j)
{
r = std::pow(
std::pow(coordinates[i]->m_x - coordinates[j]->m_x, 2) +
std::pow(coordinates[i]->m_y - coordinates[j]->m_y, 2) +
std::pow(coordinates[i]->m_z - (coordinates[j]->m_z + double(l*stack_spacing)), 2),
0.5);
//loop over q-points
for (int k = start_index; k < end_index; k++)
{
qr = results[k]->m_q * r * 1000 + 0.5;
//Lookup Condition
if (qr < sinc_lookup_size)
{
local_results[k] += 2*sinc_lookup[(int)qr] * double(n_stacks - l);
}
else
{
local_results[k] += 2* boost::math::sinc_pi(double(qr) / 1000)* double(n_stacks - l);
}
}// k index
}//j index loop
}//l index loop
//Notify user to current status


}// i index loop

//Lock mutex
boost::mutex::scoped_lock lock(io_mutex);

//Copy data from local temporary vector to results
for (unsigned int i = 0; i < results.size(); ++i)
{
results[i]->m_I += local_results[i];
}

//Unlock mutex
lock.unlock();
}
*/