/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
 *   alberto.cuoci@polimi.it   						                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "stdlib.h"
#include <vector>
#include <string>

#ifndef OpenSMOKE_QMOM_EmpiricalModels_H
#define OpenSMOKE_QMOM_EmpiricalModels_H

#include "OpenSMOKE_External_Functions.hpp"
#include "OpenSMOKE_QMOM_GordonAlgorithm.h"

class OpenSMOKE_QMOM_EmpiricalModels
{
public:

	OpenSMOKE_QMOM_EmpiricalModels(	const int N, const nucleation_models nucleation_model, const growth_models growth_model,
									const oxidation_models oxidation_model, const aggregation_models aggregation_model	);

	OpenSMOKE_QMOM_EmpiricalModels( const int N, const std::vector<std::string> names, const std::vector<int> indices, 
								    const nucleation_models nucleation_model, const growth_models growth_model,
								    const oxidation_models oxidation_model, const aggregation_models aggregation_model);

	void update(const double &T, const double &p_atm, const double &rho, const double &mu, const std::vector<double> &omega, 
				double *moments, double *S);

	void update_properties(double &m0, double &fv, double &dSoot, double &MSoot, double &ASoot);


	double* moments0();

//	void update(const double &T, const double &p_atm, const double &rho, const std::std::vector<double> &omega, 
//				const double &phiN, const double &phiM, double &s, double &S, std::std::vector<double> &Sgas);
//
//	void update_properties(	const double &rho,const double &phiN, const double &phiM,
//				double &m0, double &fv, double &dSoot, double &MSoot, double &ASoot);


	~OpenSMOKE_QMOM_EmpiricalModels(void);

private:

	void setup(	const int N, const std::vector<std::string> names, const std::vector<int> indices, 
			    const nucleation_models nucleation_model, const growth_models growth_model,
			    const oxidation_models oxidation_model, const aggregation_models aggregation_model);

	double n_nuclei_;		// number of carbon atom nuclei in primary particles [-]
	double dp_;				// diameter of primary particles [m]
	double mp_;				// primary particle molecular weight [kg/kmol]
	double A_soot_;			// surface area per unit volume [1/m]
	double m0_;				// number of particle per unit volume [#/m3]
	double fv_;				// volume fraction [-]
	double M_soot_;			// soot mass  
	double d_soot_;			// soot mass  

	// Kinetics parameters
	double A_nucleation_;	// nucleation frequency factor [1/s]
	double A_growth_;		// growth frequency factor []
	double T_nucleation_;	// nucleation activation energy [K]
	double T_growth_;		// growth activation energy [K]

	// Internal functions
	void update_weights_and_abscissas();
	void update_aggregation_kernel(const double &mu, const double &P_atm, const double &T, const double &MW);
	void nucleation_rate(const double &T, const double &c_c2h2);
	void growth_rate(const double &T, const double &c_c2h2);
	void oxidation_rate(const double &T, const double &c_o2, const double &c_oh, const double &p_o2);
	void aggregation_rate(const double &mu, const double &P_atm, const double &T, const double &MW);

	// Species indices
	unsigned int i_h2;
	unsigned int i_o2;
	unsigned int i_c2h2;
	unsigned int i_oh;
	unsigned int i_co;
	unsigned int i_h;

	// Models
	nucleation_models	nucleation_model_;
	growth_models		growth_model_;
	aggregation_models	aggregation_model_;
	oxidation_models	oxidation_model_;
	surface_functions	surface_function_;

	double aggregation_kernel(	const double Rc, const double D, const double c, const double g);
	double aggregation_kernel(	const double Rc1, const double D1, const double c1, const double g1,
								const double Rc2, const double D2, const double c2, const double g2);

private:

	OpenSMOKE_QMOM_GordonAlgorithm *gordon;

	int N_;						// number of delta dirac
	double epsilon_;			// maximum size of primary particles
	double *weights_;			// local weights
	double *abscissas_;			// local abscissas
	double *momentsN_;			// normalized moments
	double *moments_;			// moments
	double *moments0_;			// initial distribution moments
	double *pow_unit;
	
	double *Snucleation_;		// nucleation sources
	double *Sgrowth_;			// growth sources
	double *Soxidation_;		// oxidation sources
	double *Saggregation_;		// aggregation terms

	double *L3;
	double *D;
	double *c;
	double *g;
	double *Rc;
	double *Beta;

	// Units and minimum abscissa
	double unit;
	double abscissa_min;

private:

	void Summary();
	void ErrorMessage(const std::string error_message);
	void WarningMessage(const std::string warning_message);
};


#endif // OpenSMOKE_QMOM_EmpiricalModels_H
