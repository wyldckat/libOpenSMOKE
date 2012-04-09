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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


#ifndef OpenSMOKE_SootEmpiricalModels_H
#define OpenSMOKE_SootEmpiricalModels_H

#include "OpenSMOKE_External_Functions.hpp"

class OpenSMOKE_SootEmpiricalModels
{
public:

	OpenSMOKE_SootEmpiricalModels(	const nucleation_models nucleation_model, const growth_models growth_model,
					const oxidation_models oxidation_model, const aggregation_models aggregation_model	);

	OpenSMOKE_SootEmpiricalModels(const std::vector<std::string> names, const std::vector<int> indices, 
								  const nucleation_models nucleation_model, const growth_models growth_model,
								  const oxidation_models oxidation_model, const aggregation_models aggregation_model);

	void update(const double &T, const double &p_atm, const double &rho, const std::vector<double> &omega, 
				const double &phiN, const double &phiM, double &s, double &S, std::vector<double> &Sgas);

	void update_properties(	const double &rho,const double &phiN, const double &phiM,
				double &m0, double &fv, double &dSoot, double &MSoot, double &ASoot);

	void setRobustCalculations();

	double phiNStart(const double &rho);
	double phiMStart(const double &rho);

	~OpenSMOKE_SootEmpiricalModels(void);



private:

	void setup(	const std::vector<std::string> names, const std::vector<int> indices, 
			const nucleation_models nucleation_model, const growth_models growth_model,
			const oxidation_models oxidation_model, const aggregation_models aggregation_model);

	nucleation_models	nucleation_model_;
	growth_models		growth_model_;
	aggregation_models	aggregation_model_;
	oxidation_models	oxidation_model_;
	surface_functions	surface_function_;

	double n_nuclei_;			// number of carbon atom nuclei in primary particles [-]
	double dp_;					// diameter of primary particles [m]
	double mp_;					// primary particle molecular weight [kg/kmol]
	double A_soot_;				// surface area per unit volume [1/m]
	double A_soot_corrected_;	// surface area per unit volume [1/m]
	double m0_;					// number of particle per unit volume [#/m3]
	double fv_;					// volume fraction [-]

	double phiNStart_;		// initial value of phiN
	double phiMStart_;		// intitial value of phiM

	double A_nucleation_;	// nucleation frequency factor [1/s]
	double A_growth_;		// growth frequency factor []
	double T_nucleation_;	// nucleation activation energy [K]
	double T_growth_;		// growth activation energy [K]

	double Sn_;				// nucleation soot term [kmol/m3/s]
	double SN_;				// nucleation soot term [kg/m3/s]
	double SG_;				// growth source term [kg/m3/s]
	double SO_oh_;			// oxidation soot term [kg/m3/s]
	double SO_o2_;			// oxidation soot term [kg/m3/s]
	double Sc_;				// aggregation soot term [kmol/m3/s]

	double SNGas_c2h2_;		// [kg/m3/s]
	double SGGas_c2h2_;		// [kg/m3/s]
	double SNGas_h2_;		// [kg/m3/s]
	double SGGas_h2_;		// [kg/m3/s]
	double SOGas_oh_;		// [kg/m3/s]
	double SOGas_o2_;		// [kg/m3/s]
	double SOGas_co_;		// [kg/m3/s]
	double SOGas_h_;		// [kg/m3/s]


	unsigned int i_h2;
	unsigned int i_o2;
	unsigned int i_c2h2;
	unsigned int i_oh;
	unsigned int i_co;
	unsigned int i_h;

	bool _iRobustCalculations;

	void nucleation_rate(const double &T, const double &c_c2h2);
	void growth_rate(const double &T, const double &c_c2h2);
	void oxidation_rate(const double &T, const double &c_o2, const double &c_oh, const double &p_o2);
	void aggregation_rate(const double &T);

	FractionalExponents	*gamma_fv;
	FractionalExponents	*gamma_m0;

private:

	void Summary();
	void ErrorMessage(const std::string error_message);
	void WarningMessage(const std::string warning_message);
};


#endif // OpenSMOKE_SootEmpiricalModels_H
