/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci		     				       *
 *   alberto.cuoci@polimi.it                                               *
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

#include "OpenSMOKE_SootEmpiricalModels.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>


using namespace std;


OpenSMOKE_SootEmpiricalModels::OpenSMOKE_SootEmpiricalModels(	const nucleation_models nucleation_model, const growth_models growth_model,
								const oxidation_models oxidation_model, const aggregation_models aggregation_model	)
{
	std::vector<std::string> list_names(6);
	std::vector<int> list_indices(6);

	list_names[0] = "C2H2"; list_indices[0]=0;
	list_names[1] = "H2"; 	list_indices[1]=1;
	list_names[2] = "O2"; 	list_indices[2]=2;
	list_names[3] = "OH"; 	list_indices[3]=3;
	list_names[4] = "CO"; 	list_indices[4]=4;
	list_names[5] = "H"; 	list_indices[5]=5;

	setup(list_names, list_indices, nucleation_model, growth_model, oxidation_model, aggregation_model);

	_iRobustCalculations = false;
	gamma_fv = new FractionalExponents(1.e-10, 1.e11,  1.e-5);
	gamma_m0 = new FractionalExponents(1.e15,  1.e-14, 1.e-1);
}

OpenSMOKE_SootEmpiricalModels::OpenSMOKE_SootEmpiricalModels(	const vector<string> names, const vector<int> indices, 
								const nucleation_models nucleation_model, const growth_models growth_model,
								const oxidation_models oxidation_model, const aggregation_models aggregation_model
								)
{
	setup(names, indices, nucleation_model, growth_model, oxidation_model, aggregation_model);
}

void OpenSMOKE_SootEmpiricalModels::setup(	const vector<string> names, const vector<int> indices, 
						const nucleation_models nucleation_model, const growth_models growth_model,
						const oxidation_models oxidation_model, const aggregation_models aggregation_model
					 )
{
	nucleation_model_	= nucleation_model;
	growth_model_		= growth_model;
	oxidation_model_	= oxidation_model;
	aggregation_model_	= aggregation_model;

	switch(nucleation_model_)
	{
		case NUCLEATION_LIU_2001:		// Liu 2001

			n_nuclei_		= Liu_2001::n_nuclei_;
			A_nucleation_	= Liu_2001::A_nucleation_;
			T_nucleation_	= Liu_2001::T_nucleation_;
				
			break;
		
		case NUCLEATION_LIU_2002:		// Liu 2002

			n_nuclei_		= Liu_2002::n_nuclei_;
			A_nucleation_	= Liu_2002::A_nucleation_;
			T_nucleation_	= Liu_2002::T_nucleation_;

			break;

		case NUCLEATION_MOSS_1999:		// Moss 1999

			n_nuclei_		= Moss_1999::n_nuclei_;
			A_nucleation_	= Moss_1999::A_nucleation_;
			T_nucleation_	= Moss_1999::T_nucleation_;

			break;

		case NUCLEATION_WEN_2003:		// Wen 2003

			n_nuclei_		= Wen_2003::n_nuclei_;
			A_nucleation_	= Wen_2003::A_nucleation_;
			T_nucleation_	= Wen_2003::T_nucleation_;

			break;

		case NUCLEATION_LINDSTEDT_1994:	// Lindstedt 1994

			n_nuclei_		= Lindstedt_1994::n_nuclei_;
			A_nucleation_	= Lindstedt_1994::A_nucleation_;
			T_nucleation_	= Lindstedt_1994::T_nucleation_;

			break;

		case NUCLEATION_LEUNG_1991:		// Leung 1991

			n_nuclei_		= Leung_1991::n_nuclei_;
			A_nucleation_	= Leung_1991::A_nucleation_;
			T_nucleation_	= Leung_1991::T_nucleation_;

			break;

		case NUCLEATION_NONE:

			n_nuclei_ = 0.;
			
			break;
 
	}

	switch(growth_model_)
	{

		case GROWTH_LIU_2001:		// Liu 2001

			A_growth_			= Liu_2001::A_growth_;
			T_growth_			= Liu_2001::T_growth_;
			surface_function_	= Liu_2001::surface_function_;

			break;
		
		case GROWTH_LIU_2002:		// Liu 2002,2006

			A_growth_			= Liu_2002::A_growth_;
			T_growth_			= Liu_2002::T_growth_;
			surface_function_	= Liu_2002::surface_function_;
		
			break;


		case GROWTH_MOSS_1999:		// Moss 1999

			A_growth_			= Moss_1999::A_growth_;
			T_growth_			= Moss_1999::T_growth_;
			surface_function_	= Moss_1999::surface_function_;

			break;

		case GROWTH_WEN_2003:		// Wen 2003

			A_growth_			= Wen_2003::A_growth_;
			T_growth_			= Wen_2003::T_growth_;
			surface_function_	= Wen_2003::surface_function_;
		
			break;

		case GROWTH_LINDSTEDT_1994:		// Lindstedt 1994

			A_growth_			= Lindstedt_1994::A_growth_;
			T_growth_			= Lindstedt_1994::T_growth_;
			surface_function_	= Lindstedt_1994::surface_function_;

			break;

		case GROWTH_LEUNG_1991:			// Leung 1991

			A_growth_			= Leung_1991::A_growth_;
			T_growth_			= Leung_1991::T_growth_;
			surface_function_	= Leung_1991::surface_function_;

			break;

		case GROWTH_NONE:

			A_growth_ = 0.;
			T_growth_ = 0.;

			surface_function_ = SURFACE_FUNCTION_LINEAR;
			
			break;
	}

	// Soot primary particles
	dp_ = 0.27029949e-9*pow(n_nuclei_, 1./3.);		// [m]
	mp_ = OpenSMOKE::mw_c * n_nuclei_;				// [kg/kmol]

	// Indices
	unsigned int count = 0;
	for(unsigned int j=0;j<names.size();j++)
	{
		if (!names[j].compare("C2H2"))	{	i_c2h2 = indices[j];	count++; }
		if (!names[j].compare("H2"))	{	i_h2   = indices[j];	count++; }
		if (!names[j].compare("O2"))	{	i_o2   = indices[j];	count++; }
		if (!names[j].compare("OH"))	{	i_oh   = indices[j];	count++; }
		if (!names[j].compare("CO"))	{	i_co   = indices[j];	count++; }
		if (!names[j].compare("H"))		{	i_h	   = indices[j];	count++; }
	}
	if (count != 6)
		ErrorMessage("Wrong list of gas species...");

	// Initialize source terms
	Sn_		= 0.;		// [kmol/m3/s]
	SN_		= 0.;		// [kg/m3/s]
	SG_		= 0.;		// [kg/m3/s]	
	SO_o2_	= 0.;		// [kg/m3/s]
	SO_oh_	= 0.;		// [kg/m3/s]
	Sc_		= 0.;		// [kmol/m3/s]
	SNGas_c2h2_	= 0.;	// [kg/m3/s] Acetylene consumption
	SNGas_h2_	= 0.;	// [kg/m3/s] Hydrogen formation
	SGGas_c2h2_	= 0.;	// [kg/m3/s] Acetylene consumption
	SGGas_h2_	= 0.;	// [kg/m3/s] Hydrogen formation
	SOGas_o2_	= 0.;	// [kg/m3/s]	O2 consumption
	SOGas_oh_	= 0.;	// [kg/m3/s]	OH consumption
	SOGas_co_	= 0.;	// [kg/m3/s]	CO formation
	SOGas_h_	= 0.;	// [kg/m3/s]	H  formation

	Summary();
}

void OpenSMOKE_SootEmpiricalModels::setRobustCalculations()
{
	_iRobustCalculations = true;
}

double OpenSMOKE_SootEmpiricalModels::phiNStart(const double &rho)
{
	double fvStart = 1.e-14;
	double m0Start = (6./OpenSMOKE::pi/pow(dp_,3.))*fvStart;
	
	return m0Start/(rho*OpenSMOKE::Nav_kmol);
}

double OpenSMOKE_SootEmpiricalModels::phiMStart(const double &rho)
{
	double fvStart = 1.e-14;
	
	return  OpenSMOKE::rho_soot*fvStart/rho;
}

void OpenSMOKE_SootEmpiricalModels::update(const double &T, const double &p_atm, const double &rho, const vector<double> &omega, const double &phiN, const double &phiM, double &s, double &S, vector<double> &SGas)
{
	double c_tot  = p_atm*101325./(8314.*T);	// total concentration [kmol/m3]
	double mw_mix = rho/c_tot;					// total molecular weight [kg/kmol]

	m0_ = rho*OpenSMOKE::Nav_kmol*phiN;			// soot particle number density [#/m3]
	fv_ = rho/OpenSMOKE::rho_soot*phiM;			// volume fraction [-]
	
	// Soot specific area
	A_soot_ = OpenSMOKE::_36pi_to_1_over_3 * pow(m0_, 1./3.) * pow(fv_, 2./3.);	// soot specific area [1/m]

	// Correction of soot specific area (for oxidation)
	if (_iRobustCalculations == true)
		A_soot_corrected_ = OpenSMOKE::_36pi_to_1_over_3 * gamma_m0->gamma(m0_, 1./3.) * gamma_fv->gamma(fv_, 2./3.);
	else
		A_soot_corrected_ = A_soot_;

	// Concentrations
	double p_o2		= omega[i_o2]*mw_mix/OpenSMOKE::mw_o2 *p_atm;			// oxygen partial pressure [atm]
	double c_c2h2	= omega[i_c2h2]*mw_mix/OpenSMOKE::mw_c2h2 * c_tot;		// acetylene concentration [kmol/m3]
	double c_o2		= omega[i_o2]*mw_mix/OpenSMOKE::mw_o2 * c_tot;			// oxygen concentration [kmol/m3]
	double c_oh		= omega[i_oh]*mw_mix/OpenSMOKE::mw_oh * c_tot;			// OH concentration [kmol/m3]


	// Soot formation/consumption rates
	if (nucleation_model_	!= NUCLEATION_NONE)		nucleation_rate(T, c_c2h2);
	if (growth_model_		!= GROWTH_NONE)			growth_rate(T, c_c2h2);
	if (oxidation_model_	!= OXIDATION_NONE)		oxidation_rate(T, c_o2, c_oh, p_o2);
	if (aggregation_model_	!= AGGREGATION_NONE)	aggregation_rate(T);

	// Soot source terms
	s = Sn_ - Sc_;						// [kmol/m3/s]
	S = SN_ + SG_ - (SO_o2_ + SO_oh_);	// [kg/m3/s]

	// Gas phase
	SGas[i_c2h2] = 	SNGas_c2h2_ + SGGas_c2h2_;	// [kg/m3/s] Acetylene consumption
	SGas[i_h2]	= 	SNGas_h2_ + SGGas_h2_;		// [kg/m3/s] Hydrogen formation
	SGas[i_oh]	= 	SOGas_oh_;					// [kg/m3/s] OH consumption
	SGas[i_co]	= 	SOGas_co_;					// [kg/m3/s] CO formation
	SGas[i_o2]	= 	SOGas_o2_;					// [kg/m3/s] O2 consumption
	SGas[i_h]	= 	SOGas_h_;					// [kg/m3/s] H formation

}

void OpenSMOKE_SootEmpiricalModels::update_properties(const double &rho,const double &phiN, const double &phiM,
						      double &m0, double &fv, double &dSoot, double &MSoot, double &ASoot)
{
	m0	= rho*OpenSMOKE::Nav_kmol*phiN;							// soot particle number density [#/m3]
	fv	= rho/OpenSMOKE::rho_soot*phiM;							// volume fraction [-]
	dSoot	= pow(6.*fv/(OpenSMOKE::pi*m0), 1./3.);						// soot particle diameter [m]
	MSoot	= OpenSMOKE::rho_soot * fv;							// soot density [kg/m3]
	ASoot   = OpenSMOKE::_36pi_to_1_over_3 * pow(m0, 1./3.) * pow(fv, 2./3.);			// soot specific area [1/m]
}

void OpenSMOKE_SootEmpiricalModels::nucleation_rate(const double &T, const double &c_c2h2)
{
	Sn_ = A_nucleation_ * exp(-T_nucleation_/T) * c_c2h2;	// [kmol/m3/s]
	SN_ = Sn_ * mp_;										// [kg/m3/s]

	SNGas_c2h2_	= - SN_/(2.*OpenSMOKE::mw_c) * OpenSMOKE::mw_c2h2;				// [kg/m3/s] Acetylene consumption
	SNGas_h2_	=   SN_/(2.*OpenSMOKE::mw_c) * OpenSMOKE::mw_h2;				// [kg/m3/s] Hydrogen formation
}

void OpenSMOKE_SootEmpiricalModels::growth_rate(const double &T, const double &c_c2h2)
{
	if (surface_function_ == SURFACE_FUNCTION_LINEAR)
		SG_ = A_growth_ * exp(-T_growth_/T) * A_soot_* c_c2h2;			// [kg/m3/s]	
	else
		SG_ = A_growth_ * exp(-T_growth_/T) * sqrt(A_soot_)* c_c2h2;	// [kg/m3/s]	

	SGGas_c2h2_ = -SG_ / (2.*OpenSMOKE::mw_c) * OpenSMOKE::mw_c2h2;		// [kg/m3/s] Acetylene consumption
	SGGas_h2_   =  SG_ / (2.*OpenSMOKE::mw_c) * OpenSMOKE::mw_h2;		// [kg/m3/s] Hydrogen formation
}

void OpenSMOKE_SootEmpiricalModels::oxidation_rate(const double &T, const double &c_o2, const double &c_oh, const double &p_o2)
{
	switch(oxidation_model_)
	{
		case OXIDATION_LEE:		// Lee 1962
		
			SO_o2_ = Lee::A * exp(-Lee::T/T) * sqrt(T) * c_o2;		// [kg/m2/s]

			SO_oh_ = 0.;											// [kg/m3/s]

			break;
		
		case OXIDATION_NEOH:	// Neoh 1980

			SO_oh_ = Neoh::eta*Neoh::A * sqrt(T) * c_oh;			// [kg/m2/s]

			SO_o2_ = 0.;											// [kg/m3/s]
		
			break;

		case OXIDATION_NSC:		// NSC 1966

			double kA, kB, kT, kZ, chi;

			if (p_o2>1.e-12)
			{
				kA = NSC::A	* exp(-NSC::TA/T);	// [g/cm2/s/atm]
				kB = NSC::B * exp(-NSC::TB/T);	// [g/cm2/s/atm]
				kT = NSC::C * exp(-NSC::TC/T);	// [g/cm2/s]
				kZ = NSC::D * exp(+NSC::TD/T);	// [1/atm]

				chi = 1./(1.+kT/(kB*p_o2));										// [-]

				SO_o2_  = 12. * (kA*p_o2/(1.+kZ*p_o2)*chi + kB*p_o2*(1.-chi));	// [g/cm2/s]
				SO_o2_ *= 10.;													// [kg/m2/s]

				SO_oh_  = 0.;													// [kg/m3/s]
			}
			else
			{
				SO_o2_ = 0.;
				SO_oh_ = 0.;
			}

			break;


		case OXIDATION_NONE:

			SO_o2_ = 0.;			// [kg/m3/s]
			SO_oh_ = 0.;			// [kg/m3/s]
			
			break;

	}

	// Soot source terms
	SO_o2_ *= A_soot_corrected_;
	SO_oh_ *= A_soot_corrected_;

	// Gas phase
	SOGas_o2_ = - SO_o2_/OpenSMOKE::mw_c * 0.50 * OpenSMOKE::mw_o2;			// [kg/m3/s]	O2 consumption
	SOGas_oh_ = - SO_oh_/OpenSMOKE::mw_c * 1.00 * OpenSMOKE::mw_oh;			// [kg/m3/s]	OH consumption
	SOGas_co_ =   (SO_o2_+SO_oh_)/OpenSMOKE::mw_c * OpenSMOKE::mw_co;		// [kg/m3/s]	CO formation
	SOGas_h_  =   SO_oh_/OpenSMOKE::mw_c * OpenSMOKE::mw_h;					// [kg/m3/s]	H  formation
}

void OpenSMOKE_SootEmpiricalModels::aggregation_rate(const double &T)
{
	switch(aggregation_model_)
	{

		case AGGREGATION_SMOLUCHOWSKI:		// Smoluchowski
			
			if (_iRobustCalculations == true)
				Sc_ = Smoluchowski::gamma * sqrt(T) * gamma_fv->gamma(fv_, 0.166667)* pow(m0_, 11./6.);	// [kmol/m3/s]
			else
				Sc_ = Smoluchowski::gamma * sqrt(T) * pow(fv_, 0.16667)* pow(m0_, 11./6.);				// [kmol/m3/s]

			break;

		case AGGREGATION_MOSS:							// Moss

			Sc_ = Moss::gamma * sqrt(T) * (m0_*m0_);	// [kmol/m3/s]

			break;

		case AGGREGATION_NONE:

			Sc_ = 0.;
			
			break;
	}
}

OpenSMOKE_SootEmpiricalModels::~OpenSMOKE_SootEmpiricalModels(void)
{
}

void OpenSMOKE_SootEmpiricalModels::Summary()
{
	cout << endl;
	cout << "OpenSMOKE_SootEmpiricalModels      " << endl;
	cout << "Nucleation:                        " << nucleation_model_	<< endl;
	cout << "Growth:                            " << growth_model_		<< endl;
	cout << "Oxidation:                         " << oxidation_model_	<< endl;
	cout << "Aggregation:                       " << aggregation_model_ << endl;
	cout << "Diameter (primary particles):      " << dp_*1e9 << " nm" << endl;
	cout << "Mass (primary particles):          " << mp_ << " kg/kmol" << endl;
	cout << "Nuclei (primary particles):        " << n_nuclei_ << endl;
	cout << "PhiM(@1kg/m3):                     " << phiMStart(1.0) << endl;
	cout << "PhiN(@1kg/m3):                     " << phiNStart(1.0) << endl;
	cout << endl;
}

void OpenSMOKE_SootEmpiricalModels::ErrorMessage(const string error_message)
{
	cout << "OpenSMOKE_SootEmpiricalModels   " << endl;
	cout << "Error message: " << error_message << endl;
	exit(-1);
}

void OpenSMOKE_SootEmpiricalModels::WarningMessage(const string warning_message)
{
	cout << "OpenSMOKE_SootEmpiricalModels       " << endl;
	cout << "Warning message: " << warning_message << endl;
	cout << "Press enter to continue..." << endl;
	getchar();
}
