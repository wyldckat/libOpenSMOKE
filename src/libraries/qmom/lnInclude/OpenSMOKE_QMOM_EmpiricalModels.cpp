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

#include "OpenSMOKE_QMOM_EmpiricalModels.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;

OpenSMOKE_QMOM_EmpiricalModels::OpenSMOKE_QMOM_EmpiricalModels(	const int N, const nucleation_models nucleation_model, const growth_models growth_model,
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

	setup(N, list_names, list_indices, nucleation_model, growth_model, oxidation_model, aggregation_model);
}

OpenSMOKE_QMOM_EmpiricalModels::OpenSMOKE_QMOM_EmpiricalModels(	const int N, const vector<string> names, const vector<int> indices, 
																const nucleation_models nucleation_model, const growth_models growth_model,
																const oxidation_models oxidation_model, const aggregation_models aggregation_model )
{
	setup(N, names, indices, nucleation_model, growth_model, oxidation_model, aggregation_model);
}

void OpenSMOKE_QMOM_EmpiricalModels::setup(	const int N, const vector<string> names, const vector<int> indices, 
											const nucleation_models nucleation_model, const growth_models growth_model,
											const oxidation_models oxidation_model, const aggregation_models aggregation_model
											)
{
	N_ = N;

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
	dp_ = 0.27029949e-9*pow(n_nuclei_, 1./3.);			// [m]
	mp_ = OpenSMOKE::mw_c * n_nuclei_;					// [kg/kmol]
	epsilon_ = pow(4.,1./3.)*dp_;						// [m]
//	epsilon_ = 4./3.*dp_;								// [m]	dp = m3/m2
	
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

    // Memory allocation
	gordon		= new OpenSMOKE_QMOM_GordonAlgorithm(N_);
	weights_	= new double[N_];
	abscissas_	= new double[N_];
	
	moments_		= new double[2*N_];
	momentsN_		= new double[2*N_];
	moments0_		= new double[2*N_];
	Snucleation_	= new double[2*N_];
	Sgrowth_		= new double[2*N_];
	Soxidation_		= new double[2*N_];
	Saggregation_	= new double[2*N_];

	L3		= new double[N_];
	D		= new double[N_];
	c		= new double[N_];
	g		= new double[N_];
	Rc		= new double[N_];
	Beta	= new double[N_*N_];

	for(int k=0;k<2*N_;k++)
	{	
		Snucleation_[k]		= 0.;
		Sgrowth_[k]			= 0.;
		Soxidation_[k]		= 0.;
		Saggregation_[k]	= 0.;
	}

	// Initial distribution
	{
		double fv0_ = 1e-10;										// [-]
		double n0_ = 24./OpenSMOKE::pi * fv0_ / pow(epsilon_, 4.);	// [#/m4]

		for(int j=0;j<2*N_;j++)
			moments0_[j] = n0_ * pow(epsilon_, j+1)/double(j+1.);

		for(int j=0;j<2*N_;j++)
			momentsN_[j] = moments0_[j]/moments0_[0];

		gordon->Apply(momentsN_, abscissas_, weights_);

		for(int j=0;j<N_;j++)
			weights_[j] *= moments0_[0];
	}

	Summary();

	// Units and minimum abscissa
	unit = 1e-9;
	abscissa_min = dp_/1.e0;
	pow_unit = new double[2*N_];
	for(int j=0;j<2*N_;j++)
		pow_unit[j] = pow(unit,j);
}

void OpenSMOKE_QMOM_EmpiricalModels::update(const double &T, const double &p_atm, const double &rho, const double &mu, const vector<double> &omega, double *moments, double *S)
{
	// Local copy for additional calculations
	for(int j=0;j<2*N_;j++)
		moments_[j] = moments[j];

	// From unit momemts to real moments
	for(int j=0;j<2*N_;j++)
		moments_[j] *= moments0_[j];

	// Update weights and abscissas
	if (aggregation_model_ != AGGREGATION_NONE)
		update_weights_and_abscissas();
	
	// Particle number density
	m0_ = moments_[0];

	// Soot volume fraction
	fv_ = OpenSMOKE::pi/6.*moments_[3];

	// Soot specific area
	A_soot_ = OpenSMOKE::pi*moments_[2];
	A_soot_ = OpenSMOKE::_36pi_to_1_over_3 * pow(m0_, 1./3.) * pow(fv_, 2./3.);	

	// Soot specific area
	M_soot_ = OpenSMOKE::rho_soot * fv_;

	// Soot particle mean diameter
	d_soot_ = moments_[3]/moments_[2];

	// Gas properties
	double c_tot  = p_atm*101325./(8314.*T);	// total concentration [kmol/m3]
	double mw_mix = rho/c_tot;					// total molecular weight [kg/kmol]	
	
	// Concentrations
	double p_o2		= omega[i_o2]*mw_mix/OpenSMOKE::mw_o2 *p_atm;			// oxygen partial pressure [atm]
	double c_c2h2	= omega[i_c2h2]*mw_mix/OpenSMOKE::mw_c2h2 * c_tot;		// acetylene concentration [kmol/m3]
	double c_o2		= omega[i_o2]*mw_mix/OpenSMOKE::mw_o2 * c_tot;			// oxygen concentration [kmol/m3]
	double c_oh		= omega[i_oh]*mw_mix/OpenSMOKE::mw_oh * c_tot;			// OH concentration [kmol/m3]

	// Soot formation/consumption rates
	if (nucleation_model_	!= NUCLEATION_NONE)		nucleation_rate(T, c_c2h2);
	if (growth_model_		!= GROWTH_NONE)			growth_rate(T, c_c2h2);
	if (oxidation_model_	!= OXIDATION_NONE)		oxidation_rate(T, c_o2, c_oh, p_o2);
	if (aggregation_model_	!= AGGREGATION_NONE)	aggregation_rate(mu, p_atm, T, mw_mix);

	// Soot source terms
	for(int k=0;k<2*N_;k++)
		S[k] = Snucleation_[k] + Sgrowth_[k] - Soxidation_[k] + Saggregation_[k];

	// From real sources to unit sources
	for(int k=0;k<2*N_;k++)
		S[k] /= moments0_[k];
}

void OpenSMOKE_QMOM_EmpiricalModels::update_properties(double &m0, double &fv, double &dSoot, double &MSoot, double &ASoot)
{
	m0    = m0_;			//
	fv    = fv_;			//
	dSoot = d_soot_;		// soot particle diameter [m]
	MSoot = M_soot_;		// soot density [kg/m3]
	ASoot = A_soot_;		//
}

void OpenSMOKE_QMOM_EmpiricalModels::update_weights_and_abscissas()
{
	// Normalize moments
	for(int j=0;j<2*N_;j++)
		momentsN_[j] = moments_[j]/moments_[0];

	// Change unit lenghts
	for(int j=0;j<2*N_;j++)
		momentsN_[j] /= pow_unit[j];

	// Weights and abscissas (normalized)
	gordon->Apply(momentsN_, abscissas_, weights_);

	// Weights
	for(int j=0;j<N_;j++)
		weights_[j] *= moments_[0];

	// Absissas
	for(int j=0;j<N_;j++)
		abscissas_[j] *= unit;

	// Abscissas (check)
	for(int j=0;j<N_;j++)
		abscissas_[j] = (abscissas_[j] < abscissa_min) ?  abscissa_min : abscissas_[j];
}

void OpenSMOKE_QMOM_EmpiricalModels::nucleation_rate(const double &T, const double &c_c2h2)
{
	double Sn_ = A_nucleation_ * exp(-T_nucleation_/T) * c_c2h2;	// [kmol/m3/s]
	double J   = Sn_*OpenSMOKE::Nav_kmol;							// [#/m3/s]
	
	// Sources
	for(int k=0;k<2*N_;k++)
		Snucleation_[k] = pow(epsilon_,k)/double(k+1.) * J;			// 
}

void OpenSMOKE_QMOM_EmpiricalModels::growth_rate(const double &T, const double &c_c2h2)
{
	double SG_specific_ = A_growth_ * exp(-T_growth_/T) * c_c2h2;		// [kg/m2/s]

	double G;
	if (surface_function_ == SURFACE_FUNCTION_LINEAR)
		G = 2.*SG_specific_/OpenSMOKE::rho_soot;					// 
	else
		G = 2.*SG_specific_/OpenSMOKE::rho_soot/sqrt(A_soot_);	// [kg/m3/s]	

	// Sources
	Sgrowth_[0] = 0.;
	for(int k=1;k<2*N_;k++)
		Sgrowth_[k] = double(k)*moments_[k-1]*G;
}

void OpenSMOKE_QMOM_EmpiricalModels::oxidation_rate(const double &T, const double &c_o2, const double &c_oh, const double &p_o2)
{
	double SO_specific = 0.;

	switch(oxidation_model_)
	{
		case OXIDATION_LEE:		// Lee 1962
		
			SO_specific = Lee::A * exp(-Lee::T/T) * sqrt(T) * c_o2;		// [kg/m2/s]

			break;
		
		case OXIDATION_NEOH:	// Neoh 1980

			SO_specific = Neoh::eta*Neoh::A * sqrt(T) * c_oh;			// [kg/m2/s]

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

				SO_specific  = 12. * (kA*p_o2/(1.+kZ*p_o2)*chi + kB*p_o2*(1.-chi));	// [g/cm2/s]
				SO_specific *= 10.;													// [kg/m2/s]
			}

			break;

		case OXIDATION_NONE:

			SO_specific = 0.;			// [kg/m2/s]
			
			break;
	}

	double G = 2.*SO_specific/OpenSMOKE::rho_soot;					// 

	// Sources
	Soxidation_[0] = 0.;
	for(int k=1;k<2*N_;k++)
		Soxidation_[k] = double(k)*moments_[k-1]*G;
}

void OpenSMOKE_QMOM_EmpiricalModels::aggregation_rate(const double &mu, const double &P_atm, const double &T, const double &MW)
{
	update_aggregation_kernel(mu, P_atm, T, MW);

	if (N_>2)
	{
		for(int k=0;k<2*N_;k++)
		{
			double sum_a = 0.;
			for(int j=0;j<N_;j++)
			{
				double sum = 0.;
				for(int i=0;i<N_;i++)
					sum += pow( L3[j]+L3[i], double(k)/3.)*Beta[j*N_+i]*weights_[i];
				sum_a += sum*weights_[j];
			}

			double sum_b = 0.;
			for(int j=0;j<N_;j++)
			{
				double sum = 0.;
				for(int i=0;i<N_;i++)
					sum += Beta[j*N_+i]*weights_[i];
				sum_b += pow( abscissas_[j], k)*weights_[j]*sum;
			}

			Saggregation_[k] = 0.50*sum_a - sum_b;
		}
	}

	else
	{
		for(int k=0;k<2*N_;k++)
		{
			double sum_a =  pow(L3[0]+L3[0], double(k)/3.)* Beta[0]*weights_[0]*weights_[0] + 
					pow(L3[1]+L3[1], double(k)/3.)* Beta[3]*weights_[1]*weights_[1] +
					pow(L3[0]+L3[1], double(k)/3.)*(Beta[1]+Beta[2])*weights_[0]*weights_[1];

			double sum_b = 	weights_[0]*pow(abscissas_[0],k)*(Beta[0]*weights_[0]+Beta[1]*weights_[1])+
					weights_[1]*pow(abscissas_[1],k)*(Beta[2]*weights_[0]+Beta[3]*weights_[1]);

			Saggregation_[k] = 0.50*sum_a - sum_b;
		}
	}
}

double* OpenSMOKE_QMOM_EmpiricalModels::moments0()
{
	return moments0_;
}

OpenSMOKE_QMOM_EmpiricalModels::~OpenSMOKE_QMOM_EmpiricalModels(void)
{
}

void OpenSMOKE_QMOM_EmpiricalModels::Summary()
{
	cout << endl;
	cout << "OpenSMOKE_QMOM_EmpiricalModels     " << endl;
	cout << "Nucleation:                        " << nucleation_model_	<< endl;
	cout << "Growth:                            " << growth_model_		<< endl;
	cout << "Oxidation:                         " << oxidation_model_	<< endl;
	cout << "Aggregation:                       " << aggregation_model_ << endl;
	cout << "Diameter (primary particles):      " << dp_*1e9 << " nm" << endl;
	cout << "Mass (primary particles):          " << mp_ << " kg/kmol" << endl;
	cout << "Nuclei (primary particles):        " << n_nuclei_ << endl;
	cout << "Epsilon (primary particles):       " << epsilon_*1e9 << " nm" << endl;
	
	cout << "Initial dist. (csi and weights)    " << endl;
	for(int j=0;j<N_;j++)
		cout << " " << abscissas_[j]/epsilon_ << " " << weights_[j]/moments0_[0] << endl;
	cout << "Initial distribution (moments)  " << endl;
	for(int j=0;j<2*N_;j++)
		cout << " " << moments0_[j] << endl;

	cout << endl;
}

void OpenSMOKE_QMOM_EmpiricalModels::ErrorMessage(const string error_message)
{
	cout << "OpenSMOKE_QMOM_EmpiricalModels   " << endl;
	cout << "Error message: " << error_message << endl;
	exit(-1);
}

void OpenSMOKE_QMOM_EmpiricalModels::WarningMessage(const string warning_message)
{
	cout << "OpenSMOKE_QMOM_EmpiricalModels      " << endl;
	cout << "Warning message: " << warning_message << endl;
	cout << "Press enter to continue..." << endl;
	getchar();
}

double OpenSMOKE_QMOM_EmpiricalModels::aggregation_kernel(const double Rc, const double D, const double c, const double g)
{
	double A = Rc / (Rc + OpenSMOKE::_sqrt2_over_2*g); 
	double B = OpenSMOKE::_4_over_sqrt2 * D / (Rc*c);

	return OpenSMOKE::_16pi * D*Rc / (A+B);
}


double OpenSMOKE_QMOM_EmpiricalModels::aggregation_kernel(	const double Rc1, const double D1, const double c1, const double g1,
															const double Rc2, const double D2, const double c2, const double g2)
{
	double A = (Rc1+Rc2)/(Rc1+Rc2+sqrt(g1*g1+g2*g2));
	double B = 4.*(D1+D2)/(Rc1+Rc2)/sqrt(c1*c1+c2*c2);
	return OpenSMOKE::_4pi*(D1+D2)*(Rc1+Rc2)/(A+B);
}

void OpenSMOKE_QMOM_EmpiricalModels::update_aggregation_kernel(const double &mu, const double &P_atm, const double &T, const double &MW)
{
	// Mean free path [m]
	double lambda = free_path(mu, P_atm, T, MW);

	for(int j=0;j<N_;j++)
	{
		// Collisional radius [m]
		Rc[j] = 0.50*abscissas_[j];
			
		// Diameter^3 [m3]
		L3[j] = abscissas_[j]*abscissas_[j]*abscissas_[j];

		// Particle mass [kg]
		double mass = OpenSMOKE::rho_soot * L3[j];

		// Knudsen number [-]
		double Kn = 2.*lambda/Rc[j];	
		
		// Correction coefficient [-]
		double Cc = 1. + Kn*(1.257+0.40*exp(-1.10/Kn));		// Seinfeld, 2006
		
		// Particle diffusivity	[m2/s]
		D[j] = OpenSMOKE::kB * T /(OpenSMOKE::_6pi*mu*Rc[j]) * Cc;
		
		// Collisional velocity [m/s]
		c[j] = sqrt(OpenSMOKE::_8kB_over_pi * T/mass);

		// Modified mean free path	[m]
		double l = OpenSMOKE::_8_over_pi*D[j]/c[j];
		
		// Modified collisional diameter [m]
		g[j] = ( pow(2.*Rc[j]+l,3) - pow(4.*Rc[j]*Rc[j]+l*l,1.50) ) / (6.*Rc[j]*l) - 2.*Rc[j];

		/*
			cout << "j      " << j		<< endl;
			cout << "lambda " << lambda	<< endl;
			cout << "Rc     " << Rc[j]	<< endl;
			cout << "Kn     " << Kn		<< endl;
			cout << "Cc     " << Cc		<< endl;
			cout << "D      " << D[j]	<< endl;
			cout << "c      " << c[j]	<< endl;
			cout << "l      " << l		<< endl;
			cout << "g      " << g[j]	<< endl;
			cout << endl;
		*/
	}

	// Diagonal terms
	for(int j=0;j<N_;j++)
		Beta[N_*j+j] = aggregation_kernel(Rc[j], D[j], c[j], g[j]);

	// Extra diagonal terms
	for(int j=0;j<N_-1;j++)
		for(int i=j+1;i<N_;i++)
		{
			Beta[N_*j+i] = aggregation_kernel(Rc[j], D[j], c[j], g[j], Rc[i], D[i], c[i], g[i]);
			Beta[N_*i+j] = Beta[N_*j+i];
		} 
}

