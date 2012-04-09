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
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;

#include "OpenSMOKE_SootSourceTermLibrary.h"


OpenSMOKE_SootSourceTermLibrary::OpenSMOKE_SootSourceTermLibrary(void)
{
	_path_library = "constant/PDF-Library-Soot";

	nucleation_model_		= NUCLEATION_NONE;
	growth_model_			= GROWTH_NONE;
	aggregation_model_		= AGGREGATION_NONE;
	oxidation_model_		= OXIDATION_NONE;

	nucleation_path_		= "None";
	growth_path_			= "None";
	oxidation_path_			= "None";
	aggregation_path_		= "None";

	_iChiPDF = CHI_PDF_DIRAC;
	_extraction_mode = OPENSMOKE_SOOT_SOURCES_EXTRACTION_FLUCTUATIONS;

	_iRobustCalculations = false;
	gamma_fv = new FractionalExponents(1.e-10, 1.e11,  1.e-5);
	gamma_m0 = new FractionalExponents(1.e15,  1.e-14, 1.e-1);
}

void OpenSMOKE_SootSourceTermLibrary::setRobustCalculations()
{
	_iRobustCalculations = true;
}

void OpenSMOKE_SootSourceTermLibrary::SetLibraryPath(const string path_library)
{
	_path_library = path_library;
}

void OpenSMOKE_SootSourceTermLibrary::SetLogNormalChiDistribution()
{
	_iChiPDF = CHI_PDF_LOG_NORMAL;
}

void OpenSMOKE_SootSourceTermLibrary::SetNoFluctuationsExtractionMode()
{
	_extraction_mode = OPENSMOKE_SOOT_SOURCES_EXTRACTION_NO_FLUCTUATIONS;
}


void OpenSMOKE_SootSourceTermLibrary::SetNucleationModel(const nucleation_models nucleation_model)
{
	nucleation_model_ = nucleation_model;
	
	if (nucleation_model_ == NUCLEATION_LIU_2001)		nucleation_path_ = "Liu_2001";
	if (nucleation_model_ == NUCLEATION_LIU_2002)		nucleation_path_ = "Liu_2002";
	if (nucleation_model_ == NUCLEATION_LINDSTEDT_1994)	nucleation_path_ = "Lindstedt_1994";
	if (nucleation_model_ == NUCLEATION_LEUNG_1991)		nucleation_path_ = "Leung_1991";
	if (nucleation_model_ == NUCLEATION_MOSS_1999)		nucleation_path_ = "Moss_1999";
	if (nucleation_model_ == NUCLEATION_WEN_2003)		nucleation_path_ = "Wen_2003";

	if (nucleation_model_ == NUCLEATION_LIU_2001)		mp_ = OpenSMOKE::mw_c * Liu_2001::n_nuclei_;
	if (nucleation_model_ == NUCLEATION_LIU_2002)		mp_ = OpenSMOKE::mw_c * Liu_2002::n_nuclei_;
	if (nucleation_model_ == NUCLEATION_LINDSTEDT_1994)	mp_ = OpenSMOKE::mw_c * Lindstedt_1994::n_nuclei_;
	if (nucleation_model_ == NUCLEATION_LEUNG_1991)		mp_ = OpenSMOKE::mw_c * Leung_1991::n_nuclei_;
	if (nucleation_model_ == NUCLEATION_MOSS_1999)		mp_ = OpenSMOKE::mw_c * Moss_1999::n_nuclei_;
	if (nucleation_model_ == NUCLEATION_WEN_2003)		mp_ = OpenSMOKE::mw_c * Wen_2003::n_nuclei_;

}

void OpenSMOKE_SootSourceTermLibrary::SetGrowthModel(const growth_models growth_model)
{
	growth_model_ = growth_model;

	if (growth_model_ == GROWTH_LIU_2001)		growth_path_ = "Liu_2001";
	if (growth_model_ == GROWTH_LIU_2002)		growth_path_ = "Liu_2002";
	if (growth_model_ == GROWTH_LINDSTEDT_1994)	growth_path_ = "Lindstedt_1994";
	if (growth_model_ == GROWTH_LEUNG_1991)		growth_path_ = "Leung_1991";
	if (growth_model_ == GROWTH_MOSS_1999)		growth_path_ = "Moss_1999";
	if (growth_model_ == GROWTH_WEN_2003)		growth_path_ = "Wen_2003";

	if (growth_model_ == GROWTH_LIU_2001)		surface_function_ = Liu_2001::surface_function_;
	if (growth_model_ == GROWTH_LIU_2002)		surface_function_ = Liu_2002::surface_function_;
	if (growth_model_ == GROWTH_LINDSTEDT_1994)	surface_function_ = Lindstedt_1994::surface_function_;
	if (growth_model_ == GROWTH_LEUNG_1991)		surface_function_ = Leung_1991::surface_function_;
	if (growth_model_ == GROWTH_MOSS_1999)		surface_function_ = Moss_1999::surface_function_;
	if (growth_model_ == GROWTH_WEN_2003)		surface_function_ = Wen_2003::surface_function_;
}

void OpenSMOKE_SootSourceTermLibrary::SetAggregationModel(const aggregation_models aggregation_model)
{
	aggregation_model_ = aggregation_model;

	if (aggregation_model_ == AGGREGATION_SMOLUCHOWSKI)		aggregation_path_ = "Smoluchowski";
	if (aggregation_model_ == AGGREGATION_MOSS)				aggregation_path_ = "Moss";
}

void OpenSMOKE_SootSourceTermLibrary::SetOxidationModel(const oxidation_models oxidation_model)
{
	oxidation_model_ = oxidation_model;

	if (oxidation_model_ == OXIDATION_NSC)		oxidation_path_ = "NSC";
	if (oxidation_model_ == OXIDATION_LEE)		oxidation_path_ = "Lee";
	if (oxidation_model_ == OXIDATION_NEOH)		oxidation_path_ = "Neoh";
}

void OpenSMOKE_SootSourceTermLibrary::Setup()
{
	flamelets_library = new OpenSMOKE_SootSourceTermFlameletLibrary[4];

	for(int j=0;j<=3;j++)
	{
		if (_iChiPDF == CHI_PDF_LOG_NORMAL)	
			flamelets_library[j].SetLogNormalChiDistribution();
		if (_extraction_mode == OPENSMOKE_SOOT_SOURCES_EXTRACTION_NO_FLUCTUATIONS)	
			flamelets_library[j].SetNoFluctuationsExtractionMode();
	}

	if (nucleation_model_ != NUCLEATION_NONE)
	{
		flamelets_library[0].SetLibraryPath(_path_library + "/nucleation/" + nucleation_path_);
		flamelets_library[0].Read();
	}
	
	if (growth_model_ != GROWTH_NONE)
	{
		flamelets_library[1].SetLibraryPath(_path_library + "/growth/" + growth_path_);
		flamelets_library[1].Read();
	}

	if (oxidation_model_ != OXIDATION_NONE)
	{
		flamelets_library[2].SetLibraryPath(_path_library + "/oxidation/" + oxidation_path_);
		flamelets_library[2].Read();
	}

	if (aggregation_model_ != AGGREGATION_NONE)
	{
		flamelets_library[3].SetLibraryPath(_path_library + "/aggregation/" + aggregation_path_);
		flamelets_library[3].Read();
	}

}

void OpenSMOKE_SootSourceTermLibrary::GetMeanValues(const double &csi, const double &csiv2, const double &chi_st, 
													double &source_nucleation, double &source_growth,
													double &source_oxidation, double &source_aggregation)
{
	source_nucleation	= 0.;
	source_growth		= 0.;
	source_oxidation	= 0.;
	source_aggregation	= 0.;

	if (nucleation_model_ != NUCLEATION_NONE)
		flamelets_library[0].GetMeanValues(csi, csiv2, chi_st, source_nucleation);
	if (growth_model_ != GROWTH_NONE)
		flamelets_library[1].GetMeanValues(csi, csiv2, chi_st, source_growth);
	if (oxidation_model_ != OXIDATION_NONE)
		flamelets_library[2].GetMeanValues(csi, csiv2, chi_st, source_oxidation);
	if (aggregation_model_ != AGGREGATION_NONE)
		flamelets_library[3].GetMeanValues(csi, csiv2, chi_st, source_aggregation);
}

void OpenSMOKE_SootSourceTermLibrary::GetMeanValues(const double &csi, const double &csiv2, const double &chi_st, 
													const double &phiN, const double &phiM, const double &rho,
													double &s, double &S)
{
	double As_corrected;

	double mNucleation		= 0.;
	double nNucleation		= 0.;
	double mGrowth			= 0.;
	double nAggregation		= 0.;
	double mOxidation		= 0.;

	double m0 = phiN*rho*OpenSMOKE::Nav_kmol;		// soot particle number density
	double fv = phiM*rho/OpenSMOKE::rho_soot;		// soot volume fraction
	double As = OpenSMOKE::_36pi_to_1_over_3 * pow(m0, 1./3.) * pow(fv, 2./3.);
	
	// Correction of soot specific area (for oxidation)
	if (_iRobustCalculations == true)
		As_corrected = OpenSMOKE::_36pi_to_1_over_3 * gamma_m0->gamma(m0, 1./3.) * gamma_fv->gamma(fv, 2./3.);
	else
		As_corrected = As;


	// Get mean sources
	GetMeanValues(csi, csiv2, chi_st, mNucleation, mGrowth, mOxidation, nAggregation);


	// Nucleation
	nNucleation = mNucleation / mp_;

	// Growth
	if (surface_function_ == SURFACE_FUNCTION_SQUARE_ROOT)
		mGrowth *= sqrt(As);
	else
		mGrowth *= As;

	// Aggregation
	if (aggregation_model_ == AGGREGATION_SMOLUCHOWSKI)	
	{
		if (_iRobustCalculations == true)
			nAggregation *= Smoluchowski::gamma * gamma_fv->gamma(fv, 0.1666667)* pow(m0, 11./6.);	// [kmol/m3/s]
		else
			nAggregation *= Smoluchowski::gamma * pow(fv, 0.1666667)* pow(m0, 11./6.);				// [kmol/m3/s]
	}
	else if (aggregation_model_ == AGGREGATION_MOSS)
		nAggregation *= Moss::gamma * (m0*m0);	

	// Oxidation
	mOxidation *= As_corrected;

	// Sources
	s = nNucleation - nAggregation;
	S = mNucleation + mGrowth - mOxidation;
}

void OpenSMOKE_SootSourceTermLibrary::UpdateProperties(	const double &rho,const double &phiN, const double &phiM,
							double &m0, double &fv, double &dSoot, double &MSoot, double &ASoot)
{
	m0	= rho*OpenSMOKE::Nav_kmol*phiN;						// soot particle number density [#/m3]
	fv	= rho/OpenSMOKE::rho_soot*phiM;						// volume fraction [-]

	dSoot	= pow(6.*fv/(OpenSMOKE::pi*m0), 1./3.);									// soot particle diameter [m]
	MSoot	= OpenSMOKE::rho_soot * fv;												// soot density [kg/m3]
	ASoot   = OpenSMOKE::_36pi_to_1_over_3 * pow(m0, 1./3.) * pow(fv, 2./3.);		// soot specific area [1/m]

}

OpenSMOKE_SootSourceTermLibrary::~OpenSMOKE_SootSourceTermLibrary(void)
{
}

void OpenSMOKE_SootSourceTermLibrary::Summary()
{
	cout << endl;
	cout << "OpenSMOKE_SootSourceTermLibrary    " << endl;
	cout << "File name:                         " << _path_library << endl;
//	cout << "Number of Flamelets:               " << _n << endl;
//	cout << "Minimum stoic. scalar diss. rate:  " << _chi_st_min  << " Hz" << endl;
//	cout << "Maximum stoic. scalar diss. rate:  " << _chi_st_max  << " Hz" <<  endl;
	cout << endl;
}

void OpenSMOKE_SootSourceTermLibrary::ErrorMessage(const string error_message)
{
	cout << "OpenSMOKE_SootSourceTermLibrary" << endl;
	cout << "File name:     " << _path_library << endl;
	cout << "Error message: " << error_message << endl;
	exit(-1);
}

void OpenSMOKE_SootSourceTermLibrary::WarningMessage(const string warning_message)
{
	cout << "OpenSMOKE_SootSourceTermLibrary" << endl;
	cout << "File name:       " << _path_library << endl;
	cout << "Warning message: " << warning_message << endl;
	cout << "Press enter to continue..." << endl;
	getchar();
}
