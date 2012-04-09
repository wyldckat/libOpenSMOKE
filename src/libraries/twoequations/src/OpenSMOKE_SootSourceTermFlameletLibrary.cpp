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

using namespace std;


#include "OpenSMOKE_SootSourceTermFlamelet.h"
#include "OpenSMOKE_SootSourceTermFlameletLibrary.h"
#include "OpenSMOKE_ChiDistribution.hpp"

OpenSMOKE_SootSourceTermFlameletLibrary::~OpenSMOKE_SootSourceTermFlameletLibrary(void)
{
}

OpenSMOKE_SootSourceTermFlameletLibrary::OpenSMOKE_SootSourceTermFlameletLibrary()
{
	_name = "[not assigned]";
	_path_library = "constant/PDF-Library";
	_iMultiScalarDissipationRates = true;

	_iChiPDF = CHI_PDF_DIRAC;
	_extraction_mode = OPENSMOKE_SOOT_SOURCES_EXTRACTION_FLUCTUATIONS;
}

void OpenSMOKE_SootSourceTermFlameletLibrary::SetLibraryPath(const string path_library)
{
	_path_library = path_library;
}

void OpenSMOKE_SootSourceTermFlameletLibrary::SetLogNormalChiDistribution()
{
	_iChiPDF = CHI_PDF_LOG_NORMAL;
}

void OpenSMOKE_SootSourceTermFlameletLibrary::SetNoFluctuationsExtractionMode()
{
	_extraction_mode = OPENSMOKE_SOOT_SOURCES_EXTRACTION_NO_FLUCTUATIONS;
}

void OpenSMOKE_SootSourceTermFlameletLibrary::Read()
{
	int i;
	string tag;

	_name = _path_library + "/PDF-0/Source.source";

	ifstream fInput;
	fInput.open(_name.c_str(), ios::in);

	fInput >> tag;
	if (tag != "nFlamelets")
		ErrorMessage("Expected: nFlamelets, Found: " + tag);
	fInput >> _n;

	// Memory allocation
	_chi_st.resize(_n+1);	
	_flame.resize(_n+1);

	for(i=1;i<=_n;i++)
	{
		if (_extraction_mode == OPENSMOKE_SOOT_SOURCES_EXTRACTION_NO_FLUCTUATIONS)
			_flame[i].SetNoFluctuationsExtractionMode();

		_flame[i].Read(fInput, i, _name);
	//	_flame[i].ReadBinary(fInput, i, _name);
		_flame[i].Summary();
	}

	fInput.close();

	// Checking scalar dissipation rates
	for(i=1;i<=_n;i++)
		_chi_st[i] = _flame[i].chi_st();
	for(i=1;i<=_n-1;i++)
		if (_chi_st[i] >= _chi_st[i+1])
			ErrorMessage("Scalar dissipation rates are not in the right order...");
	
	_chi_st_min = _chi_st[1];
	_chi_st_max = _chi_st[_n];

 	if (_n==1)	
		_iMultiScalarDissipationRates = false;
	
	if (_n<3 && _iChiPDF == CHI_PDF_LOG_NORMAL)
		ErrorMessage("Log-normal distribution function for scalar dissipation rate can be used only with at least three flamelets...");
		
	// Prepare Scalar Dissipation rate distribution function
	if (_iChiPDF == CHI_PDF_LOG_NORMAL)
	{
		// Scalar Dissipation Rate Probability distribution function options
		//_chi_pdf.SetSigma(1.31);
		//_chi_pdf.SetLowerLimit(1.e-6);
		//_chi_pdf.SetHigherLimit(1.e3);
		//_chi_pdf.SetNumberOfPoints(40);
		//_chi_pdf.SetFixedPointRatio();
		//_chi_pdf.SetAccurateCalculation();
		//_chi_pdf.SetWeightThreshold(1.e-10);

		_chi_pdf.BuildGrid();
		_chi_pdf.AssignScalarDissipationRates(_chi_st);
		_source_favre.resize(_n+1);
	}
}

void OpenSMOKE_SootSourceTermFlameletLibrary::GetMeanValues(const double csi, const double csiv2, const double chi_st, double &source)
{
	if (_iMultiScalarDissipationRates == true)
	{
		if(_iChiPDF == CHI_PDF_LOG_NORMAL)		
			GetMeanValuesMultiScalarDissipationRatesLogNormal(csi, csiv2, chi_st, source);
		else
			GetMeanValuesMultiScalarDissipationRatesDirac(csi, csiv2, chi_st, source);
	}
	else											
		_flame[1].GetMeanValues(csi, csiv2, source);
}

void OpenSMOKE_SootSourceTermFlameletLibrary::GetMeanValuesMultiScalarDissipationRatesDirac(const double csi, const double csiv2, const double chi_st, double &source)
{
	if (chi_st <= _chi_st_min)	
	{
		_flame[1].GetMeanValues(csi, csiv2, source);
		return;
	}
	else if (chi_st >= _chi_st_max)	
	{
		_flame[_n].GetMeanValues(csi, csiv2, source);
		return;
	}
	else
	{
		int kChi = 0;
		double ratio;
		double source_1, source_2;

		// Finding kChi
		for(int k=2;k<=_n;k++)
			if (chi_st <= _chi_st[k])
			{
				kChi = k-1;
				break;
			}
			
		// Finding flamelets mean values
		_flame[kChi].GetMeanValues(csi, csiv2, source_1);
		_flame[kChi+1].GetMeanValues(csi, csiv2, source_2);
		
		// Interpolation ratio
		ratio = (chi_st - _chi_st[kChi])/(_chi_st[kChi+1] - _chi_st[kChi]);
		
		// Mean values
		source		= source_1 + ratio*(source_2-source_1);
	}
}

void OpenSMOKE_SootSourceTermFlameletLibrary::GetMeanValuesMultiScalarDissipationRatesLogNormal(const double csi, const double csiv2, const double chi_st, double &source)
{
	_chi_pdf.AssignMeanScalarDissipationRate(chi_st);

	for(int j=_chi_pdf.iXinf();j<=_chi_pdf.iXsup();j++)
	{
		_flame[j].GetMeanValues(csi, csiv2, source);
		_source_favre[j] 	= source;
	}

	source	= _chi_pdf.GetMeanValue(_source_favre);
}



void OpenSMOKE_SootSourceTermFlameletLibrary::Summary()
{
	cout << endl;
	cout << "OpenSMOKE_SootSourceTermFlameletLibrary" << endl;
	cout << "File name:                         " << _name << endl;
	cout << "Number of Flamelets:               " << _n << endl;
	cout << "Minimum stoic. scalar diss. rate:  " << _chi_st_min  << " Hz" << endl;
	cout << "Maximum stoic. scalar diss. rate:  " << _chi_st_max  << " Hz" <<  endl;
	cout << endl;
}

void OpenSMOKE_SootSourceTermFlameletLibrary::ErrorMessage(const string error_message)
{
	cout << "OpenSMOKE_SootSourceTermFlameletLibrary" << endl;
	cout << "File name:     " << _name << endl;
	cout << "Error message: " << error_message << endl;
	exit(-1);
}

void OpenSMOKE_SootSourceTermFlameletLibrary::WarningMessage(const string warning_message)
{
	cout << "OpenSMOKE_SootSourceTermFlameletLibrary" << endl;
	cout << "File name:       " << _name << endl;
	cout << "Warning message: " << warning_message << endl;
	cout << "Press enter to continue..." << endl;
	getchar();
}


	
