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

OpenSMOKE_SootSourceTermFlamelet::~OpenSMOKE_SootSourceTermFlamelet(void)
{
}

OpenSMOKE_SootSourceTermFlamelet::OpenSMOKE_SootSourceTermFlamelet()
{
	_name = "NULL";

	_extraction_mode = OPENSMOKE_SOOT_SOURCES_EXTRACTION_FLUCTUATIONS;
}

void OpenSMOKE_SootSourceTermFlamelet::SetNoFluctuationsExtractionMode()
{
	_extraction_mode = OPENSMOKE_SOOT_SOURCES_EXTRACTION_NO_FLUCTUATIONS;
}

void OpenSMOKE_SootSourceTermFlamelet::Read(ifstream &fInput, const int index, const string file_name)
{
	int i, j;
	string tag;

	_name = file_name;
	_index = index;

	// Read pressure in Pa
	fInput >> tag;
	if (tag != "pressure")
		ErrorMessage("Expected: pressure Found: " + tag);
	fInput >> _pressure_pa;

	// Read stoichiometric scalar dissipation rate
	fInput >> tag;
	if (tag != "chi")
		ErrorMessage("Expected: chi Found: " + tag);
	fInput >> _chi_st;

	// Read number of mixture fraction points
	fInput >> tag;	// 1
	if (tag != "nCsi")
		ErrorMessage("Expected: nCsi, Found: " + tag);
	fInput >> _n_csi;
	
	// Read number of mixture variance mixture fraction points
	fInput >> tag;	// 1
	if (tag != "nVariance")
		ErrorMessage("Expected: nVariance, Found: " + tag);
	fInput >> _n_variance;

	// Mixture fraction points
	fInput >> tag;
	if (tag != "mf")
		ErrorMessage("Expected: mf, Found: " + tag);
	_csi.resize(_n_csi+1);
	for(i=1;i<=_n_csi;i++)
		fInput >> _csi[i];

	// Variance of mixture fraction points
	fInput >> tag;
	if (tag != "mfv")
		ErrorMessage("Expected: mfv, Found: " + tag);
	_variance_normal.resize(_n_variance+1);
	for(j=1;j<=_n_variance;j++)
		fInput >> _variance_normal[j];


	// Memory allocation
	_source_pdf.resize(_n_variance+1);
	_source_nopdf.resize(_n_variance+1);
	for(j=1;j<=_n_variance;j++)
	{
		_source_pdf[j].resize(_n_csi+1);
		_source_nopdf[j].resize(_n_csi+1);
	}

	fInput >> tag;	// 3
	if (tag != "sources")
		ErrorMessage("Expected: sources, Found: " + tag);
	for(i=1;i<=_n_csi;i++)
		for(j=1;j<=_n_variance;j++)
		{
			fInput >> _source_pdf[j][i];
			fInput >> _source_nopdf[j][i];
		}
}

void OpenSMOKE_SootSourceTermFlamelet::ReadBinary(ifstream &fInput, const int index, const string file_name)
{
	int i, j;
	const int SIZE = 40;
	char tag[SIZE];

	_name = file_name;
	_index = index;

	// Read total number of species in the flamelet
	fInput.read(tag, SIZE);
	if (string(tag) != "pressure")
		ErrorMessage("Expected: pressure, Found: " + string(tag));
	fInput.read(reinterpret_cast < char * > (&_pressure_pa),sizeof(double));

	fInput.read(tag, SIZE);	
	if (string(tag) != "chi")
		ErrorMessage("Expected: chi, Found: " + string(tag));
	fInput.read(reinterpret_cast < char * > (&_chi_st),sizeof(double));

	fInput.read(tag, SIZE);	
	if (string(tag) != "nCsi")
		ErrorMessage("Expected: nCsi, Found: " + string(tag));
	fInput.read(reinterpret_cast < char * > (&_n_csi),sizeof(int));

	fInput.read(tag, SIZE);	
	if (string(tag) != "nVariance")
		ErrorMessage("Expected: nVariance, Found: " + string(tag));
	fInput.read(reinterpret_cast < char * > (&_n_variance),sizeof(int));

	// Mixture fraction
	fInput.read(tag, SIZE);
	if (string(tag) != "mf")
		ErrorMessage("Expected: mf, Found: " + string(tag));

	fInput.read(reinterpret_cast < char * > (&_n_csi),sizeof(int));
	_csi.resize(_n_csi+1);
	for(i=1;i<=_n_csi;i++)
		fInput.read(reinterpret_cast < char * > (&_csi[i]),sizeof(double));;

	// Variance of mixture fraction
	fInput.read(tag, SIZE);
	if (string(tag) != "mfv")
		ErrorMessage("Expected: mfv, Found: " + string(tag));

	fInput.read(reinterpret_cast < char * > (&_n_variance),sizeof(int));
	_variance_normal.resize(_n_variance+1);
	for(j=1;j<=_n_variance;j++)
		fInput.read(reinterpret_cast < char * > (&_variance_normal[j]),sizeof(double));


	// Memory allocation
	_source_pdf.resize(_n_variance+1);
	_source_nopdf.resize(_n_variance+1);
	for(j=1;j<=_n_variance;j++)
	{
		_source_pdf[j].resize(_n_csi+1);
		_source_nopdf[j].resize(_n_csi+1);
	}

	fInput.read(tag, SIZE);	// 3
	if (string(tag) != "sources")
		ErrorMessage("Expected: sources, Found: " + string(tag));
	for(i=1;i<=_n_csi;i++)
		for(j=1;j<=_n_variance;j++)
		{
			fInput.read(reinterpret_cast < char * > (&_source_pdf[j][i]),sizeof(double));
			fInput.read(reinterpret_cast < char * > (&_source_nopdf[j][i]),sizeof(double));
		}
}

void OpenSMOKE_SootSourceTermFlamelet::GetMeanValues(const double csi, const double csiv2, double &source)
{
	int iCsi = 0;
	int jCsiv2 = 0;
	
	for(int j=2;j<=_n_variance;j++)
		if (csiv2 <= _variance_normal[j])
		{
			jCsiv2 = j-1;
			break;
		}

	for(int i=2;i<=_n_csi;i++)
		if (csi <= _csi[i])
		{
			iCsi = i-1;
			break;
		}
		
	double q   = (_csi[iCsi+1]-_csi[iCsi])*(_variance_normal[jCsiv2+1]-_variance_normal[jCsiv2]);
	double q11 = (_csi[iCsi+1]-csi) * (_variance_normal[jCsiv2+1]-csiv2);
	double q21 = (csi-_csi[iCsi])   * (_variance_normal[jCsiv2+1]-csiv2);
	double q12 = (_csi[iCsi+1]-csi) * (csiv2-_variance_normal[jCsiv2]);
	double q22 = (csi-_csi[iCsi])   * (csiv2-_variance_normal[jCsiv2]);
	
	// Source (pdf)
	if (_extraction_mode == OPENSMOKE_SOOT_SOURCES_EXTRACTION_FLUCTUATIONS)	
	source	= ( 	_source_pdf[jCsiv2][iCsi]*q11   + _source_pdf[jCsiv2][iCsi+1]*q21 +
         				_source_pdf[jCsiv2+1][iCsi]*q12 + _source_pdf[jCsiv2+1][iCsi+1]*q22 ) / q;

	// Source (no pdf)
	if (_extraction_mode == OPENSMOKE_SOOT_SOURCES_EXTRACTION_NO_FLUCTUATIONS)	
		source	= (	_source_nopdf[jCsiv2][iCsi]*q11   + _source_nopdf[jCsiv2][iCsi+1]*q21 +
			     	_source_nopdf[jCsiv2+1][iCsi]*q12 + _source_nopdf[jCsiv2+1][iCsi+1]*q22 ) / q;
}

void OpenSMOKE_SootSourceTermFlamelet::Summary()
{
	cout << endl;
	cout << "OpenSMOKE_PDF_Flamelet" << endl;
	cout << "Index:                             " << _index << endl;
	cout << "File name:                         " << _name << endl;
	cout << "Mixture fraction points:           " << _n_csi << endl;
	cout << "Variance Mixture fraction points:  " << _n_variance << endl;
	cout << "Stoic. scalar dissipation rate:    " << _chi_st << " Hz"  << endl;
	cout << endl;
}

void OpenSMOKE_SootSourceTermFlamelet::ErrorMessage(const string error_message)
{
	cout << "OpenSMOKE_SootSourceTermFlamelet" << endl;
	cout << "File name:     " << _name << endl;
	cout << "Error message: " << error_message << endl;
	exit(-1);
}

void OpenSMOKE_SootSourceTermFlamelet::WarningMessage(const string warning_message)
{
	cout << "OpenSMOKE_SootSourceTermFlamelet" << endl;
	cout << "File name:     " << _name << endl;
	cout << "Warning message: " << warning_message << endl;
	cout << "Press enter to continue..." << endl;
	getchar();
}



	

