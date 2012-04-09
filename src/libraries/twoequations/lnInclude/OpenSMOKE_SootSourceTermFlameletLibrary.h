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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#ifndef OpenSMOKE_SootSourceTermFlameletLibrary_H
#define OpenSMOKE_SootSourceTermFlameletLibrary_H

#include "OpenSMOKE_SootSourceTermFlamelet.h"
#include "OpenSMOKE_ChiDistribution.hpp"
#include "OpenSMOKE_External_Functions.hpp"

class OpenSMOKE_SootSourceTermFlameletLibrary
{

public:

	OpenSMOKE_SootSourceTermFlameletLibrary();
	~OpenSMOKE_SootSourceTermFlameletLibrary(void);

	void Read();
	void Summary();
	
	void GetMeanValues(const double csi, const double csiv2, const double chi_st, double &source);

	void SetLibraryPath(const std::string path_library);
	void SetLogNormalChiDistribution();
	void SetNoFluctuationsExtractionMode();


private:

	std::string			 _path_library;						// library name
	chiPDF_kinds		 _iChiPDF;
	soot_extraction_modes	 _extraction_mode;

private:

	std::string _name;											// library name
	int    _n;												// number of flamelets

	std::vector<double> _chi_st;									// stoichiometric scalar dissipation rate [Hz]
	std::vector< OpenSMOKE_SootSourceTermFlamelet >	_flame;		// Flamelets
	OpenSMOKE_ChiDistribution				_chi_pdf;	// Scalar dissipation rate probability distribution function

	void GetMeanValuesMultiScalarDissipationRatesDirac(const double csi, const double csiv2, const double chi_st, double &source);
	void GetMeanValuesMultiScalarDissipationRatesLogNormal(const double csi, const double csiv2, const double chi_st, double &source);

private:

	bool   _iMultiScalarDissipationRates;
	
	double _chi_st_min;
	double _chi_st_max;

	std::vector<double>	_source_favre;

private:

	void ErrorMessage(const std::string error_message);
	void WarningMessage(const std::string warning_message);
};

#endif // OpenSMOKE_SootSourceTermFlameletLibrary_H
	

