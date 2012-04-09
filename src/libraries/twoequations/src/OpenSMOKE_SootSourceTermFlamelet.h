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

#ifndef OpenSMOKE_SootSourceTermFlamelet_H
#define OpenSMOKE_SootSourceTermFlamelet_H

enum soot_extraction_modes {OPENSMOKE_SOOT_SOURCES_EXTRACTION_FLUCTUATIONS, OPENSMOKE_SOOT_SOURCES_EXTRACTION_NO_FLUCTUATIONS};

class OpenSMOKE_SootSourceTermFlamelet
{

public:
	
	OpenSMOKE_SootSourceTermFlamelet();
	~OpenSMOKE_SootSourceTermFlamelet(void);

	void Read(ifstream &fInput, const int index, const std::string file_name);
	void ReadBinary(ifstream &fInput, const int index, const std::string file_name);
	void Summary();
	
	void GetMeanValues(const double csi, const double csiv2, double &source);

	void SetNoFluctuationsExtractionMode();

	inline double chi_st() {return _chi_st;}

private:

	std::string _name;
	int    _index;

	double _chi_st;				// stoichiometric scalar dissipation rate [Hz]
	double _pressure_pa;		// stoichiometric scalar dissipation rate [Hz]

	int _n_csi;					// number of mixture fraction points
	int _n_variance;			// number of mixture fraction variance points;

	std::vector<double>	_csi;						// mixture fraction points (favre)
	std::vector<double>	_variance_normal;			// normal variance points (favre)

	std::vector< std::vector<double> > _source_pdf;			// mixture fraction points (reynolds)
	std::vector< std::vector<double> > _source_nopdf;			// normal variance points (reynolds)

	soot_extraction_modes _extraction_mode;
			
private:

	void ErrorMessage(const std::string error_message);
	void WarningMessage(const std::string warning_message);

};

#endif // #ifndef OpenSMOKE_SootSourceTermFlamelet_H

	
