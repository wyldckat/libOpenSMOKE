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
 

#ifndef OpenSMOKE_SootSourceTermLibrary_H
#define OpenSMOKE_SootSourceTermLibrary_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


#include "OpenSMOKE_SootEmpiricalModels.h"
#include "OpenSMOKE_SootSourceTermFlameletLibrary.h"

class OpenSMOKE_SootSourceTermLibrary
{
public:
	OpenSMOKE_SootSourceTermLibrary(void);
	~OpenSMOKE_SootSourceTermLibrary(void);

	void SetLibraryPath(const std::string path_library);
	void SetLogNormalChiDistribution();
	void SetNoFluctuationsExtractionMode();
	void setRobustCalculations();

	void SetNucleationModel(const nucleation_models nucleation_model);
	void SetGrowthModel(const growth_models growth_model);
	void SetAggregationModel(const aggregation_models aggregation_model);
	void SetOxidationModel(const oxidation_models oxidation_model);

	void Setup();



	void GetMeanValues(	const double &csi, const double &csiv2, const double &chi_st, 
						const double &phiN, const double &phiM, const double &rho,
						double &source_phiN, double &source_phiM);

	void GetMeanValues(	const double &csi, const double &csiv2, const double &chi_st, 
						double &source_pdf_nucleation, double &source_pdf_growth,
						double &source_pdf_oxidation, double &source_pdf_aggregation);

	void UpdateProperties(	const double &rho,const double &phiN, const double &phiM,
				double &m0, double &fv, double &dSoot, double &MSoot, double &ASoot);

	FractionalExponents	*gamma_fv;
	FractionalExponents	*gamma_m0;

private:

	std::string _path_library;

	nucleation_models	nucleation_model_;
	growth_models		growth_model_;
	aggregation_models	aggregation_model_;
	oxidation_models	oxidation_model_;
	surface_functions   surface_function_;

	std::string nucleation_path_	;
	std::string growth_path_;
	std::string oxidation_path_;
	std::string aggregation_path_;

	chiPDF_kinds	 		_iChiPDF;
	soot_extraction_modes 	_extraction_mode;
	bool					_iRobustCalculations;

	double mp_;

	OpenSMOKE_SootSourceTermFlameletLibrary *flamelets_library;

private:

	void Summary();
	void ErrorMessage(const std::string error_message);
	void WarningMessage(const std::string warning_message);
};

#endif // OpenSMOKE_SootSourceTermLibrary_H
