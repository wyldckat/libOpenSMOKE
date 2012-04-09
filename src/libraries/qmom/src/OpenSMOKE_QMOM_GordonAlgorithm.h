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

#ifndef OpenSMOKE_QMOM_GordonAlgorithm_H
#define OpenSMOKE_QMOM_GordonAlgorithmH

#include "gsl/gsl_eigen.h"

class OpenSMOKE_QMOM_GordonAlgorithm
{
public:
	OpenSMOKE_QMOM_GordonAlgorithm(const int N);
	~OpenSMOKE_QMOM_GordonAlgorithm(void);

	void Apply(double *m, double *r, double *w);

private:

	int N_;
	
	double *alfa_;
	double *a_;
	double *b_;
	double *P_;

	gsl_matrix *J;
	gsl_vector *eigen_values;
	gsl_matrix *eigen_vectors;
	gsl_eigen_symmv_workspace *work;
};

#endif 
