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

#include "OpenSMOKE_QMOM_GordonAlgorithm.h"

#include <iostream>
#include <math.h>

using namespace std;

OpenSMOKE_QMOM_GordonAlgorithm::OpenSMOKE_QMOM_GordonAlgorithm(const int N)
{
	N_ = N;

	a_		= new double[N_];
	b_		= new double[N_-1];
	alfa_	= new double[2*N_];
	P_		= new double[(2*N_+1)*(2*N_+1)];

	// Matrix P
	for(int i=1;i<=(2*N_+1)*(2*N_+1);i++)
		P_[i-1] = 0.;
	P_[0]	= 1.;
	P_[1]	= 1.;

	// Vector alfa
	alfa_[0] = 0.;

	// Memory allocation
	J = gsl_matrix_alloc(N_,N_);
	eigen_values = gsl_vector_alloc(N_);
	eigen_vectors = gsl_matrix_alloc(N_,N_);
	work = gsl_eigen_symmv_alloc(N_);
}

OpenSMOKE_QMOM_GordonAlgorithm::~OpenSMOKE_QMOM_GordonAlgorithm(void)
{
}

void OpenSMOKE_QMOM_GordonAlgorithm::Apply(double *m, double *r, double *w)
{
	if (N_>2)
	{
		// Second colum
		double c = -1.;
		for(int i=2;i<=2*N_;i++)
		{
			P_[(i-1)*(2*N_+1)+1] = c*m[i-1];
			c*=-1.;
		}

		// Additional columns
		for(int j=3;j<=2*N_+1;j++)
			for(int i=1;i<=2*N_+2-j;i++)
				P_[(i-1)*(2*N_+1)+(j-1)] =	P_[(j-2)]*P_[(i)*(2*N_+1)+(j-3)] - 
											P_[(j-3)]*P_[(i)*(2*N_+1)+(j-2)];

		for(int j=2;j<=2*N_;j++)
			alfa_[j-1] = P_[j]/P_[j-1]/P_[j-2];

		for(int j=1;j<=N_;j++)
			a_[j-1] = alfa_[2*j-1]+alfa_[2*j-2];

		for(int j=1;j<=N_-1;j++)
			b_[j-1] = sqrt(alfa_[2*j]*alfa_[2*j-1]);


		// Eigenvalue evaluations
		gsl_matrix_set_zero(J);
		for(int j=0;j<N_;j++)
			gsl_matrix_set(J,j,j,a_[j]);
		for(int j=0;j<N_-1;j++)
			gsl_matrix_set(J,j,j+1,b_[j]);
		for(int j=0;j<N_-1;j++)
			gsl_matrix_set(J,j+1,j,b_[j]);

		gsl_eigen_symmv(J, eigen_values, eigen_vectors, work);
		
		for(int j=0;j<N_;j++)
		{
			r[j] = gsl_vector_get(eigen_values, j);
			w[j] = gsl_matrix_get(eigen_vectors, 0,j)*gsl_matrix_get(eigen_vectors, 0,j);
		}
	}

	else
	{
		alfa_[1] = m[1];
		alfa_[2] = (m[2]-m[1]*m[1])/m[1];
		alfa_[3] = (m[1]*m[3]-m[2]*m[2])/(m[1]*m[2]-m[1]*m[1]*m[1]);

		a_[0] = m[1];
		a_[1] = alfa_[2]+alfa_[3];

		b_[0] = sqrt( fabs(m[2]-m[1]*m[1]) );

		double sqrt_coefficient = sqrt( (a_[0]-a_[1])*(a_[0]-a_[1]) + 4.*b_[0]*b_[0] );

		r[0] = 0.50*( (a_[0]+a_[1]) - sqrt_coefficient);
		r[1] = 0.50*( (a_[0]+a_[1]) + sqrt_coefficient);
		
		{
			double nx = r[0]-a_[0];
			double ny = b_[0];
			double sum = sqrt(nx*nx+ny*ny);
			nx /= sum;
			w[1] = nx*nx;
		}
		{
			double nx = r[1]-a_[0];
			double ny = b_[0];
			double sum = sqrt(nx*nx+ny*ny);
			nx /= sum;
			w[0] = nx*nx;
		}

		// Check (to remove)
		if ( fabs((w[0]*r[0]+w[1]*r[1]) - m[1]) / fabs(m[1])> 1e-6)
		{
			cout << "Error in evaluating w and csi (2x2)" << endl;
			exit(-1);
		}
	}
}
