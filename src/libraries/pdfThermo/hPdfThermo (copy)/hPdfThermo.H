/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::hPdfThermo

Description
    Enthalpy for a mixture based on compressibility

SourceFiles
    hPdfThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hPdfThermo_H
#define hPdfThermo_H

#include "basicPdfThermo.H"
#include "basicMixture.H"

#include "OpenSMOKE_PDF_NonAdiabaticFlamelet_Library.hpp"
#include "OpenSMOKE_SootEmpiricalModels.h"
#include "OpenSMOKE_SootSourceTermLibrary.h"
#include "OpenSMOKE_QMOM_EmpiricalModels.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hPdfThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class hPdfThermo
:
    public basicPdfThermo,
    public MixtureType
{
private:


	volScalarField h_;

	volScalarField csi_;
	volScalarField csiv2_;
	volScalarField H_;

	volScalarField density_reynolds_;
	volScalarField chi_st_;
	volScalarField phiH_;
	volScalarField as_;		
	volScalarField mu_favre_;	
	volScalarField alpha_favre_;
	PtrList<volScalarField> omega_;
	PtrList<volScalarField> sootPhi_;
	PtrList<volScalarField> sootSource_;	
	PtrList<volScalarField> sootProperties_;

	PtrList<volScalarField> sootMoments_;
	PtrList<volScalarField> sootMomentsSources_;
				
	OpenSMOKE_PDF_NonAdiabaticFlamelet_Library 	flamelets_library;
	OpenSMOKE_SootEmpiricalModels			*soot_2e;
	OpenSMOKE_SootSourceTermLibrary 		*soot_2e_uncorrelated;
	OpenSMOKE_QMOM_EmpiricalModels			*soot_qmom;		

	double HFuel; 
	double HOxidizer;
	label counter;
	label propertyUpdate;
	label counter_mass_fractions;
	label massFractionsUpdate;
	label sootUpdate;
	label counter_soot;
	label nDirac;
	Switch adiabaticMode;
	Switch sootRobust;

	enum {SOOT_2E_NONE, SOOT_2E_MEAN, SOOT_2E_UNCORRELATED, SOOT_2E_CORRELATED} 	iTwoEquationModelClosure;
	enum {SOOT_QMOM_NONE, SOOT_QMOM_MEAN } 						iQMOMModelClosure;

	int soot_index_c2h2;
	int soot_index_h2;
	int soot_index_o2;
	int soot_index_oh;
	int soot_index_co;
	int soot_index_h;



private:

	hPdfThermo(const hPdfThermo<MixtureType>&);

	void calculate();

	virtual void update();
	virtual void updateMassFractions();
	virtual void updateSourceTerms_MeanClosure_2E();
	virtual void updateSourceTerms_UncorrelatedClosure_2E();
	virtual void updateSourceTerms_CorrelatedClosure_2E();
	virtual void updateSourceTerms_MeanClosure_QMOM();

	void ErrorMessage(const string message);


public:

        //- Read thermophysicalProperties dictionary
        virtual bool read();
		
	//- Mixture fraction
	//  Non-const access allowed for transport equations
	virtual volScalarField& csi()
	{
		return csi_;
	}

	//- Mixture fraction
	virtual const volScalarField& csi() const
	{
		return csi_;
	}

	//- Variance of Mixture fraction
	//  Non-const access allowed for transport equations
	virtual volScalarField& csiv2()
	{
		return csiv2_;
	}

	//- Variance of Mixture fraction
	virtual const volScalarField& csiv2() const
	{
		return csiv2_;
	}   

	//- Scalar dissipation rate [1/s]
	//  Non-const access allowed for transport equations
	virtual volScalarField& chi_st()
	{
		return chi_st_;
	}

	//- Scalar dissipation rate [1/s]
	virtual const volScalarField& chi_st() const
	{
		return chi_st_;
	}     

	//- Enthalpy [J/kg]
	//  Non-const access allowed for transport equations
	virtual volScalarField& H()
	{
		return H_;
	}

	//- Enthalpy [J/kg]
	virtual const volScalarField& H() const
	{
		return H_;
	}   

	//- Absorption coefficient [1/m]
	//  Non-const access allowed for transport equations
	virtual volScalarField& as()
	{
		return as_;
	}

	//- Absorption coefficient [1/m]
	virtual const volScalarField& as() const
	{
		return as_;
	}    

	//- phiN
	//  Non-const access allowed for transport equations
	virtual volScalarField& phiN()
	{
		return sootPhi_[0];
	}

	//- phiN
	virtual const volScalarField& phiN() const
	{
		return sootPhi_[0];
	}

	//- phiM
	//  Non-const access allowed for transport equations
	virtual volScalarField& phiM()
	{
		return sootPhi_[1];
	}

	//- phiM
	virtual const volScalarField& phiM() const
	{
		return sootPhi_[1];
	} 

	//- Sn
	//  Non-const access allowed for transport equations
	virtual volScalarField& Sn()
	{
		return sootSource_[0];
	}

	//- Sn
	virtual const volScalarField& Sn() const
	{
		return sootSource_[0];
	}

	//- Sm
	//  Non-const access allowed for transport equations
	virtual volScalarField& Sm()
	{
		return sootSource_[1];
	}

	//- Sm
	virtual const volScalarField& Sm() const
	{
		return sootSource_[1];
	}                                         

	//- Smoments
	//  Non-const access allowed for transport equations
	virtual PtrList<volScalarField>& moments()
	{
		return sootMoments_;
	}

	//- Smoments
	virtual const PtrList<volScalarField>& moments() const
	{
		return sootMoments_;
	}

	//- Smoments
	//  Non-const access allowed for transport equations
	virtual PtrList<volScalarField>& Smoments()
	{ 
		return sootMomentsSources_;
	}

	//- Smoments
	virtual const PtrList<volScalarField>& Smoments() const
	{
		return sootMomentsSources_;
	}



public:

    //- Runtime type information
    TypeName("hPdfThermo");


    // Constructors

        //- Construct from mesh
        hPdfThermo(const fvMesh&);


    //- Destructor
    virtual ~hPdfThermo();


    // Member functions

        //- Return the compostion of the mixture
        virtual basicMixture& composition()
        {
            return *this;
        }

        //- Return the compostion of the mixture
        virtual const basicMixture& composition() const
        {
            return *this;
        }

        //- Update properties
        virtual void correct();


        // Access to thermodynamic state variables

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h()
            {
                return h_;
            }

            //- Enthalpy [J/kg]
            virtual const volScalarField& h() const
            {
                return h_;
            }


        // Fields derived from thermodynamic state variables

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "hPdfThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
