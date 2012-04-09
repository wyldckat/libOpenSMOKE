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

\*---------------------------------------------------------------------------*/

#include "hPdfThermo.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::calculate()
{
        const scalarField& pCells		= p_.internalField();
	const scalarField& HCells 	= H_.internalField();
	const scalarField& csiCells 	= csi_.internalField();
	const scalarField& RhoReynolds 	= density_reynolds_.internalField();
	const scalarField& muFavre 	= mu_favre_.internalField();
	const scalarField& alphaFavre 	= alpha_favre_.internalField();
	
    scalarField& psiCells 		= psi_.internalField();
    scalarField& muCells 		= mu_.internalField();
    scalarField& alphaCells 	= alpha_.internalField();
	scalarField& phiHCells 		= phiH_.internalField();
	
    forAll(csiCells, celli)
    {
       // const typename MixtureType::thermoType& mixture_ = this->cellMixture(celli);

        psiCells[celli] 	= RhoReynolds[celli]/pCells[celli];
        
        muCells[celli] 		= muFavre[celli];
        alphaCells[celli] 	= alphaFavre[celli];
        
        phiHCells[celli]    = HCells[celli] - (HOxidizer+csiCells[celli]*(HFuel-HOxidizer));
    }
    
	// Boundaries
    forAll(T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp			= p_.boundaryField()[patchi];
        const fvPatchScalarField& pH 			= H_.boundaryField()[patchi];
        const fvPatchScalarField& pcsi 			= csi_.boundaryField()[patchi];
 		const fvPatchScalarField& pRhoReynolds 	= density_reynolds_.boundaryField()[patchi];
		const fvPatchScalarField& pmuFavre 		= mu_favre_.boundaryField()[patchi];
		const fvPatchScalarField& palphaFavre 	= alpha_favre_.boundaryField()[patchi];
        
        fvPatchScalarField& pT 		= T_.boundaryField()[patchi];
        fvPatchScalarField& ppsi 	= psi_.boundaryField()[patchi];   
        fvPatchScalarField& pphiH 	= phiH_.boundaryField()[patchi];
        fvPatchScalarField& pmu 	= mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha 	= alpha_.boundaryField()[patchi];


        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
               // const typename MixtureType::thermoType& mixture_ = this->patchFaceMixture(patchi, facei);

                //ph[facei] 		= mixture_.H(pT[facei]); // TODO
                
                ppsi[facei] 	= pRhoReynolds[facei]/pp[facei]; 
                
                pmu[facei] 		= pmuFavre[facei];
        		palpha[facei] 	= palphaFavre[facei];
                
                pphiH[facei]    = pH[facei] - (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
            }
        }
        else
        {
            forAll(pT, facei)
            {
               // const typename MixtureType::thermoType& mixture_ = this->patchFaceMixture(patchi, facei);

                ppsi[facei] 	= pRhoReynolds[facei]/pp[facei]; 
       
                pmu[facei] 		= pmuFavre[facei];
        		palpha[facei] 	= palphaFavre[facei];
                
                pphiH[facei]    = pH[facei] - (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hPdfThermo<MixtureType>::hPdfThermo(const fvMesh& mesh)
:
    basicPdfThermo(mesh),
    MixtureType(*this, mesh),

    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        hBoundaryTypes()
    ), 
    
    csi_
    (
        IOobject
        (
            "csi",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    csiv2_
    (
        IOobject
        (
            "csiv2",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    H_
    (
        IOobject
        (
            "H",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),   

    density_reynolds_
    (
        IOobject
        (
	   "density_reynolds",
	   mesh.time().timeName(),
	   mesh,
	   IOobject::NO_READ,
	   IOobject::NO_WRITE
        ),
        mesh,
		dimensionedScalar("density",dimensionSet(1,-3,0,0,0,0,0) , 0.0)
    ),      
    
    chi_st_
    (
        IOobject
        (
            "chi_st",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("chi_st",dimensionSet(0,0,-1,0,0,0,0) , 0.0)
    ),  
    
    phiH_
    (
        IOobject
        (
            "phiH",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phiH",dimensionSet(0,2,-2,0,0,0,0) , 0.0)
    ),
    
    as_
    (
        IOobject
        (
            "as",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("as",dimensionSet(0,-1,0,0,0,0,0) , 0.0)
   ),
   
    mu_favre_
    (
        IOobject
        (
            "mu_lam",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mu_lam",dimensionSet(1,-1,-1,0,0,0,0) , 0.0)
   ),
   
    alpha_favre_
    (
        IOobject
        (
            "alpha_lam",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha_lam",dimensionSet(0,2,-1,0,0,0,0) , 0.0)
   ),
   
   adiabaticMode(false),
   sootRobust(false)
{
     IOdictionary flameletsProperties_
    (
        IOobject
        (
            "flameletsProperties",
            csi_.time().constant(),
            csi_.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    	string libraryPath = flameletsProperties_.lookup("libraryPath");
    	string chiPDF = flameletsProperties_.lookup("pdf");
    	string list_of_species = flameletsProperties_.lookup("species");
    

	Switch adiabaticMode_(flameletsProperties_.lookup("adiabaticMode"));
    	adiabaticMode 		= adiabaticMode_;
   	propertyUpdate 		= readLabel(flameletsProperties_.lookup("propertyUpdate"));
   	massFractionsUpdate 	= readLabel(flameletsProperties_.lookup("massFractionsUpdate"));

	counter = propertyUpdate;
	counter_mass_fractions = massFractionsUpdate;
	counter_soot = -1;
	


	// Options
	flamelets_library.SetLibraryPath(libraryPath);
	flamelets_library.SetSpeciesToExtract(list_of_species);	
	
	// Initialize mass fraction fields
	omega_.setSize(flamelets_library.number_of_species());
	for (int j=0;j<flamelets_library.number_of_species();j++)
	{	
		std::string name_of_species = "omega_" + flamelets_library.species()[j+1];
		omega_.set(	j, new volScalarField( 	IOobject(name_of_species, 
											mesh.time().timeName(), 
											mesh, 
											IOobject::NO_READ, 
											IOobject::AUTO_WRITE),
											mesh,  
											dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0) 
										 ) 
					);
	}
	
	// Set adiabatic mode
	if (adiabaticMode == true)
	{
		Info << "Adiabatic mode..." << endl;
		flamelets_library.SetAdiabaticMode();
	}
		
	// Scalar dissipation rate distribution
	if (chiPDF == "logNormal")
	{
		Info << "Log-Normal Distribution for scalar dissipation rate..." << endl;
		
		scalar chi_sigma(readScalar(flameletsProperties_.lookup("sigma")));
    	scalar chi_points(readScalar(flameletsProperties_.lookup("points")));
    	Info << chi_sigma << " " << chi_points << endl;
    
		flamelets_library.UnsetExcludeColdFlamelets();
		flamelets_library.SetLogNormalChiDistribution();
	}
	else if (chiPDF == "dirac")
	{
		Info << "Delta-Dirac Distribution for scalar dissipation rate..." << endl;
	}
	
	flamelets_library.Read();
	flamelets_library.Summary();
	

	// Update basic fields
	Info << "Initialize basic fields... " << endl;
	HOxidizer = flamelets_library.enthalpy_f_oxidizer();
    	HFuel = flamelets_library.enthalpy_f_fuel();
	update();

	// Initialize enthalpy field
	Info << "Initialize enthalpy field (field)... " << endl;
    	scalarField& hCells = h_.internalField();
    	const scalarField& HCells = H_.internalField();
    	forAll(hCells, celli)
    	{
		hCells[celli] = HCells[celli];
    	}
    	forAll(h_.boundaryField(), patchi)
    	{
		h_.boundaryField()[patchi] == h(T_.boundaryField()[patchi], patchi);
    	}
    	hBoundaryCorrection(h_);

	Info << "Calculate (first time)... " << endl;
	calculate();
    	psi_.oldTime();   // Switch on saving old time
    
    	Info << " Fuel  enthalpy: " << HFuel 	 << " J/kg" << endl;
    	Info << " Oxid. enthalpy: " << HOxidizer << " J/kg" << endl;

	// Soot Models

	iTwoEquationModelClosure = SOOT_2E_NONE;
	iQMOMModelClosure	 = SOOT_QMOM_NONE;

	string sootMode = flameletsProperties_.lookup("sootMode");

	if (sootMode == "2E")
	{
		string soot_closure = flameletsProperties_.lookup("sootClosure");
	
	     	if (soot_closure == "Mean")			iTwoEquationModelClosure = SOOT_2E_MEAN;
		else if (soot_closure == "Uncorrelated")	iTwoEquationModelClosure = SOOT_2E_UNCORRELATED;
		else if (soot_closure == "Correlated")		iTwoEquationModelClosure = SOOT_2E_CORRELATED;
		else ErrorMessage("Wrong soot model closure...");
	}

	else if (sootMode == "QMOM")
	{
		string soot_closure = flameletsProperties_.lookup("sootClosure");
	
	     	if (soot_closure == "Mean")			iQMOMModelClosure = SOOT_QMOM_MEAN;
		else ErrorMessage("Wrong soot model closure...");
	}

	if (iTwoEquationModelClosure != SOOT_2E_NONE || iQMOMModelClosure != SOOT_QMOM_NONE)
	{

		std::vector<string> soot_properties_names(5);
		soot_properties_names[0] = "n0";
		soot_properties_names[1] = "fv";
		soot_properties_names[2] = "dSoot";
		soot_properties_names[3] = "MSoot";
		soot_properties_names[4] = "ASoot";
		sootProperties_.setSize(5);
		for (int j=0;j<5;j++)
			sootProperties_.set(	j, new volScalarField( 	IOobject(soot_properties_names[j], 
										mesh.time().timeName(), 
										mesh, 
										IOobject::NO_READ, 
										IOobject::AUTO_WRITE),
										mesh,  
										dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0) 
									) 
						);

		counter_soot = sootUpdate;
   		sootUpdate = readLabel(flameletsProperties_.lookup("sootUpdate"));

		Switch sootRobust_(flameletsProperties_.lookup("sootRobust"));
    		sootRobust = sootRobust_;

    		string nucleation_model = flameletsProperties_.lookup("nucleationModel");
    		string growth_model = flameletsProperties_.lookup("growthModel");
    		string aggregation_model = flameletsProperties_.lookup("aggregationModel");
    		string oxidation_model = flameletsProperties_.lookup("oxidationModel");

		nucleation_models	nucleation_model_ 		= NUCLEATION_NONE;
		growth_models		growth_model_			= GROWTH_NONE;
		aggregation_models	aggregation_model_		= AGGREGATION_NONE;
		oxidation_models	oxidation_model_		= OXIDATION_NONE;
					
		     if (nucleation_model == "None")		nucleation_model_ = NUCLEATION_NONE;
		else if (nucleation_model == "Liu_2001")	nucleation_model_ = NUCLEATION_LIU_2001;
		else if (nucleation_model == "Liu_2002")	nucleation_model_ = NUCLEATION_LIU_2002;
		else if (nucleation_model == "Moss_1999")	nucleation_model_ = NUCLEATION_MOSS_1999;
		else if (nucleation_model == "Wen_2003")	nucleation_model_ = NUCLEATION_WEN_2003;
    		else if (nucleation_model == "Lindstedt_1994")	nucleation_model_ = NUCLEATION_LINDSTEDT_1994;
		else if (nucleation_model == "Leung_1991")	nucleation_model_ = NUCLEATION_LEUNG_1991;
		else ErrorMessage("Wrong soot nucleation model...");

		     if (growth_model == "None")		growth_model_ = GROWTH_NONE;
		else if (growth_model == "Liu_2001")		growth_model_ = GROWTH_LIU_2001;
		else if (growth_model == "Liu_2002")		growth_model_ = GROWTH_LIU_2002;
		else if (growth_model == "Moss_1999")		growth_model_ = GROWTH_MOSS_1999;
		else if (growth_model == "Wen_2003")		growth_model_ = GROWTH_WEN_2003;
    		else if (growth_model == "Lindstedt_1994")	growth_model_ = GROWTH_LINDSTEDT_1994;
		else if (growth_model == "Leung_1991")		growth_model_ = GROWTH_LEUNG_1991;
		else ErrorMessage("Wrong soot growth model...");

		     if (aggregation_model == "None")		aggregation_model_ = AGGREGATION_NONE;
    		else if (aggregation_model == "Smoluchowski")	aggregation_model_ = AGGREGATION_SMOLUCHOWSKI;
		else if (aggregation_model == "Moss")		aggregation_model_ = AGGREGATION_MOSS;
		else ErrorMessage("Wrong soot aggregation model...");

		     if (oxidation_model == "None")		oxidation_model_ = OXIDATION_NONE;
    		else if (oxidation_model == "Lee")		oxidation_model_ = OXIDATION_LEE;
		else if (oxidation_model == "Neoh")		oxidation_model_ = OXIDATION_NEOH;
		else if (oxidation_model == "NSC")		oxidation_model_ = OXIDATION_NSC;
		else ErrorMessage("Wrong soot oxidation model...");

		if (iQMOMModelClosure == SOOT_QMOM_NONE)
		{
			sootPhi_.setSize(2);
			sootPhi_.set(	0, new volScalarField( 	IOobject("phiN", 
									 mesh.time().timeName(), 
									 mesh, 
									 IOobject::MUST_READ, 
									 IOobject::AUTO_WRITE
									),
									mesh
								) 
						);

			sootPhi_.set(	1, new volScalarField( 	IOobject("phiM", 
									 mesh.time().timeName(), 
									 mesh, 
									 IOobject::MUST_READ, 
									 IOobject::AUTO_WRITE
									),
									mesh
								) 
						);

			sootSource_.setSize(2);
			sootSource_.set(	0, new volScalarField( 	IOobject("Sn", 
												mesh.time().timeName(), 
												mesh, 
												IOobject::NO_READ, 
												IOobject::AUTO_WRITE),
												mesh,  
												dimensionedScalar("zero", dimensionSet(0,-3,-1,0,1,0,0), 0.0) 
											 ) 
						);

		
			sootSource_.set(	1, new volScalarField( 	IOobject("Sm", 
												mesh.time().timeName(), 
												mesh, 
												IOobject::NO_READ, 
												IOobject::AUTO_WRITE),
												mesh,  
												dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0,0,0), 0.0) 
											 ) 
						);

			if (iTwoEquationModelClosure == SOOT_2E_MEAN)
			{
				soot_index_c2h2 = flamelets_library.index_of_species("C2H2")-1;
				soot_index_h2 	= flamelets_library.index_of_species("H2")-1;
				soot_index_o2 	= flamelets_library.index_of_species("O2")-1;
				soot_index_oh 	= flamelets_library.index_of_species("OH")-1;
				soot_index_co 	= flamelets_library.index_of_species("CO")-1;
				soot_index_h 	= flamelets_library.index_of_species("H")-1;

				soot_2e = new OpenSMOKE_SootEmpiricalModels(nucleation_model_, growth_model_, oxidation_model_, aggregation_model_);
				
				if (sootRobust == true)	soot_2e->setRobustCalculations();

				updateMassFractions();
				updateSourceTerms_MeanClosure_2E();		
			}
			else if (iTwoEquationModelClosure == SOOT_2E_UNCORRELATED)
			{

				string sootLibraryPath = flameletsProperties_.lookup("sootLibraryPath");

				soot_2e_uncorrelated = new OpenSMOKE_SootSourceTermLibrary();
				soot_2e_uncorrelated->SetNucleationModel(nucleation_model_);
				soot_2e_uncorrelated->SetGrowthModel(growth_model_);
				soot_2e_uncorrelated->SetAggregationModel(aggregation_model_);
				soot_2e_uncorrelated->SetOxidationModel(oxidation_model_);

				soot_2e_uncorrelated->SetLibraryPath(sootLibraryPath);
				if (sootRobust == true)	soot_2e_uncorrelated->setRobustCalculations();

				//soot_2e_uncorrelated->SetLogNormalChiDistribution();
				//soot_2e_uncorrelated->SetNoFluctuationsExtractionMode();
	
				soot_2e_uncorrelated->Setup();

				updateSourceTerms_UncorrelatedClosure_2E();		
			}
		}

		if (iTwoEquationModelClosure == SOOT_2E_NONE)
		{
			soot_index_c2h2 = flamelets_library.index_of_species("C2H2")-1;
			soot_index_h2 	= flamelets_library.index_of_species("H2")-1;
			soot_index_o2 	= flamelets_library.index_of_species("O2")-1;
			soot_index_oh 	= flamelets_library.index_of_species("OH")-1;
			soot_index_co 	= flamelets_library.index_of_species("CO")-1;
			soot_index_h 	= flamelets_library.index_of_species("H")-1;

			nDirac = readLabel(flameletsProperties_.lookup("sootDirac"));
			soot_qmom = new OpenSMOKE_QMOM_EmpiricalModels(nDirac, nucleation_model_, growth_model_, oxidation_model_, aggregation_model_);

			sootMoments_.setSize(2*nDirac);
			sootMomentsSources_.setSize(2*nDirac);

			for(int j=0;j<2*nDirac;j++)
			{
				stringstream index; index << j;
				sootMoments_.set(	j, new volScalarField( 	IOobject("m"+index.str(), 
											mesh.time().timeName(), 
											mesh, 
											IOobject::MUST_READ, 
											IOobject::AUTO_WRITE
										),
										mesh
									) 
							);
				
				label unit =-3+j;
				sootMomentsSources_.set(	j, new volScalarField( 	IOobject("Sm"+index.str(), 
												mesh.time().timeName(), 
												mesh, 
												IOobject::NO_READ, 
												IOobject::NO_WRITE),
												mesh,  
												dimensionedScalar("zero", dimensionSet(0,unit,-1,0,0,0,0), 0.0) 
											 ) 
						);
			}

			Info << "Update mass..." << endl;
			updateMassFractions();
			Info << "Update sources..." << endl;
			updateSourceTerms_MeanClosure_QMOM();
			Info << "Done" << endl;	
		}
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hPdfThermo<MixtureType>::~hPdfThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering pdfThermo<MixtureType>::correct()" << endl;
    }

    // force the saving of the old-time values
    psi_.oldTime();

	if (counter == propertyUpdate)
	{
		Info << "Updating look-up table extractions..." << endl;
		update();
		counter = 0;
	}
	
	if (counter_mass_fractions == massFractionsUpdate)
	{
		Info << "Updating mass fraction extractions..." << endl;
		updateMassFractions();
		counter_mass_fractions = 0;
	}

	if (counter_soot == sootUpdate)
	{
		Info << "Updating soot source terms..." << endl;

		     if (iTwoEquationModelClosure == SOOT_2E_MEAN) 		updateSourceTerms_MeanClosure_2E();
		else if (iTwoEquationModelClosure == SOOT_2E_UNCORRELATED)	updateSourceTerms_UncorrelatedClosure_2E();
		else if (iTwoEquationModelClosure == SOOT_2E_CORRELATED)	updateSourceTerms_CorrelatedClosure_2E();

		     if (iQMOMModelClosure == SOOT_QMOM_MEAN) 			updateSourceTerms_MeanClosure_QMOM();
				

		counter_soot = 0;
	}	
	
    calculate();

	counter++;
	counter_mass_fractions++;
	counter_soot++;

    if (debug)
    {
        Info<< "exiting pdfThermo<MixtureType>::correct()" << endl;
    }
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hPdfThermo<MixtureType>::h
(
    const scalarField& T,
    const labelList& cells
) const
{
	Info << "pdfThermo h(const scalarField& T, const labelList& cells)...";

    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();
    
    const scalarField& HCells = H_.internalField();
    forAll(T, celli)
    {
        //h[celli] = this->cellMixture(cells[celli]).H(T[celli]);
        h[celli] = HCells[celli];
    }
	Info << " end.." << endl;
    return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hPdfThermo<MixtureType>::h
(
    const scalarField& T,
    const label patchi
) const
{
	Info << "pdfThermo h(const scalarField& T, const label patchi)...";

    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

	 const fvPatchScalarField& pH = H_.boundaryField()[patchi];
    forAll(T, facei)
    {
        //h[facei] = this->patchFaceMixture(patchi, facei).H(T[facei]);
        h[facei] = pH[facei];
    }

	Info << " end.." << endl;
    return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hPdfThermo<MixtureType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
	Info << "pdfThermo Cp(const scalarField& T, const label patchi)...";

    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp();

    forAll(T, facei)
    {
        cp[facei] = this->patchFaceMixture(patchi, facei).Cp(T[facei]);
    }

	Info << " end.." << endl;
    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::hPdfThermo<MixtureType>::Cp() const
{
	Info << "pdfThermo Cp()...";

    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cp = tCp();

    forAll(T_, celli)
    {
        cp[celli] = this->cellMixture(celli).Cp(T_[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        cp.boundaryField()[patchi] = Cp(T_.boundaryField()[patchi], patchi);
    }

	Info << " end.." << endl;
    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hPdfThermo<MixtureType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
	Info << "pdfThermo Cv(const scalarField& T, const label patchi)...";

    tmp<scalarField> tCv(new scalarField(T.size()));
    scalarField& cv = tCv();

    forAll(T, facei)
    {
        cv[facei] = this->patchFaceMixture(patchi, facei).Cv(T[facei]);
    }

	Info << " end.." << endl;
    return tCv;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::hPdfThermo<MixtureType>::Cv() const
{
	Info << "pdfThermo Cv()...";

    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0),
            T_.boundaryField().types()
        )
    );

    volScalarField& cv = tCv();

    forAll(T_, celli)
    {
        cv[celli] = this->cellMixture(celli).Cv(T_[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];
        fvPatchScalarField& pCv = cv.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pCv[facei] = this->patchFaceMixture(patchi, facei).Cv(pT[facei]);
        }
    }

	Info << " end.." << endl;
    return tCv;
}


template<class MixtureType>
bool Foam::hPdfThermo<MixtureType>::read()
{
    if (basicPdfThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::update()
{
	std::vector<double> extracted(7);
	double csiv2_normalized = 0.;
	double defect = 0.; 
	
	
	const scalarField& csi 		= csi_.internalField();
	const scalarField& csiv2 	= csiv2_.internalField();
	const scalarField& chi_st 	= chi_st_.internalField();
	const scalarField& HCells 	= H_.internalField();
	
	scalarField& TCells 		= T_.internalField();
	scalarField& RhoCells 		= density_reynolds_.internalField();
	scalarField& asCells 		= as_.internalField();
	scalarField& muCells 		= mu_favre_.internalField();
	scalarField& alphaCells 	= alpha_favre_.internalField();

	double small_eps = 1.e-6;
	double small_chi_st = 1.e-8;
 	
    forAll(csi, celli)
    {
    	double max_chi = max(small_chi_st,chi_st[celli]);
    	
    	if (adiabaticMode == false)	
    		defect = HCells[celli] - (HOxidizer+csi[celli]*(HFuel-HOxidizer));
    
		if (csi[celli]<=small_eps)			// Pure oxidizer
		{
			flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
		}			
		
		else if (csi[celli]>=(1.-small_eps))	// Pure fuel	
		{
			flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
		}			
		
		else									// Mixture
		{
			csiv2_normalized = csiv2[celli] / (csi[celli]*(1.-csi[celli]));	// Normalized mixture fraction variance
			
			if (csiv2_normalized >= 0.98)
				flamelets_library.GetMeanValues(csi[celli], 0.98, max_chi, defect, extracted);
			else if (csiv2_normalized < 0.)
				flamelets_library.GetMeanValues(csi[celli], 0.00, max_chi, defect, extracted);
			else
				flamelets_library.GetMeanValues(csi[celli], csiv2_normalized, max_chi, defect, extracted);	
		}			
	
		TCells[celli] 		= extracted[1];
		RhoCells[celli]		= extracted[2]; 
		asCells[celli]		= extracted[3]; 
		muCells[celli]		= extracted[4]; 
		alphaCells[celli]	= extracted[5]; 
	}

	std::vector<int> patch_type;
    	forAll(T_.boundaryField(), patchi)
    	{
		if (isA<fixedValueFvPatchScalarField>(T_.boundaryField()[patchi]))
		{
			Info << "Patch Fixed Temperature " << patchi << endl;
			patch_type.push_back(1);
		}
		else
			patch_type.push_back(0);

    	}


    if (adiabaticMode == true)
    {
	 forAll(csi_.boundaryField(), patchi)
	 {
		const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
		const fvPatchScalarField& pcsiv2 	= csiv2_.boundaryField()[patchi];
		const fvPatchScalarField& pchi_st 	= chi_st_.boundaryField()[patchi];
	
		fvPatchScalarField& pt		= T_.boundaryField()[patchi];
		fvPatchScalarField& prho 	= density_reynolds_.boundaryField()[patchi];
		fvPatchScalarField& pas 	= as_.boundaryField()[patchi];
		fvPatchScalarField& pmu 	= mu_favre_.boundaryField()[patchi];
		fvPatchScalarField& palpha 	= alpha_favre_.boundaryField()[patchi];

		forAll(pcsi, facei)
		{
	
			double max_chi = max(small_chi_st, pchi_st[facei]);
	
			if (pcsi[facei]<=0.)		// Pure oxidizer
			{
				flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
			}			
	
			else if (pcsi[facei]>=1.)	// Pure fuel	
			{
				flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
			}			
	
			else						// Mixture
			{
				csiv2_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));	// Normalized mixture fraction variance
	
				if (csiv2_normalized >= 0.98)
					flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
				else if (csiv2_normalized < 0.)
					flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
				else
					flamelets_library.GetMeanValues(pcsi[facei], csiv2_normalized, max_chi,defect, extracted);
			}			

       			pt[facei] 	= extracted[1];
			prho[facei]	= extracted[2];
			pas[facei]	= extracted[3];
			pmu[facei]	= extracted[4];
			palpha[facei]	= extracted[5];
	    	}
	    }
    }
    else
    {
	    forAll(csi_.boundaryField(), patchi)
	    {
		if (patch_type[patchi] == 0)
		{
			const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
			const fvPatchScalarField& pcsiv2 	= csiv2_.boundaryField()[patchi];
			const fvPatchScalarField& pchi_st 	= chi_st_.boundaryField()[patchi];
			const fvPatchScalarField& ph		= H_.boundaryField()[patchi];
		
				fvPatchScalarField& pt		= T_.boundaryField()[patchi];
			fvPatchScalarField& prho 	= density_reynolds_.boundaryField()[patchi];
			fvPatchScalarField& pas 	= as_.boundaryField()[patchi];
			fvPatchScalarField& pmu 	= mu_favre_.boundaryField()[patchi];
			fvPatchScalarField& palpha 	= alpha_favre_.boundaryField()[patchi];

			forAll(pcsi, facei)
			{
		
				double max_chi = max(small_chi_st, pchi_st[facei]);
			
	    			defect = ph[facei] - (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
	    		
				if (pcsi[facei]<=0.)		// Pure oxidizer
				{
					flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
				}			
		
				else if (pcsi[facei]>=1.)	// Pure fuel	
				{
					flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
				}			
		
				else						// Mixture
				{
					csiv2_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));	// Normalized mixture fraction variance
		
					if (csiv2_normalized >= 0.98)
						flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
					else if (csiv2_normalized < 0.)
						flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
					else
						flamelets_library.GetMeanValues(pcsi[facei], csiv2_normalized, max_chi,defect, extracted);
				}			
	
	       			pt[facei] 	= extracted[1];
				prho[facei]	= extracted[2];
				pas[facei]	= extracted[3];
				pmu[facei]	= extracted[4];
				palpha[facei]	= extracted[5];
		    	}
		}
	
		else if (patch_type[patchi] == 1)
		{
			const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
			const fvPatchScalarField& pcsiv2 	= csiv2_.boundaryField()[patchi];
			const fvPatchScalarField& pchi_st 	= chi_st_.boundaryField()[patchi];
			      fvPatchScalarField& ph		= H_.boundaryField()[patchi];
		
			fvPatchScalarField& pt		= T_.boundaryField()[patchi];
			fvPatchScalarField& prho 	= density_reynolds_.boundaryField()[patchi];
			fvPatchScalarField& pas 	= as_.boundaryField()[patchi];
			fvPatchScalarField& pmu 	= mu_favre_.boundaryField()[patchi];
			fvPatchScalarField& palpha 	= alpha_favre_.boundaryField()[patchi];

			forAll(pcsi, facei)
			{
		
				double max_chi = max(small_chi_st, pchi_st[facei]);
		
				if (pcsi[facei]<=0.)		// Pure oxidizer
				{
					defect = flamelets_library.GetEnthalpyDefectFromTemperature(0., 0., max_chi, pt[facei]);	
					ph[facei] = defect + HOxidizer;
					flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
				}			
		
				else if (pcsi[facei]>=1.)	// Pure fuel	
				{
					defect = flamelets_library.GetEnthalpyDefectFromTemperature(1., 0., max_chi, pt[facei]);
					ph[facei] = defect + HFuel;
					flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
				}			
		
				else				// Mixture
				{
					csiv2_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));	// Normalized mixture fraction variance
		
					if (csiv2_normalized >= 0.98)
					{
						defect = flamelets_library.GetEnthalpyDefectFromTemperature(pcsi[facei], 0.98, max_chi, pt[facei]);	
						ph[facei] = defect + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
						flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
					}					
					else if (csiv2_normalized < 0.)
					{
						defect = flamelets_library.GetEnthalpyDefectFromTemperature(pcsi[facei], 0.00, max_chi, pt[facei]);	
						ph[facei] = defect + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
						flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
					}					
					else
					{
						defect = flamelets_library.GetEnthalpyDefectFromTemperature(pcsi[facei], csiv2_normalized, max_chi, pt[facei]);	
						ph[facei] = defect + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
						flamelets_library.GetMeanValues(pcsi[facei], csiv2_normalized, max_chi,defect, extracted);
					}
				}			
	
	       			//pt[facei] 	= extracted[1];
				prho[facei]	= extracted[2];
				pas[facei]	= extracted[3];
				pmu[facei]	= extracted[4];
				palpha[facei]	= extracted[5];
		    	}
		} // patch_type 1
         } // patch cycle
    } // non adiabatic mode
}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::updateMassFractions()
{
	std::vector<double> extracted(flamelets_library.number_of_species()+1);
	double csiv2_normalized = 0.;
	double defect = 0.; 
	
	const scalarField& csi 		= csi_.internalField();
	const scalarField& csiv2 	= csiv2_.internalField();
	const scalarField& chi_st 	= chi_st_.internalField();
	const scalarField& HCells 	= H_.internalField();

	double small_eps = 1.e-6;
	double small_chi_st = 1.e-8;
 	
 	forAll(csi, celli)
    {
    	double max_chi = max(small_chi_st,chi_st[celli]);
    	
    	if (adiabaticMode == false)	
    		defect = HCells[celli] - (HOxidizer+csi[celli]*(HFuel-HOxidizer));
    
		if (csi[celli]<=small_eps)			// Pure oxidizer
		{
			flamelets_library.ExtractMeanValues(0., 0., max_chi, defect, extracted);
		}			
		
		else if (csi[celli]>=(1.-small_eps))	// Pure fuel	
		{
			flamelets_library.ExtractMeanValues(1., 0., max_chi, defect, extracted);
		}			
		
		else									// Mixture
		{
			csiv2_normalized = csiv2[celli] / (csi[celli]*(1.-csi[celli]));	// Normalized mixture fraction variance
			
			if (csiv2_normalized >= 0.98)
				flamelets_library.ExtractMeanValues(csi[celli], 0.98, max_chi, defect, extracted);
			else if (csiv2_normalized < 0.)
				flamelets_library.ExtractMeanValues(csi[celli], 0.00, max_chi, defect, extracted);
			else
				flamelets_library.ExtractMeanValues(csi[celli], csiv2_normalized, max_chi, defect, extracted);	
		}			
	
		for(int j=0;j<flamelets_library.number_of_species();j++)
			omega_[j].internalField()[celli] = extracted[j+1];
	}

    forAll(csi_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
        const fvPatchScalarField& pcsiv2 	= csiv2_.boundaryField()[patchi];
        const fvPatchScalarField& pchi_st 	= chi_st_.boundaryField()[patchi];
        const fvPatchScalarField& ph		= H_.boundaryField()[patchi];

        forAll(pcsi, facei)
        {
        
        	double max_chi = max(small_chi_st, pchi_st[facei]);
        
        	if (adiabaticMode == false)	
    			defect = ph[facei] - (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
    		
			if (pcsi[facei]<=0.)		// Pure oxidizer
			{
				flamelets_library.ExtractMeanValues(0., 0., max_chi, defect, extracted);
			}			
			
			else if (pcsi[facei]>=1.)	// Pure fuel	
			{
				flamelets_library.ExtractMeanValues(1., 0., max_chi, defect, extracted);
			}			
			
			else						// Mixture
			{
				csiv2_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));	// Normalized mixture fraction variance
			
				if (csiv2_normalized >= 0.98)
					flamelets_library.ExtractMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
				else if (csiv2_normalized < 0.)
					flamelets_library.ExtractMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
				else
					flamelets_library.ExtractMeanValues(pcsi[facei], csiv2_normalized, max_chi,defect, extracted);	
			}			
        
			for(int j=0;j<flamelets_library.number_of_species();j++)
				omega_[j].boundaryField()[patchi][facei] = extracted[j+1];
    	}
	}
	
}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::updateSourceTerms_MeanClosure_2E()
{
	const scalarField& TCells 		= T_.internalField();
	const scalarField& pCells 		= p_.internalField();
	const scalarField& rhoCells 		= density_reynolds_.internalField();
	scalarField& sootPhiNCells		= sootPhi_[0].internalField();
	scalarField& sootPhiMCells		= sootPhi_[1].internalField();
	const scalarField& sootC2H2Cells	= omega_[soot_index_c2h2].internalField();
	const scalarField& sootH2Cells		= omega_[soot_index_h2].internalField();
	const scalarField& sootO2Cells		= omega_[soot_index_o2].internalField();
	const scalarField& sootOHCells		= omega_[soot_index_oh].internalField();
	const scalarField& sootCOCells		= omega_[soot_index_co].internalField();
	const scalarField& sootHCells		= omega_[soot_index_h].internalField();

	scalarField& SnCells	= sootSource_[0].internalField();
	scalarField& SmCells	= sootSource_[1].internalField();
	
	scalarField& m0Cells	= sootProperties_[0].internalField();
	scalarField& fvCells	= sootProperties_[1].internalField();
	scalarField& dSootCells	= sootProperties_[2].internalField();
	scalarField& MSootCells	= sootProperties_[3].internalField();
	scalarField& ASootCells	= sootProperties_[4].internalField();


	std::vector<double> omega_soot(6);
	std::vector<double> Sgas(6);
 	
 	forAll(TCells, celli)
    	{	
		omega_soot[0] = sootC2H2Cells[celli];
		omega_soot[1] = sootH2Cells[celli];
		omega_soot[2] = sootO2Cells[celli];
		omega_soot[3] = sootOHCells[celli];
		omega_soot[4] = sootCOCells[celli];
		omega_soot[5] = sootHCells[celli];

		     if (sootPhiNCells[celli] <= 0.)		{	sootPhiNCells[celli] = soot_2e->phiNStart(1.0); sootPhiMCells[celli] = soot_2e->phiMStart(1.0);}
		else if (sootPhiMCells[celli] <= 0.)		{	sootPhiNCells[celli] = soot_2e->phiNStart(1.0); sootPhiMCells[celli] = soot_2e->phiMStart(1.0);}


		soot_2e->update(TCells[celli], pCells[celli]/101325., rhoCells[celli], omega_soot, sootPhiNCells[celli], sootPhiMCells[celli], SnCells[celli], SmCells[celli], Sgas);
		soot_2e->update_properties(rhoCells[celli], sootPhiNCells[celli], sootPhiMCells[celli], m0Cells[celli], fvCells[celli], dSootCells[celli], MSootCells[celli], ASootCells[celli]);

		if (celli == 100)
		{
			Info << TCells[celli] << " " << pCells[celli]/101325. << " " << rhoCells[celli] << " " << sootPhiNCells[celli] << " " << sootPhiMCells[celli] << " " << SnCells[celli] << " " << SmCells[celli] << endl;
			Info << omega_soot[2] << " " << omega_soot[3] << endl;
		}
	}

    	forAll(csi_.boundaryField(), patchi)
   	{
		const fvPatchScalarField& pt 		= T_.boundaryField()[patchi];
		const fvPatchScalarField& pp 		= p_.boundaryField()[patchi];
		const fvPatchScalarField& prho		= density_reynolds_.boundaryField()[patchi];
		fvPatchScalarField& psootPhiN 	= sootPhi_[0].boundaryField()[patchi];
		fvPatchScalarField& psootPhiM 	= sootPhi_[1].boundaryField()[patchi];
		const fvPatchScalarField& psootC2H2 	= omega_[soot_index_c2h2].boundaryField()[patchi];
		const fvPatchScalarField& psootH2 	= omega_[soot_index_h2].boundaryField()[patchi];
		const fvPatchScalarField& psootO2 	= omega_[soot_index_o2].boundaryField()[patchi];
		const fvPatchScalarField& psootOH 	= omega_[soot_index_oh].boundaryField()[patchi];
		const fvPatchScalarField& psootCO 	= omega_[soot_index_co].boundaryField()[patchi];
		const fvPatchScalarField& psootH 	= omega_[soot_index_h].boundaryField()[patchi];

		fvPatchScalarField& pSn	= sootSource_[0].boundaryField()[patchi];
		fvPatchScalarField& pSm	= sootSource_[1].boundaryField()[patchi];

		fvPatchScalarField& pm0		= sootProperties_[0].boundaryField()[patchi];
		fvPatchScalarField& pfv		= sootProperties_[1].boundaryField()[patchi];
		fvPatchScalarField& pdSoot	= sootProperties_[2].boundaryField()[patchi];
		fvPatchScalarField& pMSoot	= sootProperties_[3].boundaryField()[patchi];
		fvPatchScalarField& pASoot	= sootProperties_[4].boundaryField()[patchi];

		forAll(pt, facei)
		{
			omega_soot[0] = psootC2H2[facei];
			omega_soot[1] = psootH2[facei];
			omega_soot[2] = psootO2[facei];
			omega_soot[3] = psootOH[facei];
			omega_soot[4] = psootCO[facei];
			omega_soot[5] = psootH[facei];

		     	if (psootPhiN[facei] <= 0.)			{	psootPhiN[facei] = soot_2e->phiNStart(1.0); psootPhiM[facei] = soot_2e->phiMStart(1.0);}
			else if (psootPhiM[facei] <= 0.)		{	psootPhiN[facei] = soot_2e->phiNStart(1.0); psootPhiM[facei] = soot_2e->phiMStart(1.0);}


			soot_2e->update(pt[facei], pp[facei]/101325., prho[facei], omega_soot, psootPhiN[facei], psootPhiM[facei], pSn[facei], pSm[facei], Sgas);
			soot_2e->update_properties(prho[facei], psootPhiN[facei], psootPhiM[facei], pm0[facei], pfv[facei], pdSoot[facei], pMSoot[facei], pASoot[facei]);			 
		}
	}
	
}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::updateSourceTerms_MeanClosure_QMOM()
{
	const scalarField& TCells 		= T_.internalField();
	const scalarField& pCells 		= p_.internalField();
	const scalarField& rhoCells 		= density_reynolds_.internalField();
	const scalarField& muCells 		= mu_favre_.internalField();
	const scalarField& sootC2H2Cells	= omega_[soot_index_c2h2].internalField();
	const scalarField& sootH2Cells		= omega_[soot_index_h2].internalField();
	const scalarField& sootO2Cells		= omega_[soot_index_o2].internalField();
	const scalarField& sootOHCells		= omega_[soot_index_oh].internalField();
	const scalarField& sootCOCells		= omega_[soot_index_co].internalField();
	const scalarField& sootHCells		= omega_[soot_index_h].internalField();
	
	
	scalarField& m0Cells	= sootProperties_[0].internalField();
	scalarField& fvCells	= sootProperties_[1].internalField();
	scalarField& dSootCells	= sootProperties_[2].internalField();
	scalarField& MSootCells	= sootProperties_[3].internalField();
	scalarField& ASootCells	= sootProperties_[4].internalField();

 
	std::vector<double> omega_soot(6);
	double *moments = new double[2*nDirac];
	double *S 	= new double[2*nDirac];

 	forAll(TCells, celli)
    	{	
		omega_soot[0] = sootC2H2Cells[celli];
		omega_soot[1] = sootH2Cells[celli];
		omega_soot[2] = sootO2Cells[celli];
		omega_soot[3] = sootOHCells[celli];
		omega_soot[4] = sootCOCells[celli];
		omega_soot[5] = sootHCells[celli];

		for(int j=0;j<2*nDirac;j++)	moments[j] = sootMoments_[j].internalField()[celli];
		
		for(int j=0;j<2*nDirac;j++)
			if (moments[j]<=0.)
			{
			//	for(int k=0;k<2*nDirac;k++)
					moments[j] = 0.01; 
			//	break;
			}	

		soot_qmom->update(TCells[celli], pCells[celli]/101325., rhoCells[celli], muCells[celli], omega_soot, moments, S);

		soot_qmom->update_properties(	m0Cells[celli], fvCells[celli], dSootCells[celli], 
						MSootCells[celli], ASootCells[celli]);

		for(int j=0;j<2*nDirac;j++)	sootMomentsSources_[j].internalField()[celli] = S[j];

		if (celli == 100)
		{
			Info << "A " << TCells[celli] << " " << pCells[celli]/101325. << " " << rhoCells[celli] << " " << muCells[celli] << " " << moments[0] << " " << S[0] << moments[3] << " " << S[3] << endl;
			Info << "B " << m0Cells[celli] << " " << fvCells[celli] << " " << dSootCells[celli] << endl;
		}
	}

    	forAll(csi_.boundaryField(), patchi)
   	{
		const fvPatchScalarField& pt 		= T_.boundaryField()[patchi];
		const fvPatchScalarField& pp 		= p_.boundaryField()[patchi];
		const fvPatchScalarField& prho		= density_reynolds_.boundaryField()[patchi];
		const fvPatchScalarField& pmu		= mu_favre_.boundaryField()[patchi];
		const fvPatchScalarField& psootC2H2 	= omega_[soot_index_c2h2].boundaryField()[patchi];
		const fvPatchScalarField& psootH2 	= omega_[soot_index_h2].boundaryField()[patchi];
		const fvPatchScalarField& psootO2 	= omega_[soot_index_o2].boundaryField()[patchi];
		const fvPatchScalarField& psootOH 	= omega_[soot_index_oh].boundaryField()[patchi];
		const fvPatchScalarField& psootCO 	= omega_[soot_index_co].boundaryField()[patchi];
		const fvPatchScalarField& psootH 	= omega_[soot_index_h].boundaryField()[patchi];

		fvPatchScalarField& pm0		= sootProperties_[0].boundaryField()[patchi];
		fvPatchScalarField& pfv		= sootProperties_[1].boundaryField()[patchi];
		fvPatchScalarField& pdSoot	= sootProperties_[2].boundaryField()[patchi];
		fvPatchScalarField& pMSoot	= sootProperties_[3].boundaryField()[patchi];
		fvPatchScalarField& pASoot	= sootProperties_[4].boundaryField()[patchi];

		forAll(pt, facei)
		{
			omega_soot[0] = psootC2H2[facei];
			omega_soot[1] = psootH2[facei];
			omega_soot[2] = psootO2[facei];
			omega_soot[3] = psootOH[facei];
			omega_soot[4] = psootCO[facei];
			omega_soot[5] = psootH[facei];

			for(int j=0;j<2*nDirac;j++)	moments[j] = sootMoments_[j].boundaryField()[patchi][facei];
		
			for(int j=0;j<2*nDirac;j++)
				if (moments[j]<=0.)
				{
				//	for(int k=0;k<2*nDirac;k++)
						moments[j] = 0.01; 
				//	break;
				}	

			soot_qmom->update(pt[facei], pp[facei]/101325., prho[facei], pmu[facei], omega_soot, moments, S);

			soot_qmom->update_properties(	pm0[facei], pfv[facei], pdSoot[facei], pMSoot[facei], pASoot[facei]);

			for(int j=0;j<2*nDirac;j++)	sootMomentsSources_[j].boundaryField()[patchi][facei] = S[j];		 
		}
	}
}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::updateSourceTerms_UncorrelatedClosure_2E()
{	
	double csiv2_normalized = 0.;
	double small_eps 	= 1.e-6;
	double small_chi_st 	= 1.e-8;

	const scalarField& csiCells 	= csi_.internalField();
	const scalarField& csiv2Cells 	= csiv2_.internalField();
	const scalarField& chi_stCells 	= chi_st_.internalField();	
	const scalarField& rhoCells 	= density_reynolds_.internalField();
	scalarField& sootPhiNCells	= sootPhi_[0].internalField();
	scalarField& sootPhiMCells	= sootPhi_[1].internalField();

	scalarField& SnCells	= sootSource_[0].internalField();
	scalarField& SmCells	= sootSource_[1].internalField();
	
	scalarField& m0Cells	= sootProperties_[0].internalField();
	scalarField& fvCells	= sootProperties_[1].internalField();
	scalarField& dSootCells	= sootProperties_[2].internalField();
	scalarField& MSootCells	= sootProperties_[3].internalField();
	scalarField& ASootCells	= sootProperties_[4].internalField();


 	forAll(csiCells, celli)
    	{	
		double max_chi = max(small_chi_st,chi_stCells[celli]);

		     if (sootPhiNCells[celli] <= 0.)		{	sootPhiNCells[celli] = soot_2e->phiNStart(1.0); sootPhiMCells[celli] = soot_2e->phiMStart(1.0);}
		else if (sootPhiMCells[celli] <= 0.)		{	sootPhiNCells[celli] = soot_2e->phiNStart(1.0); sootPhiMCells[celli] = soot_2e->phiMStart(1.0);}

		if (csiCells[celli]<=small_eps)			// Pure oxidizer
			soot_2e_uncorrelated->GetMeanValues(0., 0., max_chi, sootPhiNCells[celli], sootPhiMCells[celli], rhoCells[celli],  SnCells[celli], SmCells[celli]);

		else if (csiCells[celli]>=(1.-small_eps))	// Pure fuel	
			soot_2e_uncorrelated->GetMeanValues(1., 0., max_chi, sootPhiNCells[celli], sootPhiMCells[celli], rhoCells[celli],  SnCells[celli], SmCells[celli]);			

		else									// Mixture
		{
			csiv2_normalized = csiv2Cells[celli] / (csiCells[celli]*(1.-csiCells[celli]));	// Normalized mixture fraction variance
	
			if (csiv2_normalized >= 0.98)
				soot_2e_uncorrelated->GetMeanValues(csiCells[celli], 0.98, max_chi, sootPhiNCells[celli], sootPhiMCells[celli], rhoCells[celli],  SnCells[celli], SmCells[celli]);
			else if (csiv2_normalized < 0.)
				soot_2e_uncorrelated->GetMeanValues(csiCells[celli], 0.00, max_chi, sootPhiNCells[celli], sootPhiMCells[celli], rhoCells[celli],  SnCells[celli], SmCells[celli]);
			else
				soot_2e_uncorrelated->GetMeanValues(csiCells[celli], csiv2Cells[celli], max_chi, sootPhiNCells[celli], sootPhiMCells[celli], rhoCells[celli],  SnCells[celli], SmCells[celli]);
		}			

		// Update properties
		soot_2e_uncorrelated->UpdateProperties(rhoCells[celli], sootPhiNCells[celli], sootPhiMCells[celli], m0Cells[celli], fvCells[celli], dSootCells[celli], MSootCells[celli], ASootCells[celli]);
	}

    	forAll(csi_.boundaryField(), patchi)
   	{
		const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
		const fvPatchScalarField& pcsiv2 	= csiv2_.boundaryField()[patchi];
		const fvPatchScalarField& pchi_st 	= chi_st_.boundaryField()[patchi];
		const fvPatchScalarField& prho		= density_reynolds_.boundaryField()[patchi];
		fvPatchScalarField& psootPhiN 		= sootPhi_[0].boundaryField()[patchi];
		fvPatchScalarField& psootPhiM 		= sootPhi_[1].boundaryField()[patchi];
		
		fvPatchScalarField& pSn	= sootSource_[0].boundaryField()[patchi];
		fvPatchScalarField& pSm	= sootSource_[1].boundaryField()[patchi];

		fvPatchScalarField& pm0		= sootProperties_[0].boundaryField()[patchi];
		fvPatchScalarField& pfv		= sootProperties_[1].boundaryField()[patchi];
		fvPatchScalarField& pdSoot	= sootProperties_[2].boundaryField()[patchi];
		fvPatchScalarField& pMSoot	= sootProperties_[3].boundaryField()[patchi];
		fvPatchScalarField& pASoot	= sootProperties_[4].boundaryField()[patchi];

		forAll(pcsi, facei)
		{
			double max_chi = max(small_chi_st,pchi_st[facei]);

		     		if (psootPhiN[facei] <= 0.)		{	psootPhiN[facei] = soot_2e->phiNStart(1.0); psootPhiM[facei] = soot_2e->phiMStart(1.0);}
			else 	if (psootPhiM[facei] <= 0.)		{	psootPhiN[facei] = soot_2e->phiNStart(1.0); psootPhiM[facei] = soot_2e->phiMStart(1.0);}

			if (pcsi[facei]<=small_eps)		// Pure oxidizer
				soot_2e_uncorrelated->GetMeanValues(0.0, 0.0, max_chi, psootPhiN[facei], psootPhiM[facei], prho[facei],  pSn[facei], pSm[facei]);

			else if (pcsi[facei]>=(1.-small_eps))	// Pure fuel	
				soot_2e_uncorrelated->GetMeanValues(1.0, 0.0, max_chi, psootPhiN[facei], psootPhiM[facei], prho[facei],  pSn[facei], pSm[facei]);			

			else						// Mixture
			{
				csiv2_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));	// Normalized mixture fraction variance
	
				if (csiv2_normalized >= 0.98)
					soot_2e_uncorrelated->GetMeanValues(pcsi[facei], 0.98, max_chi, psootPhiN[facei], psootPhiM[facei], prho[facei],  pSn[facei], pSm[facei]);
				else if (csiv2_normalized < 0.)
					soot_2e_uncorrelated->GetMeanValues(pcsi[facei], 0., max_chi, psootPhiN[facei], psootPhiM[facei], prho[facei],  pSn[facei], pSm[facei]);
				else
					soot_2e_uncorrelated->GetMeanValues(pcsi[facei], pcsiv2[facei], max_chi, psootPhiN[facei], psootPhiM[facei], prho[facei],  pSn[facei], pSm[facei]);
			}				

			// Update properties
			soot_2e_uncorrelated->UpdateProperties(prho[facei], psootPhiN[facei], psootPhiM[facei], pm0[facei], pfv[facei], pdSoot[facei], pMSoot[facei], pASoot[facei]);			 
		}
	}
	
}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::updateSourceTerms_CorrelatedClosure_2E()
{
}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::ErrorMessage(const string message)
{
	Info << "Class: pdfThermo" << endl;
	Info << "Error: " << message << endl;
	getchar();
	//exit(-1);
}


// ************************************************************************* //
