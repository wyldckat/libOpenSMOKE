	// Soot closures	
	enum {SOOT_2E_NONE, SOOT_2E_MEAN, SOOT_2E_UNCORRELATED, SOOT_2E_CORRELATED} 	iModelClosure;

	// Soot classes
	OpenSMOKE_SootEmpiricalModels			*soot_2e;
	OpenSMOKE_SootSourceTermLibrary 		*soot_2e_uncorrelated;

	// Soot scalar fields
	PtrList<volScalarField> sootProperties_;
	PtrList<volScalarField> omegaGas_;
	PtrList<volScalarField> flameletsData_;


	// Soot dictionary
	IOdictionary sootPostProcessing_
	(
		IOobject
		(
		    "sootPostProcessing",
		    U.time().constant(),
		    U.db(),
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		)
	);

	// Prandtl number for mixture fraction, variance of mixture fraction and enthalpy diffusion terms
	scalar sigmat_ = readScalar(sootPostProcessing_.lookup("sigmat"));
	dimensionedScalar sigmat("sigmat", dimensionSet(0,0,0,0,0,0,0), sigmat_);

	// Soot robust
	Switch sootRobust(sootPostProcessing_.lookup("sootRobust"));
	
	// Soot update
	label sootUpdate = readLabel(sootPostProcessing_.lookup("sootUpdate"));

	// Closure (turbulent fluctuations)
	Foam::string soot_closure = sootPostProcessing_.lookup("sootClosure");
	
							iModelClosure = SOOT_2E_NONE;
	if 	(soot_closure == "Mean")		iModelClosure = SOOT_2E_MEAN;
	else if (soot_closure == "Uncorrelated")	iModelClosure = SOOT_2E_UNCORRELATED;
	else if (soot_closure == "Correlated")		iModelClosure = SOOT_2E_CORRELATED;
	else ErrorMessage("Wrong soot model closure...");

	
	// Soot properties
	{	
		std::vector<std::string> soot_properties_names(5);
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
	}


	nucleation_models	nucleation_model_ 		= NUCLEATION_NONE;
	growth_models		growth_model_			= GROWTH_NONE;
	aggregation_models	aggregation_model_		= AGGREGATION_NONE;
	oxidation_models	oxidation_model_		= OXIDATION_NONE;

	{
		Foam::string nucleation_model 	= sootPostProcessing_.lookup("nucleationModel");
		Foam::string growth_model 	= sootPostProcessing_.lookup("growthModel");
		Foam::string aggregation_model 	= sootPostProcessing_.lookup("aggregationModel");
		Foam::string oxidation_model	= sootPostProcessing_.lookup("oxidationModel");
			
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
	}

	if (iModelClosure == SOOT_2E_MEAN)
	{
		soot_2e = new OpenSMOKE_SootEmpiricalModels(nucleation_model_, growth_model_, oxidation_model_, aggregation_model_);
		
		if (sootRobust == true)	
			soot_2e->setRobustCalculations();	

		// Gas mass fractions
		std::vector<std::string> omega_gas_names(6);
		omega_gas_names[0] = "omega_C2H2";
		omega_gas_names[1] = "omega_H2";
		omega_gas_names[2] = "omega_O2";
		omega_gas_names[3] = "omega_OH";
		omega_gas_names[4] = "omega_CO";
		omega_gas_names[5] = "omega_H";

		omegaGas_.setSize(6);
		for (int j=0;j<6;j++)
			omegaGas_.set(	j, new volScalarField( 	IOobject(
										omega_gas_names[j], 
										mesh.time().timeName(), 
										mesh, 
										IOobject::MUST_READ, 
										IOobject::AUTO_WRITE),
										mesh
									) );
	
	}

	else if (iModelClosure == SOOT_2E_UNCORRELATED)
	{
		Foam::string sootLibraryPath = sootPostProcessing_.lookup("libraryPath");

		soot_2e_uncorrelated = new OpenSMOKE_SootSourceTermLibrary();
		
		soot_2e_uncorrelated->SetNucleationModel(nucleation_model_);
		soot_2e_uncorrelated->SetGrowthModel(growth_model_);
		soot_2e_uncorrelated->SetAggregationModel(aggregation_model_);
		soot_2e_uncorrelated->SetOxidationModel(oxidation_model_);

		soot_2e_uncorrelated->SetLibraryPath(sootLibraryPath);
		if (sootRobust == true)	
			soot_2e_uncorrelated->setRobustCalculations();

		//soot_2e_uncorrelated->SetLogNormalChiDistribution();
		//soot_2e_uncorrelated->SetNoFluctuationsExtractionMode();

		soot_2e_uncorrelated->Setup();

		// Flamelets data names
		std::vector<std::string> flamelets_data_names(3);
		flamelets_data_names[0] = "csi";
		flamelets_data_names[1] = "csiv2";
		flamelets_data_names[2] = "chi_st";

		flameletsData_.setSize(3);
		for (int j=0;j<3;j++)
			flameletsData_.set(	j, new volScalarField( 	IOobject(
											flamelets_data_names[j], 
											mesh.time().timeName(), 
											mesh, 
											IOobject::MUST_READ, 
											IOobject::AUTO_WRITE),
											mesh
										) );
	}

	int soot_counter = sootUpdate;

