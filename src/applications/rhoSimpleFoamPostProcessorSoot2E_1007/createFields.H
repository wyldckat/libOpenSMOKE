	Info<< "Preparing field rho (density)..." << endl;
	volScalarField rho
	(
		IOobject
		(
		    "rho",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Reading field U (velocity)..." << endl;
	volVectorField U
	(
		IOobject
		(
		    "U",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Preparing field T (temperature)..." << endl;
	volScalarField T	
	(
		IOobject
		(
		    "T",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Preparing field p (pressure)..." << endl;
	volScalarField p
	(
		IOobject
		(
		    "p",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Preparing field mut (turbulent viscosity)..." << endl;
	volScalarField mut
	(
		IOobject
		(
		    "mut",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

#   	include "compressibleCreatePhi.H"


	bool iTwoEquationModel 	= true;

	Info<< "Preparing field phiN (number of soot particles per unit volume)..." << endl;
	volScalarField phiN
	(
		IOobject
		(
		    "phiN",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Preparing field phiM (soot mass fraction)..." << endl;
	volScalarField phiM
	(
		IOobject
		(
		    "phiM",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Preparing field Sn (source of number of soot particles per unit volume)..." << endl;
	volScalarField Sn
	(
		IOobject
		(
		    "Sn",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		mesh,  
		dimensionedScalar("zero", dimensionSet(0,-3,-1,0,1,0,0), 0.0) 
	);

	Info<< "Preparing field Sm (source of soot mass fraction)..." << endl;
	volScalarField Sm
	(
		IOobject
		(
		    "Sm",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		mesh,  
		dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
	);

/*	Info<< "Preparing field omega_C2H2..." << endl;
	volScalarField omega_C2H2
	(
		IOobject
		(
		    "omega_C2H2",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Preparing field omega_H2..." << endl;
	volScalarField omega_H2
	(
		IOobject
		(
		    "omega_H2",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Preparing field omega_O2..." << endl;
	volScalarField omega_O2
	(
		IOobject
		(
		    "omega_O2",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Preparing field omega_OH..." << endl;
	volScalarField omega_OH
	(
		IOobject
		(
		    "omega_OH",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Preparing field omega_CO..." << endl;
	volScalarField omega_CO
	(
		IOobject
		(
		    "omega_CO",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Preparing field omega_H..." << endl;
	volScalarField omega_H
	(
		IOobject
		(
		    "omega_H",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	);*/


