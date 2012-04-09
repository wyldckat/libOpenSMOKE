{	
	if (adiabaticMode == true)
	{
		fvScalarMatrix hEqn
		(
			fvm::div(phi, H)
	      	      - fvm::Sp(fvc::div(phi), H)		
                      - fvm::laplacian(turbulence->mut()/sigmat, H)
		);
		
		hEqn.relax();
		hEqn.solve(mesh.solver("H"));
		
		Qrad = 4. * sigmaSB * as * (pow(thermo.T(),4) - pow(Tenv,4));
	}
	else
	{
		
		Qrad = 4. * sigmaSB * as * (pow(thermo.T(),4) - pow(Tenv,4));
		
		fvScalarMatrix hEqn
		(
			fvm::div(phi, H)
	      		- fvm::Sp(fvc::div(phi), H)		
          		- fvm::laplacian(turbulence->mut()/sigmat, H) + Qrad
		);
		
		hEqn.relax();
		hEqn.solve(mesh.solver("H"));
	}
	
	Info << "Correct thermodynamics" << endl;
	thermo.correct();
}
