if (iTwoEquationModel == true)
{

	//  Solving equation of phiN
	fvScalarMatrix phiN_Eqn
	(
		(
			fvm::div(phi, phiN)
	      	      - fvm::Sp(fvc::div(phi), phiN)		
              	      - fvm::laplacian(mut/sigmat, phiN)
		      - Sn
		)
    	);
	phiN_Eqn.relax();
	phiN_Eqn.solve(mesh.solver("phiN"));
	   

	//  Solving equation of phiM
	fvScalarMatrix phiM_Eqn
	(
		(
			fvm::div(phi, phiM)
	      	      - fvm::Sp(fvc::div(phi), phiM)		
              	      - fvm::laplacian(mut/sigmat, phiM)
		      - Sm
		)
    	);
	phiM_Eqn.relax();
	phiM_Eqn.solve(mesh.solver("phiM"));

}
