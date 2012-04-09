if (soot_counter == sootUpdate)
{	
	if (iModelClosure == SOOT_2E_MEAN)
	{
		Info << "Updating soot source terms (mean closure)..." << endl;

		std::vector<double> omega_soot(6);
		std::vector<double> Sgas(6);
	
		{
			const scalarField& TCells 		= T.internalField();
			const scalarField& pCells 		= p.internalField();
			const scalarField& rhoCells 		= rho.internalField();
			      scalarField& sootPhiNCells	= phiN.internalField();
			      scalarField& sootPhiMCells	= phiM.internalField();
	
			const scalarField& sootC2H2Cells	= omegaGas_[0].internalField();
			const scalarField& sootH2Cells		= omegaGas_[1].internalField();
			const scalarField& sootO2Cells		= omegaGas_[2].internalField();
			const scalarField& sootOHCells		= omegaGas_[3].internalField();
			const scalarField& sootCOCells		= omegaGas_[4].internalField();
			const scalarField& sootHCells		= omegaGas_[5].internalField();

			scalarField& SnCells	= Sn.internalField();
			scalarField& SmCells	= Sm.internalField();
	
			scalarField& m0Cells	= sootProperties_[0].internalField();
			scalarField& fvCells	= sootProperties_[1].internalField();
			scalarField& dSootCells	= sootProperties_[2].internalField();
			scalarField& MSootCells	= sootProperties_[3].internalField();
			scalarField& ASootCells	= sootProperties_[4].internalField();

		 	
		 	forAll(TCells, celli)
		    	{	
				omega_soot[0] = sootC2H2Cells[celli];
				omega_soot[1] = sootH2Cells[celli];
				omega_soot[2] = sootO2Cells[celli];
				omega_soot[3] = sootOHCells[celli];
				omega_soot[4] = sootCOCells[celli];
				omega_soot[5] = sootHCells[celli];

				     if (sootPhiNCells[celli] <= 0.)		{	sootPhiNCells[celli] = soot_2e->phiNStart(1.0); sootPhiMCells[celli] = soot_2e->phiMStart(1.0);	}
				else if (sootPhiMCells[celli] <= 0.)		{	sootPhiNCells[celli] = soot_2e->phiNStart(1.0); sootPhiMCells[celli] = soot_2e->phiMStart(1.0);	}


				soot_2e->update(TCells[celli], pCells[celli]/101325., rhoCells[celli], omega_soot, sootPhiNCells[celli], sootPhiMCells[celli], SnCells[celli], SmCells[celli], Sgas);
				soot_2e->update_properties(rhoCells[celli], sootPhiNCells[celli], sootPhiMCells[celli], m0Cells[celli], fvCells[celli], dSootCells[celli], MSootCells[celli], ASootCells[celli]);

				if (celli == 100)
					Info << TCells[celli] << " " << rhoCells[celli] << " " << sootPhiNCells[celli] << " " << sootPhiMCells[celli] << " " << SnCells[celli] << " " << SmCells[celli] << endl;
			}
		}

	    	forAll(T.boundaryField(), patchi)
	   	{
			const fvPatchScalarField& pt 		= T.boundaryField()[patchi];
			const fvPatchScalarField& pp 		= p.boundaryField()[patchi];
			const fvPatchScalarField& prho		= rho.boundaryField()[patchi];
			      fvPatchScalarField& psootPhiN 	= phiN.boundaryField()[patchi];
			      fvPatchScalarField& psootPhiM 	= phiM.boundaryField()[patchi];
		
			const fvPatchScalarField& psootC2H2 	= omegaGas_[0].boundaryField()[patchi];
			const fvPatchScalarField& psootH2 	= omegaGas_[1].boundaryField()[patchi];
			const fvPatchScalarField& psootO2 	= omegaGas_[2].boundaryField()[patchi];
			const fvPatchScalarField& psootOH 	= omegaGas_[3].boundaryField()[patchi];
			const fvPatchScalarField& psootCO 	= omegaGas_[4].boundaryField()[patchi];
			const fvPatchScalarField& psootH 	= omegaGas_[5].boundaryField()[patchi];

			fvPatchScalarField& pSn	= Sn.boundaryField()[patchi];
			fvPatchScalarField& pSm	= Sm.boundaryField()[patchi];

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

	else if (iModelClosure == SOOT_2E_UNCORRELATED)
	{
		Info << "Updating soot source terms (uncorrelated closure)..." << endl;

		double csiv2_normalized = 0.;
		double small_eps 	= 1.e-6;
		double small_chi_st 	= 1.e-8;

		{
			const scalarField& csiCells 		= flameletsData_[0].internalField();
			const scalarField& csiv2Cells 		= flameletsData_[1].internalField();
			const scalarField& chi_stCells 		= flameletsData_[2].internalField();	
			const scalarField& rhoCells 		= rho.internalField();
			      scalarField& sootPhiNCells	= phiN.internalField();
			      scalarField& sootPhiMCells	= phiM.internalField();

			scalarField& SnCells	= Sn.internalField();
			scalarField& SmCells	= Sm.internalField();
	
			scalarField& m0Cells	= sootProperties_[0].internalField();
			scalarField& fvCells	= sootProperties_[1].internalField();
			scalarField& dSootCells	= sootProperties_[2].internalField();
			scalarField& MSootCells	= sootProperties_[3].internalField();
			scalarField& ASootCells	= sootProperties_[4].internalField();

		 	forAll(csiCells, celli)
		    	{	
				double max_chi = max(small_chi_st,chi_stCells[celli]);

				     if (sootPhiNCells[celli] <= 0.)		{	sootPhiNCells[celli] = 2e-11; sootPhiMCells[celli] = 2.3e-15; }
				else if (sootPhiMCells[celli] <= 0.)		{	sootPhiNCells[celli] = 2e-11; sootPhiMCells[celli] = 2.3e-15; }

				if (csiCells[celli]<=small_eps)			// Pure oxidizer
					soot_2e_uncorrelated->GetMeanValues(0., 0., max_chi, sootPhiNCells[celli], sootPhiMCells[celli], rhoCells[celli],  SnCells[celli], SmCells[celli]);

				else if (csiCells[celli]>=(1.-small_eps))	// Pure fuel	
					soot_2e_uncorrelated->GetMeanValues(1., 0., max_chi, sootPhiNCells[celli], sootPhiMCells[celli], rhoCells[celli],  SnCells[celli], SmCells[celli]);			

				else						// Mixture
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
		}

		forAll(flameletsData_[0].boundaryField(), patchi)
		{
			const fvPatchScalarField& pcsi 		= flameletsData_[0].boundaryField()[patchi];
			const fvPatchScalarField& pcsiv2 	= flameletsData_[1].boundaryField()[patchi];
			const fvPatchScalarField& pchi_st 	= flameletsData_[2].boundaryField()[patchi];
			const fvPatchScalarField& prho		= rho.boundaryField()[patchi];
			      fvPatchScalarField& psootPhiN 	= phiN.boundaryField()[patchi];
			      fvPatchScalarField& psootPhiM 	= phiM.boundaryField()[patchi];

			fvPatchScalarField& pSn	= Sn.boundaryField()[patchi];
			fvPatchScalarField& pSm	= Sm.boundaryField()[patchi];

			fvPatchScalarField& pm0		= sootProperties_[0].boundaryField()[patchi];
			fvPatchScalarField& pfv		= sootProperties_[1].boundaryField()[patchi];
			fvPatchScalarField& pdSoot	= sootProperties_[2].boundaryField()[patchi];
			fvPatchScalarField& pMSoot	= sootProperties_[3].boundaryField()[patchi];
			fvPatchScalarField& pASoot	= sootProperties_[4].boundaryField()[patchi];

			forAll(pcsi, facei)
			{
				double max_chi = max(small_chi_st,pchi_st[facei]);

					if (psootPhiN[facei] <= 0.)		{	psootPhiN[facei] = 2e-11; psootPhiM[facei] = 2.3e-15; }
				else 	if (psootPhiM[facei] <= 0.)		{	psootPhiN[facei] = 2e-11; psootPhiM[facei] = 2.3e-15; }

				if (pcsi[facei]<=small_eps)		// Pure oxidizer
					soot_2e_uncorrelated->GetMeanValues(0.0, 0.0, max_chi, psootPhiN[facei], psootPhiM[facei], prho[facei],  pSn[facei], pSm[facei]);

				else if (pcsi[facei]>=(1.-small_eps))	// Pure fuel	
					soot_2e_uncorrelated->GetMeanValues(1.0, 0.0, max_chi, psootPhiN[facei], psootPhiM[facei], prho[facei],  pSn[facei], pSm[facei]);			

				else					// Mixture
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

	soot_counter = 0;
}

soot_counter++;


