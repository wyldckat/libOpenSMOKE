/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    rhoSimpleFoam

Description
    Steady-state solver for turbulent flow of compressible fluids

\*---------------------------------------------------------------------------*/

#include "OpenSMOKE_SootEmpiricalModels.h"
#include "OpenSMOKE_SootSourceTermLibrary.h"

#include "fvCFD.H"
#include "fixedGradientFvPatchFields.H"

#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void ErrorMessage(const Foam::string message);

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "sootTwoEquationsOptions.H"
#   include "initContinuityErrs.H"

    simpleControl simple(mesh);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        {
#		include "sootTwoEquationsUpdate.H"
#		include "sootTwoEquationsEqn.H"
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

void ErrorMessage(const Foam::string message)
{
	Info << "rhoSimpleFoamPostProcessor2E" << endl;
	Info << "Error: " << message << endl;
	getchar();
}

// ************************************************************************* //
