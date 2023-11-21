/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoFoam

Description
    Transient solver for incompressible flow.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pisoControl.H"

//~ #include "fvIOoptionList.H"
//~ #include "IOporosityModelList.H"
//~ #include "IOMRFZoneList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pisoControl piso(mesh);

#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
    #include "CourantNo.H"							    
    #include "createTimeControls.H"						    
    #include "setInitialDeltaT.H"

    //~ #include "createFvOptions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"


        // Pressure-velocity PISO corrector
        {
            // Momentum predictor
			phi=(0.5*(3.0*phi.oldTime()-phi.oldTime().oldTime()));
			
            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              + turbulence->divDevReff()
              - fvc::reconstruct
				(
					fvc::interpolate(rhok)*(g & mesh.Sf())
				)
            );
            
            //~ volScalarField p_rgh = p + gh;
            
			solve(UEqn == -fvc::grad(p));
			//~ solve(UEqn == -fvc::grad(p_rgh));
            U.correctBoundaryConditions();
            
            dimensionedScalar dt=runTime.deltaT();

			U += dt*fvc::grad(p);
			U.correctBoundaryConditions();
            
            phi = (fvc::interpolate(U) & mesh.Sf());
            adjustPhi(phi, U, p);
            
            fvScalarMatrix pEqn
			(
				fvm::laplacian(dt, p) == fvc::div(phi)
				//~ fvm::laplacian(dt, p_rgh) == fvc::div(phi)
			);
			
			pEqn.setReference(pRefCell, pRefValue);
    
			pEqn.solve();
			
			//~ p = p_rgh - gh;
			
#           include "continuityErrs.H"

			U -= dt*fvc::grad(p);

			U.correctBoundaryConditions();
			p.correctBoundaryConditions();
			
#       	include "TEqn.H"
			T.correctBoundaryConditions();
        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
