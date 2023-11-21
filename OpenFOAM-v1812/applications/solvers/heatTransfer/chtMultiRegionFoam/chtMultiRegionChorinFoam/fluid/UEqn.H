    // Solve the Momentum equation

    //~ MRF.correctBoundaryVelocity(U);

    //~ UEqn =
    //~ (
        //~ fvm::div(phi, U)
      //~ + MRF.DDt(rho, U)
      //~ + turb.divDevRhoReff(U)
     //~ ==
        //~ fvOptions(rho, U)
    //~ );

    //~ UEqn.relax();

    //~ fvOptions.constrain(UEqn);

    //~ if (momentumPredictor)
    //~ {
        //~ solve
        //~ (
            //~ UEqn
        //~ ==
            //~ fvc::reconstruct
            //~ (
                //~ (
                //~ - ghf*fvc::snGrad(rho)
                //~ - fvc::snGrad(p_rgh)
                //~ )*mesh.magSf()
            //~ )
        //~ );

        //~ fvOptions.correct(U);
    //~ }
//####################################################################//
{
            // Momentum predictor
			phi=(0.5*(3.0*phi.oldTime()-phi.oldTime().oldTime()));
			
            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              + turb.divDevReff()
            );

            solve(UEqn == -fvc::grad(p_rgh));
            U.correctBoundaryConditions();
            
            dimensionedScalar dt=runTime.deltaT();

			U += dt*fvc::grad(p_rgh);
			U.correctBoundaryConditions();
            
            phi = (fvc::interpolate(U) & mesh.Sf());
            adjustPhi(phi, U, p_rgh);
            
            fvScalarMatrix pEqn
			(
				fvm::laplacian(dt, p_rgh) == fvc::div(phi)
			);
    
			pEqn.setReference(pRefCell, pRefValue);
    
			pEqn.solve();

#           include "continuityErrs.H"

			U -= dt*fvc::grad(p_rgh);
			
			U.correctBoundaryConditions();
			p_rgh.correctBoundaryConditions();
			
			volScalarField kappaEff
			(
				"kappaEff",
				turbulence->nu()/Pr + turbulence->nut()/Prt
			);

			fvScalarMatrix TEqn
			(
				fvm::ddt(T)
				+ fvm::div(phi, T)
				- fvm::laplacian(kappaEff, T)
			);

			TEqn.relax();
			TEqn.solve();
			T.correctBoundaryConditions();
        }

