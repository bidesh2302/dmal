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
phi=(0.5*(3.0*phi.oldTime()-phi.oldTime().oldTime()));
UEqn =
(
	fvm::ddt(rho,U)
	+ fvm::div(phi, U)
	+ turb.divDevRhoReff(U)
 ==
	fvc::reconstruct
	(
		(
		- ghf*fvc::snGrad(rho)
		- fvc::snGrad(p_rgh)
		)*mesh.magSf()
	)
);

solve(UEqn == -fvc::grad(p_rgh));

dimensionedScalar dt=runTime.deltaT();

U += dt*fvc::grad(p_rgh/rho);

U.correctBoundaryConditions();

phi = (fvc::interpolate(rho*U) & mesh.Sf());
adjustPhi(phi, U, p_rgh);

fvScalarMatrix pEqn
(
	fvm::laplacian(dt, p_rgh) == fvc::div(phi)
);

pEqn.setReference(pRefCell, pRefValue);

pEqn.solve();

#include "continuityErrs.H"

U -= dt*fvc::grad(p_rgh/rho);

U.correctBoundaryConditions();

p = p_rgh + rho*gh;

p.correctBoundaryConditions();

//~ dpdt = fvc::ddt(p);

p_rgh = p -  rho*gh;

