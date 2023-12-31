/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Global
    surfaceCourantNo

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Description
    Calculates and outputs the mean and maximum Courant Numbers for the
    Finite Area method.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
scalar velMag = 0.0;
if (aMesh.nInternalEdges())
{
    edgeScalarField edgeSpeed
    (
        mag(fac::interpolate(Us) & aMesh.Le()/aMesh.magLe())
      + Foam::sqrt(max(dimensionedScalar("gh_min", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0), fac::interpolate(R*c*h*boundaryCell)*fac::interpolate(gn)))
    );

    velMag = max(edgeSpeed.internalField()).value();

    CoNum = max(aMesh.edgeInterpolation::deltaCoeffs()*edgeSpeed.internalField()).value()*runTime.deltaT().value();

    meanCoNum = sum(aMesh.edgeInterpolation::deltaCoeffs()*edgeSpeed.internalField()*aMesh.magLe()).value()/gSum(aMesh.magLe())*runTime.deltaT().value();
}


Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum
    << " velocity max magnitude: " << velMag << endl;

// ************************************************************************* //
