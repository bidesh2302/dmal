/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | www.cfmesh.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
     Copyright (C) 2014-2017 Creative Fields, Ltd.
-------------------------------------------------------------------------------
Author
     Franjo Juretic (franjo.juretic@c-fields.com)

License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Writes the mesh in fpma format readable by AVL's CfdWM

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMeshGenModifier.H"
#include "meshSurfaceOptimizer.H"
#include "triSurf.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "meshSurfaceEngine.H"
#include "labelledTri.H"

using namespace Foam;
using namespace Foam::Module;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    polyMeshGen pmg(runTime);
    pmg.read();

    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fileName surfaceFile = meshDict.lookup("surfaceFile");
    if (Pstream::parRun())
        surfaceFile = ".."/surfaceFile;

    triSurf surf(runTime.path()/surfaceFile);

    meshOctree mo(surf);

    meshOctreeCreator(mo, meshDict).createOctreeBoxes();

    meshSurfaceEngine mse(pmg);
    meshSurfaceOptimizer surfOpt(mse, mo);

    surfOpt.optimizeSurface(10);

    pmg.write();

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
