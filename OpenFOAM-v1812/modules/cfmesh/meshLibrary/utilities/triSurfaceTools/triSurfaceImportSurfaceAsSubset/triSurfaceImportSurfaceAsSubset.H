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

Class
    Foam::Module::triSurfaceImportSurfaceAsSubset

Description
    Creates a subset in the master surface consisting of facets which are
    near the other surface

SourceFiles
    triSurfaceImportSurfaceAsSubset.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceImportSurfaceAsSubset_H
#define triSurfaceImportSurfaceAsSubset_H

#include "triSurf.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

class meshOctree;

/*---------------------------------------------------------------------------*\
               Class triSurfaceImportSurfaceAsSubset Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceImportSurfaceAsSubset
{
    // Private data

        //- reference to triSurf
        triSurf& surf_;

        //- pointer to meshOctree, needed for searching on the master surface
        meshOctree* octreePtr_;


    // Private member functions

        void createOctree(const triSurf&, meshOctree&);

        //- Disallow default bitwise copy construct
        triSurfaceImportSurfaceAsSubset(const triSurfaceImportSurfaceAsSubset&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceImportSurfaceAsSubset&);


public:

    //- Construct from triSurf
    triSurfaceImportSurfaceAsSubset(triSurf& surface);

    //- Destructor
    ~triSurfaceImportSurfaceAsSubset();


    // Member Functions

        //- finds the nearest faces in the surface to the import surf
        //  and creates a subset
        void addSurfaceAsSubset
        (
            const triSurf& importSurf,
            const word& subsetName,
            const scalar angleTol = 5.*M_PI/180.
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
