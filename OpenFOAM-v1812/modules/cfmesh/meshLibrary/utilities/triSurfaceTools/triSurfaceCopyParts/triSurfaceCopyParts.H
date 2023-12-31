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
    Foam::Module::triSurfaceCopyParts

Description
    Copies user-selected surface parts into a new surface mesh

SourceFiles
    triSurfaceCopyParts.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceCopyParts_H
#define triSurfaceCopyParts_H

#include "triSurf.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

/*---------------------------------------------------------------------------*\
                     Class triSurfaceCopyParts Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceCopyParts
{
    // Private data

        //- reference to triSurf
        const triSurf& surf_;


    // Private member functions

        //- mark facets which shall be copied
        void markFacetsForCopying(const wordList&, boolList&) const;

        //- copies data to a new surface mesh
        void copySurfaceMesh(const boolList&, triSurf&) const;

        //- Disallow default bitwise copy construct
        triSurfaceCopyParts(const triSurfaceCopyParts&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceCopyParts&);


public:

    //- Construct from octree
    triSurfaceCopyParts(const triSurf& surface);

    //- Destructor
    ~triSurfaceCopyParts() = default;


    // Member Functions

        //- copies selected patches/subsets to an already created mesh
        void copySurface(const wordList&, triSurf&) const;

        //- creates a new surface mesh and copies selected patches/subsets
        triSurf* copySurface(const wordList&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
