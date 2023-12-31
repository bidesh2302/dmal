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
    Foam::Module::triSurfaceMetaData

Description
    Generates meta data of the surface mesh.

    It contains:
    - The number of points
    - The number of triangles
    - The number of patches
    - The number of feature edges
    - The number of subsets, names, and the number of elements in each subset

SourceFiles
    triSurfaceMetaData.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMetaData_H
#define triSurfaceMetaData_H

#include "triSurf.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

/*---------------------------------------------------------------------------*\
                     Class triSurfaceMetaData Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMetaData
{
    // Private data

        //- const reference to triSurf
        const triSurf& surf_;

        //- patches/subsets for removal
        dictionary metaDict_;


    // Private member functions

        //- create dictionary with meta data
        void createMetaData();

        //- Disallow default bitwise copy construct
        triSurfaceMetaData(const triSurfaceMetaData&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceMetaData&);


public:

    //- Construct from triSurf
    triSurfaceMetaData(const triSurf& surface);

    //- Destructor
    ~triSurfaceMetaData() = default;


    // Member Functions

        //- return a constant reference to meta data
        const dictionary& metaData() const
        {
            return metaDict_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
