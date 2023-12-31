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
    Foam::Module::symmetryPlaneOptimisation

Description
    Smoothing of symmetry planes in the mesh such that all points
    are in the plane.

SourceFiles
    symmetryPlaneOptimisation.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryPlaneOptimisation_H
#define symmetryPlaneOptimisation_H

#include "DynList.H"
#include "polyMeshGenModifier.H"
#include "boundBox.H"
#include "labelLongList.H"
#include "boolList.H"
#include "plane.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

/*---------------------------------------------------------------------------*\
                  Class symmetryPlaneOptimisation Declaration
\*---------------------------------------------------------------------------*/

class symmetryPlaneOptimisation
{
    // Private data

        //- reference to the mesh
        polyMeshGen& mesh_;

        //- symmetry planes in the mesh
        std::map<label, plane> symmetryPlanes_;


    // Private member functions

        //- detect symmetry planes
        void detectSymmetryPlanes();

        //- point-planes addressing
        bool pointInPlanes(VRWGraph&) const;


public:

    //- Construct from mesh
    symmetryPlaneOptimisation(polyMeshGen& mesh);


    //- Destructor
    ~symmetryPlaneOptimisation() = default;


    // Member Functions

        //- move vertices to the symmetry planes
        void optimizeSymmetryPlanes();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
