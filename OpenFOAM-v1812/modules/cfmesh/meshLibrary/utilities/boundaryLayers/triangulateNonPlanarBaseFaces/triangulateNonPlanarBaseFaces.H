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
    Foam::Module::triangulateNonPlanarBaseFaces

Description
    Splits selected boundary layer cells into triangular prisms

SourceFiles
    triangulateNonPlanarBaseFaces.C
    triangulateNonPlanarBaseFacesFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef triangulateNonPlanarBaseFaces_H
#define triangulateNonPlanarBaseFaces_H

#include "polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class dictionary;

namespace Module
{

/*---------------------------------------------------------------------------*\
                Class triangulateNonPlanarBaseFaces Declaration
\*---------------------------------------------------------------------------*/

class triangulateNonPlanarBaseFaces
{
    // Private data

        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- boolean list containing information about inverted cells
        boolList invertedCell_;

        //- classification of faces in the mesh
        boolList decomposeFace_;

        //- relative deviation compared to the layer thickness
        scalar tol_;


    // Private member functions

        //- find faces with non-planarity greater than the required
        //- layer thickness
        bool findNonPlanarBoundaryFaces();

        //- decompose existing faces
        void decomposeBoundaryFaces();

        //- decompose adjacent cells into pyramids
        void decomposeCellsIntoPyramids();

        //- Disallow bitwise copy construct
        triangulateNonPlanarBaseFaces
        (
            const triangulateNonPlanarBaseFaces&
        );

        //- Disallow bitwise assignment
        void operator=(const triangulateNonPlanarBaseFaces&);


public:

    //- Construct from mesh
    triangulateNonPlanarBaseFaces(polyMeshGen& mesh);

    //- Destructor
    ~triangulateNonPlanarBaseFaces() = default;


    // Public member functions

        //- set the relative tolerance between the requested boundary
        //- layer thickness and the deviation from planarity
        //- the default is 1.0
        void setRelativeTolerance(const scalar tol);

        //- splits inverted boundary layer prisms into triangular prisms
        void triangulateLayers();


    // Static member functions

        static void readSettings
        (
            const dictionary&,
            triangulateNonPlanarBaseFaces&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
