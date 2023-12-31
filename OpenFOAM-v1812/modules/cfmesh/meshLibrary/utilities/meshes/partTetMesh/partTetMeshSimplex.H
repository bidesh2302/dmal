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
    Foam::Module::partTetMeshSimplex

Description
    A simplex which is used for smoothing purposes

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef partTetMeshSimplex_H
#define partTetMeshSimplex_H

#include "partTetMesh.H"
#include "parPartTet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

// Forward declarations
class partTetMesh;

/*---------------------------------------------------------------------------*\
                     Class partTetMeshSimplex Declaration
\*---------------------------------------------------------------------------*/

class partTetMeshSimplex
{
    // Private data

        //- points making the simplex
        DynList<point, 128> pts_;

        //- tets making the simplex
        DynList<partTet, 128> tets_;


public:

    // Constructors

        //- Construct from partTetMeshSimplex and point label
        partTetMeshSimplex(const partTetMesh& tm, const label pI);

        //- Construct from the list of parPartTet and point label
        partTetMeshSimplex
        (
            const DynList<parPartTet>& pt,
            const label gpI
        );

        //- contruct from a list of points and a list of partTets
        //- ad a point label
        partTetMeshSimplex
        (
            const DynList<point, 128>& pts,
            const DynList<partTet, 128>& tets,
            const label pointI
        );


    //- Destructor
    ~partTetMeshSimplex() = default;


    // Member functions

        //- return points
        inline DynList<point, 128>& pts()
        {
            return pts_;
        }

        //- return points
        inline const DynList<point, 128>& pts() const
        {
            return pts_;
        }

        //- return tets
        inline const DynList<partTet, 128>& tets() const
        {
            return tets_;
        }

        //- return centre point coordinates
        inline const point& centrePoint() const
        {
            return pts_[tets_[0][3]];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
