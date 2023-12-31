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
    Foam::Module::sphereRefinement

Description
    Mark objects contained in a sphere

SourceFiles
    sphereRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef sphereRefinement_H
#define sphereRefinement_H

#include "objectRefinement.H"
#include "point.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

/*---------------------------------------------------------------------------*\
                      Class sphereRefinement Declaration
\*---------------------------------------------------------------------------*/

class sphereRefinement
:
    public objectRefinement
{
    // Private data

        //- centre of the sphere
        point centre_;

        //- radius of the sphere
        scalar radius_;


public:

    //- Runtime type information
    TypeName("sphere");


    // Constructors

        //- Null construct
        sphereRefinement();

        //- Construct from name, cell size, centre and radius
        sphereRefinement
        (
            const word& name,
            const scalar cellSize,
            const direction additionalRefLevels,
            const point& centre,
            const scalar radius
        );

        //- Construct from dictionary
        sphereRefinement(const word& name, const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<objectRefinement> clone
        (
            const sphereRefinement& sr
        ) const
        {
            return autoPtr<objectRefinement>
            (
                new sphereRefinement
                (
                    sr.name(),
                    sr.cellSize(),
                    sr.additionalRefinementLevels(),
                    sr.centre_,
                    sr.radius_
                )
            );
        }


    // Member Functions

        //- check if a boundBox intersects or is inside the object
        bool intersectsObject(const boundBox&) const;

        //- Return as dictionary of entries
        dictionary dict(bool ignoreType = false) const;


        // Write

            //- Write
            void write(Ostream&) const;

            //- Write dictionary
            void writeDict(Ostream&, bool subDict = true) const;


    // Member Operators

        //- assign from dictionary
        void operator=(const dictionary&);


    // IOstream Operators

        Ostream& operator<<(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
