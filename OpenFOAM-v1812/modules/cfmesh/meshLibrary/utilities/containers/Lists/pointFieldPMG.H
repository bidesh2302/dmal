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
    Foam::Module::pointFieldPMG

Description
    This is a container with additional size to prevent re-allocation
    every time it is resized

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef pointFieldPMG_H
#define pointFieldPMG_H

#include "regIOobject.H"
#include "pointField.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
class pointFieldPMG;
inline Ostream& operator<<(Ostream&, const pointFieldPMG&);
inline Istream& operator>>(Istream&, pointFieldPMG&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class pointFieldPMG
:
    public regIOobject,
    public pointField
{
    // Private data

        //- number of used elements
        label nElmts_;


    // Disallow bitwise assignment

        void operator=(const pointFieldPMG&);

        pointFieldPMG(const pointFieldPMG&);


    // Disallow transfer from pointField

        void transfer(pointField&);


public:

    //- Declare type name
    TypeName("vectorField");


    // Constructors

        //- construct from IOobject
        inline pointFieldPMG(const IOobject&);

        //- construct from IOobject and size
        inline pointFieldPMG(const IOobject&, const label);

        //- construct from IOobject and pointField
        inline pointFieldPMG(const IOobject&, const pointField&);


    //- Destructor
    ~pointFieldPMG() = default;


    // Member functions

        //- return the number of used elements
        inline label size() const;

        //- set the number of used elements
        inline void setSize(const label nElmts);

        //- reserve the memory for the specified number of elements
        //- the request is ignored if the specified value is smaller than
        //- the current number of elements
        //- does not change the size
        inline void reserve(const label capacity);

        //- set the size to zero
        inline void clear();

        //- add a point at the end of the list
        inline void append(const point&);

        //- return an element with bound checking
        inline point& newElmt(const label);

        //- read/write the list onto disk
        inline bool writeData(Ostream&) const;


    // Member operators

        inline void operator=(const pointField&);

        inline friend Ostream& operator<<(Ostream&, const pointFieldPMG&);

        inline friend Istream& operator>>(Istream&, pointFieldPMG&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointFieldPMGI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
