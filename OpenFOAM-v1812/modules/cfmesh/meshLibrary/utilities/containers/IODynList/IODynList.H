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
    Foam::Module::IODynList

Description
    An IODynList of a given type is a DynList of that type which supports
    automated input and output.

SourceFiles
    IODynList.C

\*---------------------------------------------------------------------------*/

#ifndef IODynList_H
#define IODynList_H

#include "DynList.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

/*---------------------------------------------------------------------------*\
                          Class IODynList Declaration
\*---------------------------------------------------------------------------*/

template<class T, class IndexType = label>
class IODynList
:
    public regIOobject,
    public DynList<T, IndexType>
{
public:

    //- Runtime type information
    TypeName("DynList");


    // Constructors

        //- Construct from IOobject
        IODynList(const IOobject&);

        //- Construct from IOobject and size of IODynList
        IODynList(const IOobject&, const IndexType);

        //- Construct from IOobject and a List
        IODynList(const IOobject&, const DynList<T, IndexType>&);


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const IODynList<T, IndexType>&);

        void operator=(const DynList<T, IndexType>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IODynList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
