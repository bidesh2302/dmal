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
    Foam::Module::faceIOGraph

Description
    A graph of faces which support automated output

SourceFiles
    faceIOGraph.C

\*---------------------------------------------------------------------------*/

#ifndef faceIOGraph_H
#define faceIOGraph_H

#include "VRWGraph.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

/*---------------------------------------------------------------------------*\
                         Class faceIOGraph Declaration
\*---------------------------------------------------------------------------*/

class faceIOGraph
:
    public regIOobject,
    public VRWGraph
{
public:

    //- Runtime type information
    TypeName("faceIOGraph");

    // Constructors

        //- Construct from IOobject
        faceIOGraph(const IOobject&);

        //- Construct from IOobject and size of faceIOGraph
        faceIOGraph(const IOobject&, const label);

        //- Construct from IOobject and a VRWGraph
        faceIOGraph(const IOobject&, const VRWGraph&);


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const faceIOGraph&);
        void operator=(const VRWGraph&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
