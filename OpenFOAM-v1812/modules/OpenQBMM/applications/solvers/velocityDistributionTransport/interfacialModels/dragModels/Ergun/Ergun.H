/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
-------------------------------------------------------------------------------
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
    Foam::dragModels::Ergun

Description
    Ergun drag law.

    Reference
    \verbatim
        "Eulerian Two-Phase Flow Theory Applied to Fluidization"
        Enwald, H., Peirano, E., Almstedt A-E.
        International Journal of Multiphase Flow, 
        Volume 22, Pages 21-66, 1996
    \endverbatim

SourceFiles
    Ergun.C

\*---------------------------------------------------------------------------*/

#ifndef Ergun_H
#define Ergun_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class Ergun Declaration
\*---------------------------------------------------------------------------*/

class Ergun
:
    public dragModel
{
public:

    //- Runtime type information
    TypeName("Ergun");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Ergun
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Ergun();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe(const label, const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
