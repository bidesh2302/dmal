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

\*---------------------------------------------------------------------------*/

#include "tetCreatorOctree.H"
#include "meshOctree.H"
#include "demandDrivenData.H"

//#define DEBUGTets

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Module::meshOctreeCubeCoordinates
Foam::Module::tetCreatorOctree::edgeCoordinates_[12][4]=
    {
        {
            // edge 0
            meshOctreeCubeCoordinates(0, 0, 1, 0),
            meshOctreeCubeCoordinates(0, -1, 1, 0),
            meshOctreeCubeCoordinates(0, -1, 0, 0),
            meshOctreeCubeCoordinates(0, 0, 0, 0)
        },
        {
            // edge 1
            meshOctreeCubeCoordinates(0, 0, 0, 0),
            meshOctreeCubeCoordinates(0, -1, 0, 0),
            meshOctreeCubeCoordinates(0, -1, -1, 0),
            meshOctreeCubeCoordinates(0, 0, -1, 0)
        },
        {
            // edge 2
            meshOctreeCubeCoordinates(0, 1, 0, 0),
            meshOctreeCubeCoordinates(0, 0, 0, 0),
            meshOctreeCubeCoordinates(0, 0, -1, 0),
            meshOctreeCubeCoordinates(0, 1, -1, 0)
        },
        {
            // edge 3
            meshOctreeCubeCoordinates(0, 1, 1, 0),
            meshOctreeCubeCoordinates(0, 0, 1, 0),
            meshOctreeCubeCoordinates(0, 0, 0, 0),
            meshOctreeCubeCoordinates(0, 1, 0, 0)
        },
        {
            // edge 4
            meshOctreeCubeCoordinates(0, 0, 0, 0),
            meshOctreeCubeCoordinates(0, 0, 1, 0),
            meshOctreeCubeCoordinates(1, 0, 1, 0),
            meshOctreeCubeCoordinates(1, 0, 0, 0)
        },
        {
            // edge 5
            meshOctreeCubeCoordinates(0, 0, -1, 0),
            meshOctreeCubeCoordinates(0, 0, 0, 0),
            meshOctreeCubeCoordinates(1, 0, 0, 0),
            meshOctreeCubeCoordinates(1, 0, -1, 0)
        },
        {
            // edge 6
            meshOctreeCubeCoordinates(-1, 0, -1, 0),
            meshOctreeCubeCoordinates(-1, 0, 0, 0),
            meshOctreeCubeCoordinates(0, 0, 0, 0),
            meshOctreeCubeCoordinates(0, 0, -1, 0)
        },
        {
            // edge 7
            meshOctreeCubeCoordinates(-1, 0, 0, 0),
            meshOctreeCubeCoordinates(-1, 0, 1, 0),
            meshOctreeCubeCoordinates(0, 0, 1, 0),
            meshOctreeCubeCoordinates(0, 0, 0, 0)
        },
        {
            // edge 8
            meshOctreeCubeCoordinates(-1, 0, 0, 0),
            meshOctreeCubeCoordinates(-1, -1, 0, 0),
            meshOctreeCubeCoordinates(0, -1, 0, 0),
            meshOctreeCubeCoordinates(0, 0, 0, 0)
        },
        {
            // edge 9
            meshOctreeCubeCoordinates(0, 0, 0, 0),
            meshOctreeCubeCoordinates(0, -1, 0, 0),
            meshOctreeCubeCoordinates(1, -1, 0, 0),
            meshOctreeCubeCoordinates(1, 0, 0, 0)
        },
        {
            // edge 10
            meshOctreeCubeCoordinates(0, 1, 0, 0),
            meshOctreeCubeCoordinates(0, 0, 0, 0),
            meshOctreeCubeCoordinates(1, 0, 0, 0),
            meshOctreeCubeCoordinates(1, 1, 0, 0)
        },
        {
            // edge 11
            meshOctreeCubeCoordinates(-1, 1, 0, 0),
            meshOctreeCubeCoordinates(-1, 0, 0, 0),
            meshOctreeCubeCoordinates(0, 0, 0, 0),
            meshOctreeCubeCoordinates(0, 1, 0, 0)
        }
    };

const Foam::label Foam::Module::tetCreatorOctree::faceCentreHelper_[3][4] =
    {
        {3, 5, 2, 4},
        {5, 1, 4, 0},
        {1, 3, 0, 2}
    };


void Foam::Module::tetCreatorOctree::createTets()
{
    createPointsAndAddressing();

    createTetsFromFacesWithCentreNode();

    createTetsAroundSplitEdges();

    createTetsAroundEdges();

    createTetsFromSplitFaces();

    clearOut();
    sortedLeaves_.setSize(0);

    created_ = true;
}


void Foam::Module::tetCreatorOctree::clearOut()
{
    sortedLeaves_.clear();
    deleteDemandDrivenData(subNodeLabelsPtr_);
    deleteDemandDrivenData(cubeLabelPtr_);
    deleteDemandDrivenData(faceCentreLabelPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::tetCreatorOctree::tetCreatorOctree
(
    const meshOctree& octree,
    const IOdictionary& meshDict
)
:
    octreeCheck_(octree, meshDict, true),
    tetPoints_(),
    tets_(),
    sortedLeaves_(),
    subNodeLabelsPtr_(nullptr),
    cubeLabelPtr_(nullptr),
    faceCentreLabelPtr_(nullptr),
    created_(false)
{
    createTets();

    clearOut();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Module::tetCreatorOctree::~tetCreatorOctree()
{
    clearOut();
}


// ************************************************************************* //
