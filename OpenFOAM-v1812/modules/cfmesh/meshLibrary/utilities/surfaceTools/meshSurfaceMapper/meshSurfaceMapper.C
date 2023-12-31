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

#include "meshSurfaceMapper.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "triSurf.H"
#include "triSurfacePartitioner.H"
#include "demandDrivenData.H"
#include "meshOctree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceMapper::createMeshSurfacePartitioner() const
{
    surfaceEnginePartitionerPtr_ = new meshSurfacePartitioner(surfaceEngine_);
}


void Foam::Module::meshSurfaceMapper::createTriSurfacePartitioner() const
{
    surfPartitionerPtr_ = new triSurfacePartitioner(meshOctree_.surface());
}


void Foam::Module::meshSurfaceMapper::clearOut()
{
    if (deletePartitioner_)
        deleteDemandDrivenData(surfaceEnginePartitionerPtr_);
    deleteDemandDrivenData(surfPartitionerPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::meshSurfaceMapper::meshSurfaceMapper
(
    const meshSurfaceEngine& mse,
    const meshOctree& octree
)
:
    surfaceEngine_(mse),
    meshOctree_(octree),
    surfaceEnginePartitionerPtr_(nullptr),
    deletePartitioner_(true),
    surfPartitionerPtr_(nullptr)
{
    if (Pstream::parRun())
    {
        // allocate bpAtProcs and other addressing
        // this is done here to prevent possible deadlocks
        surfaceEngine_.bpAtProcs();
    }
}


Foam::Module::meshSurfaceMapper::meshSurfaceMapper
(
    const meshSurfacePartitioner& mPart,
    const meshOctree& octree
)
:
    surfaceEngine_(mPart.surfaceEngine()),
    meshOctree_(octree),
    surfaceEnginePartitionerPtr_(&mPart),
    deletePartitioner_(false),
    surfPartitionerPtr_(nullptr)
{
    if (Pstream::parRun())
    {
        // allocate bpAtProcs and other addressing
        // this is done here to prevent possible deadlocks
        surfaceEngine_.bpAtProcs();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Module::meshSurfaceMapper::~meshSurfaceMapper()
{
    clearOut();
}


// ************************************************************************* //
