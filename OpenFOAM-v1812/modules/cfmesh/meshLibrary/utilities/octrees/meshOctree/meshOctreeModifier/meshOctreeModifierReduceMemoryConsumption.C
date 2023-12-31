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

#include "meshOctreeModifier.H"
#include "triSurf.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeModifier::reduceMemoryConsumption()
{
    //createListOfLeaves();

    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    forAll(octree_.dataSlots_, slotI)
    {
        // deleting triangles
        VRWGraph& containedTriangles =
            octree_.dataSlots_[slotI].containedTriangles_;
        label nElmts = containedTriangles.size();
        boolList mayDeleteData(nElmts, true);
        forAll(leaves, leafI)
        {
            const meshOctreeCube& oc = *leaves[leafI];

            if (
                oc.hasContainedElements() &&
                oc.slotPtr() == &octree_.dataSlots_[slotI]
            )
                mayDeleteData[oc.containedElements()] = false;
        }

        for (label i = 0; i < nElmts; ++i)
            if (mayDeleteData[i])
                containedTriangles.setRowSize(i, 0);
        containedTriangles.optimizeMemoryUsage();

        // deleting edges
        VRWGraph& containedEdges =
            octree_.dataSlots_[slotI].containedEdges_;
        nElmts = containedEdges.size();
        mayDeleteData.setSize(nElmts);
        mayDeleteData = true;
        forAll(leaves, leafI)
        {
            const meshOctreeCube& oc = *leaves[leafI];

            if
            (
                oc.hasContainedEdges()
             && oc.slotPtr() == &octree_.dataSlots_[slotI]
            )
                mayDeleteData[oc.containedEdges()] = false;
        }

        for (label i = 0; i < nElmts; ++i)
            if (mayDeleteData[i])
                containedEdges.setRowSize(i, 0);

        containedEdges.optimizeMemoryUsage();
    }
}


// ************************************************************************* //
