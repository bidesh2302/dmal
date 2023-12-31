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

#include "meshSurfaceEdgeExtractor.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::meshSurfaceEdgeExtractor::meshSurfaceEdgeExtractor
(
    polyMeshGen& mesh,
    const meshOctree& octree,
    const labelList& pointRegion
)
:
    mesh_(mesh),
    nPoints_(mesh.points().size()),
    boundaryCell_(mesh.cells().size(), false),
    nFacesInCell_(mesh.cells().size(), direction(0)),
    meshOctree_(octree),
    pointRegions_(pointRegion.size())
{
    forAll(pointRegion, pointI)
        if (pointRegion[pointI] != -1)
            pointRegions_.append(pointI, pointRegion[pointI]);

    createEdgeVertices();

    removeOldBoundaryFaces();

    createBoundaryFaces();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceEdgeExtractor::removeOldBoundaryFaces()
{
    const labelList neighbour_ = mesh_.neighbour();
    polyMeshGenModifier meshModifier_(mesh_);
    cellListPMG& cells_ = meshModifier_.cellsAccess();

    forAll(cells_, cellI)
    {
        const cell& c = cells_[cellI];

        cell newC(c);

        forAll(c, fI)
            if (neighbour_[c[fI]] != -1)
            {
                boundaryCell_[cellI] = true;
                newC[nFacesInCell_[cellI]++] = c[fI];
            }

        if (nFacesInCell_[cellI] < direction(c.size()))
        {
            newC.setSize(nFacesInCell_[cellI]);

            cells_[cellI] = newC;
        };
    }

    PtrList<boundaryPatch>& boundaries = meshModifier_.boundariesAccess();
    boundaries.setSize(1);
    boundaries[0].patchSize() = 0;
    meshModifier_.facesAccess().setSize(boundaries[0].patchStart());
}


// ************************************************************************* //
