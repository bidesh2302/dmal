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

Description
    cfMesh utility to merge the supplied list of patches onto a single
    patch.

Author
    Ivor Clifford<ivor.clifford@psi.ch>

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "autoPtr.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "demandDrivenData.H"
#include "Pair.H"

using namespace Foam;
using namespace Foam::Module;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Find the supplied list of patch names and return a list of patch Ids
void getPatchIds
(
    const triSurf& origSurf,
    const wordList& patchNames,
    DynamicList<label>& patchIds
)
{
    const geometricSurfacePatchList& origPatches = origSurf.patches();

    // Create patch name map
    HashSet<word> patchNameHash(patchNames);

    // Find selected patches
    label nFound = 0;
    forAll(origPatches, patchI)
    {
        if (patchNameHash.found(origPatches[patchI].name()))
        {
            patchIds.append(patchI);
            nFound++;
        }
    }

    if (nFound != patchNames.size())
    {
        WarningInFunction
            << "Not all supplied patch names were found on the surface mesh"
            << endl;
    }
}


// Copy all face subsets from one triSurf to another
void copyFaceSubsets
(
    const triSurf& origSurf,
    triSurf& newSurf
)
{
    DynList<label> subsetIds;
    origSurf.facetSubsetIndices(subsetIds);

    forAll(subsetIds, subsetI)
    {
        label newSubsetId = newSurf.addFacetSubset
        (
            origSurf.facetSubsetName(subsetI)
        );

        labelList origFaces;
        origSurf.facetsInSubset(subsetI, origFaces);

        forAll(origFaces, faceI)
        {
            newSurf.addFacetToSubset
            (
                newSubsetId,
                origFaces[faceI]
            );
        }
    }
}


// Copy all edge subsets from one triSurf to another
void copyEdgeSubsets
(
    const triSurf& origSurf,
    triSurf& newSurf
)
{
    DynList<label> subsetIds;
    origSurf.edgeSubsetIndices(subsetIds);

    forAll(subsetIds, subsetI)
    {
        label newSubsetId = newSurf.addEdgeSubset
        (
            origSurf.edgeSubsetName(subsetI)
        );

        labelList origEdges;
        origSurf.edgesInSubset(subsetI, origEdges);

        forAll(origEdges, faceI)
        {
            newSurf.addEdgeToSubset
            (
                newSubsetId,
                origEdges[faceI]
            );
        }
    }
}


// Copy all point subsets from one triSurf to another
void copyPointSubsets
(
    const triSurf& origSurf,
    triSurf& newSurf
)
{
    DynList<label> subsetIds;
    origSurf.pointSubsetIndices(subsetIds);

    forAll(subsetIds, subsetI)
    {
        label newSubsetId = newSurf.addPointSubset
        (
            origSurf.pointSubsetName(subsetI)
        );

        labelList origPoints;
        origSurf.pointsInSubset(subsetI, origPoints);

        forAll(origPoints, faceI)
        {
            newSurf.addPointToSubset
            (
                newSubsetId,
                origPoints[faceI]
            );
        }
    }
}


// Merge the supplied list of patchIds onto a new patch
autoPtr<triSurf> mergeSurfacePatches
(
    const triSurf& origSurf, // Surface
    const UList<label>& patchIds, // Ids of patches to merge
    const word& newPatchName,       // Name of new (merged) patch // Name of new (merged) patch
    bool keepPatches // Keep the original patches - they will be emptied
)
{
    const geometricSurfacePatchList& origPatches = origSurf.patches();
    const LongList<labelledTri>& origFacets = origSurf.facets();

    label newPatchId = origPatches.size();

    // Determine new patch type
    word newPatchType = origPatches[patchIds[0]].geometricType();

    // Create patch addressing
    List<DynamicList<label>> patchAddr(origPatches.size()+1);

    forAll(origFacets, faceI)
    {
        patchAddr[origFacets[faceI].region()].append(faceI);
    }

    // Move selected patches to new patch
    forAll(patchIds, patchI)
    {
        patchAddr[newPatchId].append(patchAddr[patchIds[patchI]]);
        patchAddr[patchIds[patchI]].clear();
    }

    // Create new facets list
    LongList<labelledTri> newFacets(origFacets.size());
    labelList newFaceAddr(origFacets.size(), -1);

    label patchCount = 0;
    label faceI = 0;
    forAll(patchAddr, patchI)
    {
        const UList<label>& addr = patchAddr[patchI];

        if (addr.size())
        {
            forAll(addr, i)
            {
                newFacets[faceI] = origFacets[addr[i]];
                newFacets[faceI].region() = patchCount;

                newFaceAddr[addr[i]] = faceI;

                faceI++;
            }
        }

        if (addr.size() || keepPatches)
        {
            patchCount++;
        }
    }

    // Create new patch list
    geometricSurfacePatchList newPatches(patchCount);

    patchCount = 0;
    forAll(origPatches, patchI)
    {
        // Only add patches if they contain faces
        if (patchAddr[patchI].size())
        {
            newPatches[patchCount] = origPatches[patchI];
            newPatches[patchCount].index() = patchCount;
        }

        if (patchAddr[patchI].size() || keepPatches)
        {
            patchCount++;
        }
    }

    // Add new patch if it contains faces
    if (patchAddr[patchAddr.size()-1].size())
    {
        newPatches[patchCount] = geometricSurfacePatch
        (
            newPatchName,
            patchCount,
            newPatchType
        );
    }
    if (patchAddr[patchAddr.size()-1].size() || keepPatches)
    {
        patchCount++;
    }

    // Create new surface
    autoPtr<triSurf> newSurf
    (
        new triSurf
        (
            newFacets,
            newPatches,
            origSurf.featureEdges(),
            origSurf.points()
        )
    );

    // Transfer face subsets
    copyFaceSubsets(origSurf, newSurf());
    newSurf->updateFacetsSubsets(newFaceAddr);

    // Transfer feature edge subsets
    copyEdgeSubsets(origSurf, newSurf());

    // Transfer point subsets
    copyPointSubsets(origSurf, newSurf());

    // Done
    return newSurf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "(cfmesh)\n"
        "Merge the supplied list of patches onto a single patch."
    );

    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
    argList::validArgs.append("new patch");
    argList::addOption("patchNames", "list of names");
    argList::addOption("patchIds", "list of patchIds");
    argList::addOption("patchIdRange", "(start end)");
    argList::addOption("output", "file name (default overwrite)");
    argList::addBoolOption("keep");
    argList args(argc, argv);

    // Process commandline arguments
    fileName inFileName(args[1]);

    word newPatchName(args[2]);

    fileName outFileName(inFileName);

    args.readIfPresent("output", outFileName);

    const bool keepPatches = args.found("keep");

    // Read original surface
    triSurf origSurf(inFileName);

    // Get patch ids
    DynamicList<label> patchIds;

    if (args.found("patchNames"))
    {
        if (args.found("patchIds"))
        {
            FatalErrorInFunction
                << "Cannot specify both patch names and ids" << nl
                << Foam::abort(FatalError);
        }

        wordList patchNames(args.getList<word>("patchNames"));

        getPatchIds
        (
            origSurf,
            patchNames,
            patchIds
        );
    }
    else if (args.found("patchIds"))
    {
        patchIds.append(args.getList<label>("patchIds"));
    }

    labelPair idRange;
    if (args.readIfPresent("patchIdRange", idRange))
    {
        for (label id = idRange.first(); id <= idRange.second(); ++id)
        {
            patchIds.append(id);
        }
    }

    if (!patchIds.size())
    {
        FatalErrorInFunction
            << "No patches specified" << nl
            << Foam::abort(FatalError);
    }

    // Merge patches
    autoPtr<triSurf> newSurf = mergeSurfacePatches
    (
        origSurf,
        patchIds,
        newPatchName,
        keepPatches
    );

    // Write new surface mesh
    newSurf->writeSurface(outFileName);

    Info<< "Original surface patches: " << origSurf.patches().size() << nl
        << "Final surface patches: " << newSurf->patches().size() << nl
        << "Surface written to " << outFileName <<  endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
