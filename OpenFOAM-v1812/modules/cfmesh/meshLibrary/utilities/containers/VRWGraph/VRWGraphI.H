/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | www.cfmesh.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 Creative Fields, Ltd.
    Copyright (C) 2021 OpenCFD Ltd.
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

inline void Foam::Module::VRWGraph::checkIndex
(
    const label i,
    const label j
) const
{
    if (i < 0 || i >= rows_.size())
    {
        FatalErrorInFunction
            << "Row index " << i
            << " is not in range [0.." << rows_.size() << ")"
            << abort(FatalError);
    }

    if (j < 0 || j >= rows_[i].size())
    {
        FatalErrorInFunction
            << "Column index " << j
            << " is not in range [0.." << rows_[i].size() << ")"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::Module::VRWGraph::VRWGraph()
:
    data_(),
    rows_()
{}


inline Foam::Module::VRWGraph::VRWGraph
(
    const label size
)
:
    data_(),
    rows_(size)
{
    for (label rowI = 0; rowI < size; ++rowI)
    {
        rows_[rowI].start() = INVALIDROW;
        rows_[rowI].size() = NONE;
    }
}


inline Foam::Module::VRWGraph::VRWGraph
(
    const label nRows,
    const label nColumnsInRow
)
:
    data_(nRows*nColumnsInRow),
    rows_(nRows)
{
    for (label rowI = 0; rowI < nRows; ++rowI)
    {
        rows_[rowI].start() = rowI*nColumnsInRow;
        rows_[rowI].size() = nColumnsInRow;
    }
}


inline Foam::Module::VRWGraph::VRWGraph
(
    const label nRows,
    const label nColumnsInRow,
    const label t
)
:
    data_(nRows*nColumnsInRow, t),
    rows_(nRows)
{
    for (label rowI = 0; rowI < nRows; ++rowI)
    {
        rows_[rowI].start() = rowI*nColumnsInRow;
        rows_[rowI].size() = nColumnsInRow;
    }
}


inline Foam::Module::VRWGraph::VRWGraph
(
    const VRWGraph& ol
)
:
    data_(ol.data_),
    rows_(ol.rows_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::Module::VRWGraph::size() const
{
    return rows_.size();
}


inline Foam::label Foam::Module::VRWGraph::sizeOfRow(const label rowI) const
{
    return rows_[rowI].size();
}


inline void Foam::Module::VRWGraph::setSize(const label size)
{
    if (size > rows_.size())
    {
        rowElement rowInfo(INVALIDROW, NONE);

        for (label i = rows_.size(); i < size; ++i)
        {
            rows_.append(rowInfo);
        }
    }
    else
    {
        rows_.setSize(size);
    }
}


void Foam::Module::VRWGraph::setSizeAndColumnWidth
(
    const label newNumRows,
    const label rcWidth
)
{
    if (rows_.size() != 0)
    {
        FatalErrorInFunction
            << "This function should be used for empty graphs, only!"
            << exit(FatalError);
    }

    data_.setSize(newNumRows*rcWidth);
    data_ = FREEENTRY;

    rows_.setSize(newNumRows);
    label start(0);

    for (label i = 0; i < newNumRows; ++i)
    {
        rows_[i].start() = start;
        rows_[i].size() = 0;
        data_[start] = FREESTART;

        start += rcWidth;
    }
}


template<class ListType>
inline void Foam::Module::VRWGraph::setSizeAndRowSize(const ListType& l)
{
    //- set the size of graph rows
    const label nRows = l.size();
    rows_.setSize(nRows);

    label start(0);
    for (label rowI = 0; rowI < nRows; ++rowI)
    {
        rows_[rowI].size() = l[rowI];

        if (rows_[rowI].size() != NONE)
        {
            rows_[rowI].start() = start;
        }
        else
        {
            rows_[rowI].start() = INVALIDROW;
        }

        start += rows_[rowI].size();
    }

    data_.setSize(start);
}


inline void Foam::Module::VRWGraph::setRowSize
(
    const label rowI,
    const label newSize
)
{
    # ifdef FULLDEBUG
    if ((rowI < 0) || (rowI >= rows_.size()))
    {
        FatalErrorInFunction
            << "Row index " << Foam::label(rowI)
            << " is not in range " << Foam::label(0)
            << " and " << rows_.size() << abort(FatalError);
    }
    # endif

    const label start = rows_[rowI].start();
    if (start == INVALIDROW)
    {
        if (newSize > 0)
        {
            rows_[rowI].start() = data_.size();
            for (label i = 0; i < newSize; ++i)
            {
                data_.append(NONE);
            }
            rows_[rowI].size() = newSize;
        }
    }
    else if (newSize > rows_[rowI].size())
    {
        //- check if there is some unused space after the last element
        bool foundUnused(true);

        for (label i = rows_[rowI].size(); i < newSize; ++i)
        {
            const label j = start + i;
            if
            (
                (j >= data_.size()) ||
                (data_[j] != FREEENTRY) ||
                (data_[j] == FREESTART)
            )
            {
                foundUnused = false;
                break;
            }
        }

        if (foundUnused)
        {
            //- row can be extended without copying
            for (label i = rows_[rowI].size(); i < newSize; ++i)
            {
                data_[start + i] = NONE;
            }
        }
        else
        {
            //- row is copied at the end of the data list
            rows_[rowI].start() = data_.size();
            for (label i = 0; i < rows_[rowI].size(); ++i)
            {
                data_.append(data_[start + i]);
                data_[start + i] = FREEENTRY;
            }
            for (label i = rows_[rowI].size(); i < newSize; ++i)
            {
                data_.append(NONE);
            }
        }

        rows_[rowI].size() = newSize;
    }
    else if (newSize < rows_[rowI].size())
    {
        for (label i = newSize; i < rows_[rowI].size(); ++i)
        {
            data_[start + i] = FREEENTRY;
        }
        rows_[rowI].size() = newSize;
        if (newSize == 0)
        {
            rows_[rowI].start() = INVALIDROW;
        }
    }
}


inline void Foam::Module::VRWGraph::clear()
{
    data_.setSize(0);
    rows_.setSize(0);
}


template<class ListType>
inline void Foam::Module::VRWGraph::appendList
(
    const ListType& l
)
{
    if (l.size() == 0)
    {
        rows_.append(rowElement(INVALIDROW, 0));
        return;
    }

    rowElement rowInfo(data_.size(), l.size());
    const label size = l.size();
    for (label elI = 0; elI < size; ++elI)
    {
        data_.append(l[elI]);
    }
    rows_.append(rowInfo);
}


inline void Foam::Module::VRWGraph::append(const label rowI, const label el)
{
    rowElement& re = rows_[rowI];

    if (re.start() == INVALIDROW)
    {
        re.start() = data_.size();
        re.size() = 1;
        data_.append(el);
    }
    else
    {
        const label oldStart = re.start();
        const label oldSize = re.size();
         ++re.size();

        if (oldStart + oldSize < data_.size())
        {
            if
            (
                (data_[oldStart + oldSize] == FREEENTRY) ||
                (data_[oldStart + oldSize] == FREESTART)
            )
            {
                data_[oldStart + oldSize] = el;
            }
            else
            {
                re.start() = data_.size();
                for (label i = 0; i < oldSize; ++i)
                {
                    data_.append(data_[oldStart + i]);
                    data_[oldStart + i] = FREEENTRY;
                }
                data_.append(el);
            }
        }
        else
        {
            data_.append(el);
        }
    }
}


inline void Foam::Module::VRWGraph::appendUniq
(
    const label rowI,
    const label el
)
{
    if (!found(rowI, el))
    {
        append(rowI, el);
    }
}


template<class ListType>
inline void Foam::Module::VRWGraph::setRow
(
    const label rowI,
    const ListType& l
)
{
    this->setRowSize(rowI, l.size());
    const label start = rows_[rowI].start();
    const label size = l.size();
    for (label elI = 0; elI < size; ++elI)
    {
        data_[start + elI] = l[elI];
    }
}


inline void Foam::Module::VRWGraph::mergeGraphs
(
    const List<VRWGraph>& graphParts
)
{
    const label nGraphs = graphParts.size();
    const label nRows = graphParts[0].size();
    forAll(graphParts, i)
    {
        if (nRows != graphParts[i].size())
        {
            FatalErrorInFunction
                << "Cannot merge graphs" << abort(FatalError);
        }
    }

    //- find the number of elements in each row
    labelLongList nElmtsInRow(nRows);
    for (label rowI = 0; rowI < nRows; ++rowI)
    {
        label sum(0);
        for (label i = 0; i < nGraphs; ++i)
        {
            sum += graphParts[i].sizeOfRow(rowI);
        }

        nElmtsInRow[rowI] = sum;
    }

    setSizeAndRowSize(nElmtsInRow);

    //- Finally, assemble the merged graph
    for (label rowI = 0; rowI < nRows; ++rowI)
    {
        forAll(graphParts, i)
        {
            const VRWGraph& gp = graphParts[i];
            for (label j = 0; j < gp.sizeOfRow(rowI); ++j)
            {
                this->operator()(rowI, --nElmtsInRow[rowI]) = gp(rowI, j);
            }
        }
    }
}


template<class GraphType>
inline void Foam::Module::VRWGraph::reverseAddressing
(
    const label nRows,
    const GraphType& origGraph
)
{
    const label origSize = origGraph.size();
    labelLongList nElmtsInRow(nRows);

    for (label rowI = 0; rowI < nRows; ++rowI)
    {
        nElmtsInRow[rowI] = 0;
    }

    for (label rowI = 0; rowI < origSize; ++rowI)
    {
        const label rowSize = origGraph[rowI].size();

        for (label i = 0; i < rowSize; ++i)
        {
            ++nElmtsInRow[origGraph[rowI][i]];
        }
    }

    setSizeAndRowSize(nElmtsInRow);
    nElmtsInRow = 0;

    //- finally fill in the data
    for (label rowI = 0; rowI < origSize; ++rowI)
    {
        const label rowSize = origGraph[rowI].size();

        for (label i = 0; i < rowSize; ++i)
        {
            const label el = origGraph[rowI][i];
            this->operator()(el, nElmtsInRow[el]++) = rowI;
        }
    }
}


template<class GraphType>
inline void Foam::Module::VRWGraph::reverseAddressing
(
    const GraphType& origGraph
)
{
    const label size = origGraph.size();
    label maxValue(-1);

    for (label rowI = 0; rowI < size; ++rowI)
    {
        const label rowSize = origGraph[rowI].size();
        for (label i = 0; i < rowSize; ++i)
        {
            maxValue = Foam::max(maxValue, origGraph[rowI][i]);
        }
    }

    ++maxValue;
    reverseAddressing(maxValue, origGraph);
}


inline void Foam::Module::VRWGraph::reverseAddressing
(
    const label nRows,
    const VRWGraph& origGraph
)
{
    const label origSize = origGraph.size();
    labelLongList nElmtsInRow(nRows);

    for (label rowI = 0; rowI < nRows; ++rowI)
    {
        nElmtsInRow[rowI] = 0;
    }

    for (label rowI = 0; rowI < origSize; ++rowI)
    {
        const label rowSize = origGraph.sizeOfRow(rowI);

        for (label i = 0; i < rowSize; ++i)
        {
            ++nElmtsInRow[origGraph(rowI, i)];
        }
    }

    setSizeAndRowSize(nElmtsInRow);
    nElmtsInRow = 0;

    //- finally fill in the data
    for (label rowI = 0; rowI < origSize; ++rowI)
    {
        const label rowSize = origGraph.sizeOfRow(rowI);

        for (label i = 0; i < rowSize; ++i)
        {
            const label el = origGraph(rowI, i);
            this->operator()(el, nElmtsInRow[el]++) = rowI;
        }
    }
}


inline void Foam::Module::VRWGraph::reverseAddressing
(
    const VRWGraph& origGraph
)
{
    const label size = origGraph.size();
    label maxValue(-1);

    for (label rowI = 0; rowI < size; ++rowI)
    {
        const label rowSize = origGraph.sizeOfRow(rowI);
        for (label i = 0; i < rowSize; ++i)
        {
            maxValue = Foam::max(maxValue, origGraph(rowI, i));
        }
    }

    ++maxValue;
    reverseAddressing(maxValue, origGraph);
}


inline bool Foam::Module::VRWGraph::found
(
    const label rowI,
    const label e
) const
{
    return (this->find(rowI, e) >= 0);
}


inline Foam::label Foam::Module::VRWGraph::find
(
    const label rowI,
    const label e
) const
{
    const label start = rows_[rowI].start();
    if (start == INVALIDROW)
    {
        return -1;
    }

    const label size = rows_[rowI].size();

    for (label i = 0; i < size; ++i)
    {
        if (data_[start + i] == e)
        {
            return i;
        }
    }

    return -1;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::label Foam::Module::VRWGraph::operator()
(
    const label i,
    const label j
) const
{
    #ifdef FULLDEBUG
    checkIndex(i, j);
    #endif

    return data_[rows_[i].start() + j];
}


inline Foam::label& Foam::Module::VRWGraph::operator()
(
    const label i, const label j
)
{
    #ifdef FULLDEBUG
    checkIndex(i, j);
    #endif

    return data_[rows_[i].start() + j];
}


inline Foam::Module::constRow
Foam::Module::VRWGraph::operator[](const label i) const
{
    return constRow(*this, i);
}


inline Foam::Module::row Foam::Module::VRWGraph::operator[](const label i)
{
    return row(*this, i);
}


inline void Foam::Module::VRWGraph::operator=
(
    const VRWGraph& l
)
{
    data_ = l.data_;
    rows_ = l.rows_;
}


// ************************************************************************* //
