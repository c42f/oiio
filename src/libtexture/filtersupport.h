/*
  Copyright 2009 Larry Gritz and the other authors and contributors.
  All Rights Reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:
  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of the software's owners nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  (This is the Modified BSD License)
*/

/// \file
///
/// \brief Utilities for dealing with filter supports, including wrapping the
/// support using some wrap mode.

#ifndef OPENIMAGEIO_FILTERSUPPORT_H
#define OPENIMAGEIO_FILTERSUPPORT_H

#include "fmath.h"

OIIO_NAMESPACE_ENTER
{

//------------------------------------------------------------------------------
/// \brief Represent a filter support region in integer raster coordinates.
///
/// A support is represented by two integers, "start" and "end".  These specify
/// the range [start,end) where the left end is *inclusive* and the right end
/// *exclusive*.
///
/// \todo Clean up this class for a more uniform interface (make things more
/// uniform?)
///
struct FilterSupport1D
{
    /// Start of the support, inclusive
    int start;
    /// End of the support range, exclusive.
    int end;
    /// Trivial constructor
    FilterSupport1D(int start, int end);
    /// \brief Truncate the support into the interval [rangeStart, rangeEnd)
    ///
    /// Note that this may result in empty supports.
    ///
    void truncate(int rangeStart, int rangeEnd);
    /// Return the number of points in the support.
    int range() const;
    /// Return true if the support is empty.
    bool isEmpty() const;
    /// Return true if the support covers part of the given range.
    bool intersectsRange(int rangeStart, int rangeEnd) const;
    /// Return true if the support is wholly inside the given range.
    bool inRange(int rangeStart, int rangeEnd) const;
};

/// Return the intersection of two support regions.
FilterSupport1D intersect(const FilterSupport1D s1, const FilterSupport1D s2);

//------------------------------------------------------------------------------
/// \brief Hold filter support area.
///
/// The end markers are *exclusive*, so the support of a filter is inside the
/// rectangular region [sx.start, ..., endX-1] x [startY, ..., endY-1].
///
struct FilterSupport
{
    FilterSupport1D sx; ///< support in x-direction
    FilterSupport1D sy; ///< support in y-direction
    /// Trivial constructor.
    FilterSupport(int startX = 0, int endX = 0, int startY = 0, int endY = 0);
    /// Constructor from two 1D supports
    FilterSupport(const FilterSupport1D s1, const FilterSupport1D s2);
    /// Return area of the support in number of pixels.
    int area() const;
    /// Return true if the support is an empty set.
    bool isEmpty() const;
    /// Return true if the support covers part of the given range.
    bool intersectsRange(int startX, int endX, int startY, int endY) const;
    /// Return true if the support is wholly inside the given range.
    bool inRange(int startX, int endX, int startY, int endY) const;
};

/// Return the intersection of two supports
FilterSupport intersect(const FilterSupport& s1, const FilterSupport& s2);

//==============================================================================
// Implementation details
//==============================================================================
// FilterSupport1D

inline FilterSupport1D::FilterSupport1D(int start, int end)
    : start(start), end(end)
{ }

inline void FilterSupport1D::truncate(int rangeStart, int rangeEnd)
{
    if(start < rangeStart)
        start = rangeStart;
    if(end > rangeEnd)
        end = rangeEnd;
}

inline int FilterSupport1D::range() const
{
    return end - start;
}

inline bool FilterSupport1D::isEmpty() const
{
    return start >= end;
}

inline bool FilterSupport1D::intersectsRange(
        int rangeStart, int rangeEnd) const
{
    return !isEmpty() && !(rangeStart >= end || rangeEnd <= start);
}

inline bool FilterSupport1D::inRange(int rangeStart, int rangeEnd) const
{
    return start >= rangeStart && end <= rangeEnd;
}


inline FilterSupport1D intersect(const FilterSupport1D s1, const FilterSupport1D s2)
{
    return FilterSupport1D(std::max(s1.start, s2.start), std::min(s1.end, s2.end));
}

//------------------------------------------------------------------------------
// FilterSupport

inline FilterSupport::FilterSupport(int startX, int endX,
        int startY, int endY)
    : sx(startX, endX),
    sy(startY, endY)
{ }

inline FilterSupport::FilterSupport(const FilterSupport1D sx,
        const FilterSupport1D sy)
    : sx(sx),
    sy(sy)
{ }

inline int FilterSupport::area() const
{
    return sx.range()*sy.range();
}

inline bool FilterSupport::isEmpty() const
{
    return sx.isEmpty() || sy.isEmpty();
}

inline bool FilterSupport::inRange(
        int startX, int endX, int startY, int endY) const
{
    return sx.inRange(startX, endX) && sy.inRange(startY, endY);
}

inline bool FilterSupport::intersectsRange(int startX, int endX,
        int startY, int endY) const
{
    return sx.intersectsRange(startX, endX) && sy.intersectsRange(startY, endY);
}

inline FilterSupport intersect(const FilterSupport& s1, const FilterSupport& s2)
{
    return FilterSupport(intersect(s1.sx, s2.sx), intersect(s1.sy, s2.sy));
}

}
OIIO_NAMESPACE_EXIT

#endif // OPENIMAGEIO_FILTERSUPPORT_H
