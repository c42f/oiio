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
/// \brief Sampling quadrilateral struct definition

#ifndef OPENIMAGEIO_SAMPLEQUAD_H
#define OPENIMAGEIO_SAMPLEQUAD_H

#include <cmath>
#include <ImathVec.h>

OIIO_NAMESPACE_ENTER
{


//------------------------------------------------------------------------------
/// \brief 2D parallelogram over which to sample a texture
///
/// The parallelogram is defined by its center point along with vectors along
/// its sides.  As a picture:
///
/// \verbatim
///
///         _-------
///         /|    /
///     s2 /  c  /
///       /     /
///      .----->
///        s1
///
/// \endverbatim
///
/// This is a special case of SqSampleQuad, holding less information, but more
/// specifically adapted to some types of filtering such as EWA.
///
struct SamplePllgram
{
    /// center point for the sample
    Imath::V2f c;
    /// first side of parallelogram
    Imath::V2f s1;
    /// second side of parallelogram
    Imath::V2f s2;

    /// Trivial constructor
    SamplePllgram(const Imath::V2f& c, const Imath::V2f& s1, const Imath::V2f s2);

    /// \brief Remap the parallelogram to lie in the box [0,1]x[0,1]
    ///
    /// The remapping occurs by translating the center by an integer lattice
    /// point (n,m) such that c lies inside [0,1]x[0,1].
    ///
    /// \param xPeriodic - if true, remap in the x direction.
    /// \param yPeriodic - if true, remap in the y direction.
    ///
    void remapPeriodic(bool xPeriodic, bool yPeriodic);

    /// \brief Scale the parallelogram about its center point.
    ///
    /// \param xWidth - amount to expand sample quad in the x direction
    /// \param yWidth - amount to expand sample quad in the y direction
    ///
    void scaleWidth(float xWidth, float yWidth);
};



//==============================================================================
// Implementation details
//==============================================================================
// SamplePllgram implementation
inline SamplePllgram::SamplePllgram(const Imath::V2f& c, const Imath::V2f& s1,
        const Imath::V2f s2)
    : c(c),
    s1(s1),
    s2(s2)
{ }

inline void SamplePllgram::remapPeriodic(bool xPeriodic, bool yPeriodic)
{
    if(xPeriodic || yPeriodic)
    {
        if(c.x < 0 || c.y < 0 || c.x >= 1 || c.y >= 1)
            c -= Imath::V2f(std::floor(c.x), std::floor(c.y));
    }
}

inline void SamplePllgram::scaleWidth(float xWidth, float yWidth)
{
    if(xWidth != 1 || yWidth != 1)
    {
        s1.x = s1.x*xWidth;
        s1.y = s1.y*yWidth;
        s2.x = s1.x*xWidth;
        s2.y = s1.y*yWidth;
    }
}

//------------------------------------------------------------------------------

}
OIIO_NAMESPACE_EXIT

#endif // SAMPLEQUAD_H_INCLUDED
