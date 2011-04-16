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
/// \brief Utilities for working with Elliptical Weighted Average filters

#include "ewafilter.h"

OIIO_NAMESPACE_ENTER
{

namespace {

/// \brief Clamp the filter anisotropy to a maximum value.
///
/// \param covariance - Covariance matrix for a gaussian filter
/// \param minorAxisWidth - return the width of the minor ellipse axis in here.
/// \param minorAxisDirection - return direction of minor axis of ellipse
/// \param maxAspectRatio - maximum allowable aspect ratio for the filter.
/// \param logEdgeWeight - log of the filter weight at the filter edge.
///
inline void clampEccentricity(Matrix2D& covariance, float& minorAxisWidth,
                              Imath::V2f& minorAxisDirection,
                              float maxAspectRatio, float logEdgeWeight)
{
    // Eigenvalues of the covariance matrix are the squares of the lengths of
    // the semi-major and semi-minor axes of the filter support ellipse.
    // eigenvalues() guarentees that eig1 >= eig2.
    float eig1 = 1, eig2 = 1;
    covariance.eigenvalues(eig1, eig2);
    Matrix2D R = covariance.orthogDiagonalize(eig1, eig2);
    // by construction, covariance = R * D * R^T, where
    // D = Matrix2D(eig1, eig2)
    if(maxAspectRatio*maxAspectRatio*eig2 < eig1)
    {
        // Clamp eccentricity
        if(1000000*eig2 < eig1)
        {
            // Special case for numerical stability in the case of absolutely
            // extreme anisotropy: expand minor axis so that it's always
            // overblurred.
            eig2 = eig1/(maxAspectRatio*maxAspectRatio);
        }
        else
        {
            // Normal case: choose a length scale equal to the geometric mean
            // of lengths, and expand long axis + contract short axis about
            // this length to get a max anisotropic filter.  This results in
            // some blurring along the short axis, and some aliasing along the
            // long axis, but seems like a good compromise.
            float g = sqrt(eig2*eig1);
            eig1 = g*maxAspectRatio;
            eig2 = g/maxAspectRatio;
        }
        covariance = R * Matrix2D(eig1, eig2) * R.transpose();
    }
    minorAxisDirection = Imath::V2f(R.b, R.d);
    minorAxisWidth = std::sqrt(8*eig2*logEdgeWeight);
}

} // unnamed namespace

namespace detail {
// Lookup-table for exp(-x) for filter weighting.  20 Data points is probably
// about right, though it's possible to use even less.  Using 10
NegExpTable negExpTable(20, 6);

} // namespace detail

void EwaFilterFactory::computeFilter(
        const Imath::V2f& dstdx, const Imath::V2f& dstdy,
        float baseResS, float baseResT, const Matrix2D& blurVariance,
        float maxAspectRatio)
{
    // Get Jacobian of the inverse texture warp, and scale it by the resolution
    // of the base texture.
    Matrix2D invJ = Matrix2D (dstdx.x, dstdy.x,
                              dstdx.y, dstdy.y);

    // Reconstruction filter variance (conceptually the filter which
    // reconstructs a continuous image from the underlying discrete samples)
    //
    // A variance of 1/(2*pi) gives a filter with centeral weight 1, but in
    // practise this is slightly too small (resulting in a little bit of
    // aliasing).  Therefore it's adjusted up slightly.  In his original work,
    // Heckbert appears to recommend using 1 for these variances, but this
    // results in excess blurring.
    const float reconsVar = 1.3/(2*M_PI);
    // "Prefilter" variance (antialiasing filter which is used immediately
    // before resampling onto the discrete grid)
    const float prefilterVar = 1.3/(2*M_PI);

    // Construct the covariance matrix, coVar.
    //
    // The inverse Jacobian is combined with the blur variance matrix first.
    // The transpose is in the opposite position from Heckbert's thesis since
    // we're using a column-vector rather than row-vector convention.
    //
    // The RI spec (version 3.2, page 149) says that blur is "an additional
    // area to be added to the texture area filtered in both the s and t
    // directions, expressed in units of texture coordinates".  Though from
    // their example it's obvious that they meant length rather than area.
    //
    // The closure of gaussian filters under convolution provides a very nice
    // way of incorporating extra filter blurring in a systematic way:  We
    // imagine that any blur is an extra filter convoluted with the
    // reconstructed image.  For gaussian filters this convolution corresponds
    // to adding covariance matrices.
    //
    // Note the brackets around (invJ*invJ.transpose()) are important for
    // numerical stability to ensure symmetry in the result matrix.
    Matrix2D coVar = prefilterVar*(invJ*invJ.transpose()) + blurVariance;
    // Transform variance matrix to put it in base texture _raster_ units
    // if S = Matrix2D(baseResS, baseResT) then this is just
    // coVar = S*coVar*S;
    coVar.a *= baseResS*baseResS;
    coVar.b *= baseResS*baseResT;
    coVar.c *= baseResS*baseResT;
    coVar.d *= baseResT*baseResT;
    // Finally, add in the reconstruction variance.
    coVar += Matrix2D(reconsVar);

    // Clamp the eccentricity of the filter.
    clampEccentricity(coVar, m_minorAxisWidth, m_minorAxis, maxAspectRatio,
                      m_logEdgeWeight);
    // ensure that minor axis x-component is positive.
    if(m_minorAxis.x < 0)
        m_minorAxis = -m_minorAxis;
    // Get the quadratic form
    m_quadForm = 0.5*coVar.inv();
}

}
OIIO_NAMESPACE_EXIT

