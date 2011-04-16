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
/// \brief Utilities for working with Elliptical Weighted Average filters

#ifndef OPENIMAGEIO_EWAFILTER_H
#define OPENIMAGEIO_EWAFILTER_H

#include <vector>
#include <cmath>

#include <OpenEXR/ImathFun.h>
#include <OpenEXR/ImathVec.h>

#include "dassert.h"

#include "filtersupport.h"
#include "matrix2d.h"

OIIO_NAMESPACE_ENTER
{

//------------------------------------------------------------------------------
/// \brief A filter functor for evaluating 2D gaussian EWA filter weights.
///
/// This filter functor class is conveniently constructed by the
/// EwaFilterFactory class.  It's basically a 2D gaussian filter with a
/// cutoff, evaluated as described in the class documentation for
/// EwaFilterFactory, and repeated briefly here:
///
/// \verbatim
///   Q(x)  = (x-c)^T * Q * (x-c)
///   W(x)  =  / exp[-Q(x)]    for  Q(x) < logEdgeWeight
///            \ 0             elsewhere
/// \endverbatim
///
/// Where c is the filter center, Q is the quadratic form matrix and x is the
/// point at which the filter is being evaluated.  The total fraction of the
/// ideal filter weight which is lost outside the support is equal to
/// exp(-logEdgeWeight).
///
class EwaFilter
{
    public:
        /// Construct a 2D gaussian filter with the given coefficients.
        ///
        /// \param quadForm - quadratic form matrix
        /// \param filterCenter - the filter is centered on this point.
        /// \param logEdgeWeight - log of the filter weight at the edge cutoff.
        /// \param minorAxis - vector in minor axis direction
        /// \param boundOffset - displacement from center to right bound line
        ///
        EwaFilter(Matrix2D quadForm, Imath::V2f filterCenter,
                  float logEdgeWeight, Imath::V2f minorAxis,
                  Imath::V2f boundOffset);

        /// EWA filters are never pre-noramlized; return false.
        static bool isNormalized() { return false; }

        /// \brief Evaluate the filter at the given point in image space.
        ///
        /// \param x
        /// \param y - these parameters are the position in raster space which
        ///            the filter weight is to be calculated at.  (So for an
        ///            image of size WxH, x and y would normally lie somewhere
        ///            in the ranges [0,W] and [0,H] respectively, thought they
        ///            don't have to)
        ///
        float operator()(float x, float y) const;
        /// Get the extent of the filter in integer raster coordinates.
        FilterSupport support() const;

        /// Translate centre of filter in the image plane.
        void translate(Imath::V2f trans);

        /// First x coordinate inside filter at height y.
        ///
        /// Diagonally anisotropic filters don't fit well into an axis-aligned
        /// bounding box as shown here (x's represent nonzero filter values):
        ///
        ///   +---------+
        ///   |      / x|
        ///   |     / x |
        ///   |    /xx /|
        ///   |   /xxx/ |
        ///   |  /xxx/  |
        ///   | /xxx/   |
        ///   |/ xx/    |
        ///   | x /     |
        ///   |x /      |
        ///   +---------+
        ///
        /// This means that the filter inside the bounding box provided by
        /// support() contains a lot of zeros, so it's inefficient to iterate
        /// over this entire region.  Instead, we can define a bounding slab
        /// by a pair of parallel lines, and iterate between those instead.
        /// xbegin() and xend() calculate these lines at a given height y.
        int xbegin(float y) const
        {
            return Imath::ceil(m_leftBound_x0 + m_bound_dxdy*y);
        }
        /// One-past-last x coordinate inside filter at height y.
        int xend(float y) const
        {
            return Imath::ceil(m_rightBound_x0 + m_bound_dxdy*y);
        }

    private:
        /// Quadratic form matrix
        Matrix2D m_quadForm;
        /// Center point of the gaussian filter function
        Imath::V2f m_filterCenter;
        /// The log of the filter weight at the filter edge cutoff.
        float m_logEdgeWeight;
        /// Coefficients for bounding slab.
        float m_bound_dxdy;
        float m_leftBound_x0;
        float m_rightBound_x0;
};

//------------------------------------------------------------------------------
/// \brief A class encapsulating Elliptically Weighted Average (EWA) filter
/// weight computation.
///
/// EWA filtering is based on the convolution of several gaussian filters and
/// composition with the linear approximation to the image warp at the
/// sampling point.  The original theory was developed by Paul Heckbert et al.
/// and is well-explained in his masters thesis, "Fundamentals of Texture
/// Mapping and Image Warping", which may be found at
/// http://www.cs.cmu.edu/~ph/.  "Physically Based Rendering" also has a small
/// section mentioning EWA.
///
/// The derivation of an EWA filter may be broken into four conceptual stages:
///   1) Reconstruct a continuous image from the discrete samples using a
///      gaussian *reconstruction filter*.  The reconstruction filter support
///      should be wide enough to avoid samples "falling between" the discrete
///      points of the source image.
///   2) Apply a texture blur filter to this continuous image.
///   3) Warp the 2D image into another part of the 2D plane with an
///      arbitrary transformation.  This effect of the transformation on the
///      filter kernel is approximated by the local linear approximation (ie,
///      the Jacobian).
///   4) Filter the resulting image with a gaussian "prefilter" before
///      converting back to discrete samples.  The support of the prefilter
///      should be wide enough to remove aliasing.
///
/// The neat thing about EWA is that the full filter which you get after
/// putting the steps 1-4 together is just another gaussian filter acting on
/// the original image.  This means we can implement it by simple iteration
/// over a box in the source image which contains the filter support.  The
/// result may be computed deterministically and hence suffers from no
/// sampling error.  The inclusion of the linear transformation implies good
/// anisotropic filtering characteristics.
///
/// A 2D gaussian filter may be defined by a quadratic form,
///
/// \verbatim
///   Q(x,y) = a*x^2 + (b+c)*x*y + d*y^2,
/// \endverbatim
///
/// such that the filter weights are given by
///
/// \verbatim
///   W(x,y) = exp(-Q(x,y)).
/// \endverbatim
///
/// This filter has infinite support; in practise we need to truncate it at
/// some point, ideally such that we leave only a fixed fraction of the filter
/// weight outside the truncated region.  To do this, we choose a cutoff
/// parameter, C for the weight function, and set all weights less than that
/// cutoff to zero.  The newly truncated filter weight, W', looks like:
///
/// \verbatim
///   W'(x,y) = W(x,y)     for W(x,y) >= C
///           = 0          for W(x,y) < C
/// \endverbatim
///
/// A curious and convenient feature of such gaussian filters in two
/// dimensions (and only in two dimensions!) is that the fraction of the
/// *total* weight which we're ignoring by doing this is just (1-C).  (To see
/// this is true, do the gaussian integrals inside and outside the cutoff.)
/// In practise we work with the quantity logEdgeWeight = -ln(C) below.
////
class EwaFilterFactory
{
    public:
        /// \brief Perform EWA filter weight setup
        ///
        /// This initializes the filter to be used over a texture of
        /// resolution baseResS x baseResT.  This is assumed to be the maximum
        /// resolution texture over which we will want to use the filter. The
        /// filter can be adjusted for other lower resolutions using the
        /// function adjustTextureScale().
        ///
        /// \param st - texture coordinates
        /// \param dstdx,dstdy - derivatives of texture coordinates:
        ///            approximate preimage of an output pixel box under the
        ///            image warp.  dstdx and dstdy give the linear
        ///            approximation to the image warp at the centre (that is,
        ///            they represent the Jacobian of the mapping).
        /// \param baseResS - width of the base texture (used to determine a
        ///            minimum reconstruction filter variance)
        /// \param baseResT - height of the base texture (used to determine a
        ///            minimum reconstruction filter variance)
        /// \param blurVariance - Variance matrix for additional filter blur
        ///            (see ewaBlurMatrix() )
        /// \param tBlur - Additional filter blur in the t-direction
        /// \param logEdgeWeight - Related to the total fraction of the ideal
        ///            filter weight, which is equal to exp(-logEdgeWeight).
        /// \param maxAspectRatio - maximum anisotropy at which the filter will
        ///            be clamped.
        ///
        EwaFilterFactory(const Imath::V2f& st,
                         const Imath::V2f& dstdx, const Imath::V2f& dstdy,
                         float baseResS, float baseResT,
                         const Matrix2D& blurVariance,
                         float logEdgeWeight = 3, 
                         float maxAspectRatio = 20);

        /// \brief Create an EWA filter for a given miplevel
        ///
        /// \param width - width of the miplevel
        /// \param height - height of the miplevel
        ///
        EwaFilter createFilter(int width, int height) const;

        /// Get the width of the filter along the minor axis of the ellipse
        float minorAxisWidth() const;
    private:
        /// \brief Compute and cache EWA filter coefficients
        ///
        /// This function initializes the appropriate matrix for the quadratic
        /// form,
        ///   Q = [a b]
        ///       [c d]
        /// which represents the EWA filter over the quadrilateral given by
        /// sampleQuad.  Q is cached in m_quadForm.  The minimum width along
        /// the minor axis of the filter is cached in m_minorAxisWidth.
        ///
        /// For parameters, see the corresponding ones in the
        /// EwaFilterFactory constructor.
        ///
        ///
        void computeFilter(const Imath::V2f& dstdx, const Imath::V2f& dstdy,
                           float baseResS, float baseResT,
                           const Matrix2D& blurVariance, float maxAspectRatio);

        /// Base resolution of the texture
        int m_baseResS;
        int m_baseResT;
        /// Quadratic form matrix
        Matrix2D m_quadForm;
        /// Center point of the gaussian filter function
        Imath::V2f m_filterCenter;
        /// The log of the filter weight at the filter edge cutoff.
        float m_logEdgeWeight;
        /// Width of the semi-minor axis of the elliptical filter
        float m_minorAxisWidth;
        /// Direction of minor axis of ellipse.  Always has m_minorAxis.x > 0.
        Imath::V2f m_minorAxis;
};


/// \brief Compute the blur variance matrix for axis-aligned blur.
///
/// The returned matrix gives an appropriate blur ellipse aligned with
/// the x and y axes.
///
/// \param sBlur - blur in x-direction
/// \param tBlur - blur in y-direction
///
Matrix2D ewaBlurMatrix(float sBlur, float tBlur);


//==============================================================================
// Implementation details
//==============================================================================
// EwaFilterFactory implementation
inline EwaFilterFactory::EwaFilterFactory(const Imath::V2f& st,
                                          const Imath::V2f& dstdx,
                                          const Imath::V2f& dstdy,
                                          float baseResS, float baseResT,
                                          const Matrix2D& blurVariance,
                                          float logEdgeWeight, 
                                          float maxAspectRatio)
    : m_baseResS(baseResS),
    m_baseResT(baseResT),
    m_quadForm(0),
    m_filterCenter(st),
    m_logEdgeWeight(logEdgeWeight),
    m_minorAxisWidth(0)
{
    // Scale the filterCenter up to the dimensions of the base texture.
    m_filterCenter.x = m_filterCenter.x*baseResS;
    m_filterCenter.y = m_filterCenter.y*baseResT;
    // compute and cache the filter
    computeFilter(dstdx, dstdy, baseResS, baseResT, blurVariance,
                  maxAspectRatio);
}

inline Matrix2D ewaBlurMatrix(float sBlur, float tBlur)
{
    if(sBlur > 0 || tBlur > 0)
    {
        // The factor blurScale gives an amount of blur which is roughly
        // consistent with that used by PRMan (and 3delight) for the same
        // scenes.
        const float blurScale = 0.5f;
        float sStdDev = sBlur*blurScale;
        float tStdDev = tBlur*blurScale;
        return Matrix2D(sStdDev*sStdDev, tStdDev*tStdDev);
    }
    else
        return Matrix2D(0);
}

inline EwaFilter EwaFilterFactory::createFilter(int width, int height) const
{
    // Displacement from filter center to edge of filter support on minor axis.
    Imath::V2f offset = 0.5f*m_minorAxisWidth*m_minorAxis;
    if (width == m_baseResS && height == m_baseResT) {
        // Special case for base level.  Yep, this actually can have a minor
        // but measurable impact on speed in some cases (3% say).
        return EwaFilter(m_quadForm, m_filterCenter - Imath::V2f(0.5f),
                         m_logEdgeWeight, m_minorAxis, offset);
    }
    // Transform the filter coefficients that are relevant to the base
    // texture into filter coeffs for the miplevel.
    Imath::V2f scale = Imath::V2f(float(width) / m_baseResS,
                                  float(height) / m_baseResT);
    Imath::V2f isc = Imath::V2f(1.0f)/scale;
    // Compute new filter center.  Note well that we define the center of the
    // first pixel as occuring a distance of 0.5 pixels from the NDC origin,
    // hence the offset of -0.5.
    Imath::V2f newCenter = scale*m_filterCenter - Imath::V2f(0.5f);
    // The strange-looking matrix that is passed into the EwaFilter
    // constructor below is simply the hand-written version of the following M:
    //
    // Matrix2D S(1/isc.x, 1/isc.y);
    // M = scaleMatrix*m_quadForm*scaleMatrix;
    //
    // Note that axes of the ellipse transform like normals.
    return EwaFilter(
        Matrix2D(isc.x*isc.x*m_quadForm.a, isc.x*isc.y*m_quadForm.b,
                 isc.x*isc.y*m_quadForm.c, isc.y*isc.y*m_quadForm.d),
        newCenter, m_logEdgeWeight, m_minorAxis*isc, offset*scale
    );
}

inline float EwaFilterFactory::minorAxisWidth() const
{
    return m_minorAxisWidth;
}


//------------------------------------------------------------------------------
namespace detail {

/// A lookup table for std::exp(-x).
///
/// The lookup is via operator() which does linear interpolation between
/// tabulated values.
///
class NegExpTable
{
    private:
        std::vector<float> m_values;
        float m_invRes;
        float m_rangeMax;
    public:
        /// \brief Construct the lookup table
        /// \param numPoints - number of points in the table.
        /// \param rangeMax - maximum value of the input variable that the
        ///                   table should be computed for.  Inputs larger
        ///                   than or equal to this will return 0.
        NegExpTable(int numPoints, float rangeMax)
            : m_values(),
            m_invRes((numPoints-1)/rangeMax),
            m_rangeMax(rangeMax)
        {
            float res = 1/m_invRes;
            m_values.resize(numPoints);
            for(int i = 0; i < numPoints; ++i)
            {
                m_values[i] = exp(-i*res);
            }
        }

        /// \brief Look up an approximate exp(-x) for x > 0
        ///
        /// This does linear interpolation between x values.
        /// \param x - 
        float operator()(float x) const
        {
            /// \todo Optimization: Possibly may remove some of the checks here.
            if(x >= m_rangeMax)
                return 0;
            float xRescaled = x*m_invRes;
            int index = Imath::floor(xRescaled);
            DASSERT(index >= 0);
            float interp = xRescaled - index;
            return (1-interp)*m_values[index] + interp*m_values[index+1];
        }
};
extern NegExpTable negExpTable;

} // namespace detail


//------------------------------------------------------------------------------
// EwaFilter implementation
inline EwaFilter::EwaFilter(Matrix2D quadForm, Imath::V2f filterCenter,
                            float logEdgeWeight, Imath::V2f minorAxis,
                            Imath::V2f boundOffset)
    : m_quadForm(quadForm),
    m_filterCenter(filterCenter),
    m_logEdgeWeight(logEdgeWeight)
{
    if(fabs(minorAxis.x) > 0.001*fabs(minorAxis.y))
    {
        // slope of bounding lines
        m_bound_dxdy = -minorAxis.y/minorAxis.x;
        // p1,p2 = points on left,right bounding lines
        Imath::V2f p1 = m_filterCenter - boundOffset;
        Imath::V2f p2 = m_filterCenter + boundOffset;
        // constant offsets for left and right bounding lines.
        m_leftBound_x0  = p1.x - m_bound_dxdy*p1.y;
        m_rightBound_x0 = p2.x - m_bound_dxdy*p2.y;
    }
    else
    {
        // Special case: bounding slab is basically parallel to the x-axis.  In
        // this case, give up doing a careful bounding job - the rectangular
        // filter support will work well.
        m_bound_dxdy = 0;
        m_leftBound_x0 = INT_MIN/4;
        m_rightBound_x0 = INT_MAX/4;
    }
}

inline float EwaFilter::operator()(float x, float y) const
{
    x -= m_filterCenter.x;
    y -= m_filterCenter.y;
    // evaluate quadratic form
    float q = m_quadForm.a*x*x + (m_quadForm.b+m_quadForm.c)*x*y
              + m_quadForm.d*y*y;
    // Check whether we're inside the filter cutoff; if so use a lookup table
    // to get the filter weight.  Using a lookup table rather than directly
    // using std::exp() results in very large speedups, since the filter
    // weights are needed inside the inner loop
    if(q < m_logEdgeWeight)
        return detail::negExpTable(q);
    return 0;
}

inline FilterSupport EwaFilter::support() const
{
    float detQ = m_quadForm.det();
    // Compute filter radii
    float sRad = std::sqrt(m_quadForm.d*m_logEdgeWeight/detQ);
    float tRad = std::sqrt(m_quadForm.a*m_logEdgeWeight/detQ);
    return FilterSupport(
            Imath::ceil(m_filterCenter.x-sRad),     // startX
            Imath::floor(m_filterCenter.x+sRad)+1,  // endX
            Imath::ceil(m_filterCenter.y-tRad),     // startY
            Imath::floor(m_filterCenter.y+tRad)+1   // endY
        );
}

inline void EwaFilter::translate(Imath::V2f trans)
{
    m_filterCenter += trans;
    m_leftBound_x0  += trans.x - m_bound_dxdy*trans.y;
    m_rightBound_x0 += trans.x - m_bound_dxdy*trans.y;
}

}
OIIO_NAMESPACE_EXIT

#endif // OPENIMAGEIO_EWAFILTER_H
