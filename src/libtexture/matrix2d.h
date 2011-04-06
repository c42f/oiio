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
/// A minimal 2D matrix class.

#ifndef OPENIMAGEIO_MATRIX2D_H
#define OPENIMAGEIO_MATRIX2D_H

#include <cmath>

#include "dassert.h"

OIIO_NAMESPACE_ENTER
{


/// A minimal 2D matrix class.
struct Matrix2D
{
    /// \brief Matrix components:
    ///
    /// The matrix looks like
    ///
    ///   [a b]
    ///   [c d]
    ///
    float a;
    float b;
    float c;
    float d;

    /// Construct a multiple of the identity.
    inline Matrix2D(float diag);
    /// Construct a diagonal matrix.
    inline Matrix2D(float a, float d);
    /// Construct a general matrix.
    inline Matrix2D(float a, float b, float c, float d);

    /// \name Addition operators
    //@{
    inline Matrix2D operator+(const Matrix2D& rhs) const;
    inline Matrix2D operator+(float f) const;
    friend inline Matrix2D operator+(float f, const Matrix2D& mat);

    inline Matrix2D& operator+=(const Matrix2D& rhs);
    //@}

    /// \name Multiplication operators
    //@{
    inline Matrix2D operator*(const Matrix2D& rhs) const;
    inline Matrix2D operator*(float mult) const;
    friend inline Matrix2D operator*(float mult, const Matrix2D& mat);
    //@}

    /// \name Properties and transformations of the matrix.
    //@{
    /// Return the inverse
    inline Matrix2D inv() const;
    /// Return the determinant.
    inline float det() const;
    /// Return the matrix transpose
    inline Matrix2D transpose() const;
    /// \brief Get the eigenvalues of the matrix
    ///
    /// If the eigenvalues are complex rather than real, trigger an assert in
    /// debug mode.  In release mode, return the real part of the eigenvalues.
    ///
    /// The eigenvalues are guarenteed to be ordered such that that l1 >= l2.
    ///
    /// \param l1 - first eigenvalue (output variable)
    /// \param l2 - second eigenvalue (output variable)
    ///
    inline void eigenvalues(float& l1, float& l2) const;

    /// \brief Attempt to get a matrix which orthogonally diagonalizes this
    /// matrix.
    ///
    /// That is, find a matrix R such that if A is the matrix to diagonalize,
    ///
    ///   A = R * D * R^T
    ///
    /// where R is an orthogonal matrix, and D is the diagonal matrix,
    /// diag(l1,l2).
    ///
    /// \param l1 - first eigenvalue
    /// \param l2 - the second eigenvalue.
    ///
    /// The somewhat nasty behaviour of feeding eigenvalues are back into this
    /// function is used for efficiency, since we may want to compute the
    /// eigenvalues seperately, but at the same time avoid computing them
    /// again when computing the diagonalizing matrix.
    ///
    /// The behaviour of this function is undefined if the matrix is not
    /// symmetric.
    ///
    inline Matrix2D orthogDiagonalize(float l1, float l2) const;
    //@}
};

//==============================================================================
// Implentation details
//==============================================================================
// Matrix2D implementation

// Constructors
inline Matrix2D::Matrix2D(float diag)
    : a(diag), b(0), c(0), d(diag)
{ }
inline Matrix2D::Matrix2D(float a, float d)
    : a(a), b(0), c(0), d(d)
{ }
inline Matrix2D::Matrix2D(float a, float b, float c, float d)
    : a(a), b(b), c(c), d(d)
{ }

// Addition
inline Matrix2D Matrix2D::operator+(const Matrix2D& rhs) const
{
    return Matrix2D(a+rhs.a, b+rhs.b, c+rhs.c, d+rhs.d);
}
inline Matrix2D Matrix2D::operator+(float f) const
{
    return Matrix2D(a+f, b, c, d+f);
}
inline Matrix2D operator+(float f, const Matrix2D& mat)
{
    return mat+f;
}
inline Matrix2D& Matrix2D::operator+=(const Matrix2D& rhs)
{
    a += rhs.a;
    b += rhs.b;
    c += rhs.c;
    d += rhs.d;
    return *this;
}

// Multiplication
inline Matrix2D Matrix2D::operator*(const Matrix2D& rhs) const
{
    return Matrix2D(
            a*rhs.a + b*rhs.c, a*rhs.b + b*rhs.d,
            c*rhs.a + d*rhs.c, c*rhs.b + d*rhs.d
            );
}
inline Matrix2D Matrix2D::operator*(float mult) const
{
    return Matrix2D(a*mult, b*mult, c*mult, d*mult);
}
inline Matrix2D operator*(float mult, const Matrix2D& mat)
{
    return mat*mult;
}

// Inverse, determinant, etc.
inline Matrix2D Matrix2D::inv() const
{
    // There's a simple formula for the inverse of a 2D matrix.  We use this here.
    float D = det();
    DASSERT(D != 0);
    if(D != 0)
        return Matrix2D(d/D, -b/D, -c/D, a/D);
    else
        return Matrix2D(1);
}
inline float Matrix2D::det() const
{
    return a*d - b*c;
}
inline Matrix2D Matrix2D::transpose() const
{
    return Matrix2D(a, c, b, d);
}

inline void Matrix2D::eigenvalues(float& l1, float& l2) const
{
    // Special-case formula for eigenvalues of a 2D matrix.  It simply boils
    // down to solving the quadratic equation
    //
    // l^2 - Tr(A)*l + det(A) = 0
    //
    // for l.
    float firstTerm = (a+d)*0.5;
    float secondTerm = (a-d)*(a-d) + 4*c*b;
    DASSERT(secondTerm >= -std::numeric_limits<float>::epsilon());
    // For robustness, set secondTerm = 0 if it's negative.  This will get the
    // real part of the result.
    if(secondTerm < 0)
        secondTerm = 0;
    secondTerm = std::sqrt(secondTerm)*0.5;
    l1 = firstTerm + secondTerm;
    l2 = firstTerm - secondTerm;
}

inline Matrix2D Matrix2D::orthogDiagonalize(float l1, float l2) const
{
    // As usual, we construct the matrix from the two orthonormal eigenvectors.
    // These eigenvectors only exist if the matrix is symmetric, so assert
    // symmetry:
    DASSERT(std::fabs((b - c)) <= 1e-5*std::fabs(c) ||
            std::fabs((b - c)) <= 1e-5*std::fabs(b) );
    if(l1 == l2)
    {
        // Special (easy) case for degenerate eigenvalues
        return Matrix2D(1, 0,
                          0, 1);
    }
    // Calculate first eigenvector.  [u1, u2] and [v1,v2] are two alternatives
    // for the eigenvector - we take the one with the larger length for
    // maximum numerical stability.
    float u1 = b;
    float u2 = l1-a;
    float uLen2 = u1*u1 + u2*u2;
    float v1 = l1-d;
    float v2 = c;
    float vLen2 = v1*v1 + v2*v2;
    if(vLen2 > uLen2)
    {
        u1 = v1;
        u2 = v2;
        uLen2 = vLen2;
    }
    float invLenU = 1/std::sqrt(uLen2);
    u1 *= invLenU;
    u2 *= invLenU;
    return Matrix2D(u1, -u2,
                      u2, u1);
}


}
OIIO_NAMESPACE_EXIT

#endif // OPENIMAGEIO_MATRIX2D_H
