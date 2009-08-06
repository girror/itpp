/*!
 * \file
 * \brief Definitions of LU factorisation functions
 * \author Tony Ottosson
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2009  (see AUTHORS file for a list of contributors)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * -------------------------------------------------------------------------
 */

#ifndef LU_H
#define LU_H

#include <itpp/base/mat.h>


namespace itpp
{


/*! \addtogroup matrixdecomp
 */
//!@{
/*!
  \brief LU factorisation of real matrix

  The LU factorization of the real matrix \f$\mathbf{X}\f$ of size \f$n \times n\f$ is given
  by
  \f[
  \mathbf{X} = \mathbf{P}^T \mathbf{L} \mathbf{U} ,
  \f]
  where \f$\mathbf{L}\f$ and \f$\mathbf{U}\f$ are lower and upper triangular matrices
  and \f$\mathbf{P}\f$ is a permutation matrix.

  The interchange permutation vector \a p is such that \a k and \a p(k) should be
  changed for all \a k. Given this vector a permutation matrix can be constructed using the
  function
  \code
  bmat permutation_matrix(const ivec &p)
  \endcode

  If \a X is an \a n by \a n matrix \a lu(X,L,U,p) computes the LU decomposition.
  \a L is a lower triangular, \a U an upper triangular matrix.
  \a p is the interchange permutation vector such that \a k and \a p(k) should be
  changed for all \a k.

  Returns true is calculation succeeds. False otherwise.
*/
bool lu(const mat &X, mat &L, mat &U, ivec &p);


/*!
  \brief LU factorisation of real matrix

  The LU factorization of the complex matrix \f$\mathbf{X}\f$ of size \f$n \times n\f$ is given
  by
  \f[
  \mathbf{X} = \mathbf{P}^T \mathbf{L} \mathbf{U} ,
  \f]
  where \f$\mathbf{L}\f$ and \f$\mathbf{U}\f$ are lower and upper triangular matrices
  and \f$\mathbf{P}\f$ is a permutation matrix.

  The interchange permutation vector \a p is such that \a k and \a p(k) should be
  changed for all \a k. Given this vector a permutation matrix can be constructed using the
  function
  \code
  bmat permutation_matrix(const ivec &p)
  \endcode

  If \a X is an \a n by \a n matrix \a lu(X,L,U,p) computes the LU decomposition.
  \a L is a lower triangular, \a U an upper triangular matrix.
  \a p is the interchange permutation vector such that elements \a k and row \a p(k) should be
  interchanged.

  Returns true is calculation succeeds. False otherwise.
*/
bool lu(const cmat &X, cmat &L, cmat &U, ivec &p);


//! Makes swapping of vector b according to the interchange permutation vector p.
void interchange_permutations(vec &b, const ivec &p);

//! Make permutation matrix P from the interchange permutation vector p.
bmat permutation_matrix(const ivec &p);
//!@}

} // namespace itpp

#endif // #ifndef LU_H
