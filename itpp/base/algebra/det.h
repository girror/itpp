/*!
 * \file
 * \brief Definitions of determinant calculations
 * \author Tony Ottosson
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2008  (see AUTHORS file for a list of contributors)
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

#ifndef DET_H
#define DET_H

#include <itpp/base/mat.h>


namespace itpp
{

/*!
  \brief Determinant of real square matrix.
  \ingroup determinant

  Calculate determinant of the real matrix \f$\mathbf{X}\f$

  Uses LU-factorisation.
  \f[
  \det(\mathbf{X}) = \det(\mathbf{P}^T \mathbf{L}) \det(\mathbf{U}) = \det(\mathbf{P}^T) \prod(\mathrm{diag}(\mathbf{U}))
  \f]
  and the determinant of the permuation matrix is \f$ \pm 1\f$ depending on the number of row permutations
*/
double det(const mat &X);


/*!
  \brief Determinant of complex square matrix.
  \ingroup determinant

  Calculate determinant of the complex matrix \f$\mathbf{X}\f$

  Uses LU-factorisation.
  \f[
  \det(\mathbf{X}) = \det(\mathbf{P}^T \mathbf{L}) \det(\mathbf{U}) = \det(\mathbf{P}^T) \prod(\mathrm{diag}(\mathbf{U}))
  \f]
  and the determinant of the permuation matrix is \f$ \pm 1\f$ depending on the number of row permutations
*/
std::complex<double> det(const cmat &X);


} // namespace itpp

#endif // #ifndef DET_H
