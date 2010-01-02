/*!
 * \file
 * \brief Definitions of some specific functions useful in communications
 * \author Tony Ottosson and Erik G. Larsson
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
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

#ifndef COMMFUNC_H
#define COMMFUNC_H

#include <itpp/base/mat.h>
#include <itpp/base/vec.h>


namespace itpp
{

/*!
  \brief Generate Gray code of blocklength m.
  \ingroup misccommfunc

  The codes are contained as binary codewords {0,1} in the rows of the
  returned matrix.
  See also the \c gray() function in \c math/scalfunc.h.
*/
bmat graycode(int m);

/*!
  \brief Calculate the Hamming distance between \a a and \a b
  \ingroup misccommfunc
*/
int hamming_distance(const bvec &a, const bvec &b);

/*!
  \brief Calculate the Hamming weight of \a a
  \ingroup misccommfunc
*/
int weight(const bvec &a);

/*!
 * \brief Compute the water-filling solution
 * \ingroup misccommfunc
 *
 * This function computes the solution of the water-filling problem
 * \f[
 * \max_{p_0,...,p_{n-1}} \sum_{i=0}^{n-1} \log\left(1+p_i\alpha_i\right)
 * \f]
 * subject to
 * \f[
 * \sum_{i=0}^{n-1} p_i \le P
 * \f]
 *
 * \param alpha vector of \f$\alpha_0,...,\alpha_{n-1}\f$ gains (must have
 * strictly positive elements)
 * \param P power constraint
 * \return vector of power allocations \f$p_0,...,p_{n-1}\f$
 *
 * The computational complexity of the method is \f$O(n^2)\f$ at most
 */
vec waterfilling(const vec& alpha, double P);

} // namespace itpp

#endif // #ifndef COMMFUNC_H
