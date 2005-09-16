/*!
 * \file 
 * \brief Definitions of a vector quantizer training functions
 * \author Thomas Eriksson
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
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
 * -------------------------------------------------------------------------
 */

#ifndef VQTRAIN_H
#define VQTRAIN_H

#include <itpp/base/mat.h>
#include <itpp/base/array.h>

namespace itpp {

  //!
  double kmeansiter(Array<vec> &DB, mat &codebook);
  //!
  mat kmeans(Array<vec> &DB, int SIZE, int NOITER=9999, bool VERBOSE=true);
  //!
  mat lbg(Array<vec> &DB, int SIZE, int NOITER=9999, bool VERBOSE=true);


  /*! 
    \ingroup sourcecoding
    \brief Function for vector quantization training

    The following code illustrates how the VQ can be trained.

    \code
    VQ Quantizer;
    mat A;
    Array<vec> database;

    // read vectors into database somehow
    ...

    // train a vq
    A = vqtrain(database, 1024, 1000000);
    Quantizer.set_codebook(A);
    \endcode
  */
  mat vqtrain(Array<vec> &DB, int SIZE, int NOITER, double STARTSTEP=0.2, bool VERBOSE=true);
  //!
  vec sqtrain(const vec &inDB, int SIZE);

  //!
  ivec bitalloc(const vec& variances, int nobits); 
} // namespace itpp

#endif // #ifndef VQTRAIN_H
