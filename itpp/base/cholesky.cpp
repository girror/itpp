/*!
 * \file 
 * \brief Implementation of Cholesky factorisation functions
 * \author Tony Ottosson
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

#include <algorithm>
#include <cassert>
#include <itpp/config.h>
#include <itpp/itconfig.h>
#include <itpp/base/cholesky.h>

namespace itpp { 

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

#include <itpp/base/lapack.h>

  bool chol(const mat &X, mat &F)
  {

    char uplo='U';
    int n, lda, info;
    n = lda = X.rows();

    F = X; // input matrix is overwritten

    dpotrf_(&uplo, &n, F._data(), &lda, &info);

    // Set lower part to zero
    for (int i=0; i<n; i++)
      for(int j=i+1; j<n; j++)
	F(j,i) = 0;

    return (info==0);
  }

  bool chol(const cmat &X, cmat &F)
  {
    char uplo='U';
    int n, lda, info;
    n = lda = X.rows();

    F = X; // input matrix is overwritten

    zpotrf_(&uplo, &n, F._data(), &lda, &info);

    // Set lower part to zero
    for (int i=0; i<n; i++)
      for(int j=i+1; j<n; j++)
	F(j,i) = 0;

    return (info==0);
  }

#else // HAVE_LAPACK or HAVE_MKL

  bool chol(const mat &X, mat &F)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for chol() to exist");
    return false;
  }

  bool chol(const cmat &X, cmat &F)
  {

    it_error("You need to compile IT++ with LAPACK or MKL for chol() to exist");
    return false;
  }

#endif // HAVE_LAPACK or HAVE_MKL

  cmat chol(const cmat &X)
  {
    cmat F;
    if (!chol(X, F)) {
      it_warning("cholesky factorization didn't succeed");
    }

    return F;
  }

  mat chol(const mat &X)
  {
    mat F;
    if (!chol(X, F)) {
      it_warning("cholesky factorization didn't succeed");
    }

    return F;
  }

} // namespace itpp
