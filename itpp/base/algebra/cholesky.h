/*!
 * \file 
 * \brief Definitions of Cholesky factorisation functions
 * \author Tony Ottosson
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2006  (see AUTHORS file for a list of contributors)
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

#ifndef CHOLESKY_H
#define CHOLESKY_H

#include <itpp/base/mat.h>


namespace itpp {

  /*! \addtogroup matrixdecomp
   */
  //!@{

  /*!
    \brief Cholesky factorisation of real symmetric and positive definite matrix

    The Cholesky factorisation of a real symmetric positive-definite matrix \f$\mathbf{X}\f$
    of size \f$n \times n\f$ is given by
    \f[
    \mathbf{X} = \mathbf{F}^T \mathbf{F}
    \f]
    where \f$\mathbf{F}\f$ is an upper trangular \f$n \times n\f$ matrix.

    Returns true if calcuation succeeded. False otherwise.
  */
  bool chol(const mat &X, mat &F);

  /*!
    \brief Cholesky factorisation of real symmetric and positive definite matrix

    The Cholesky factorisation of a real symmetric positive-definite matrix \f$\mathbf{X}\f$
    of size \f$n \times n\f$ is given by
    \f[
    \mathbf{X} = \mathbf{F}^T \mathbf{F}
    \f]
    where \f$\mathbf{F}\f$ is an upper trangular \f$n \times n\f$ matrix.
  */
  mat chol(const mat &X);


  /*!
    \brief Cholesky factorisation of complex hermitian and positive-definite matrix

    The Cholesky factorisation of a hermitian positive-definite matrix \f$\mathbf{X}\f$
    of size \f$n \times n\f$ is given by
    \f[
    \mathbf{X} = \mathbf{F}^H \mathbf{F}
    \f]
    where \f$\mathbf{F}\f$ is an upper trangular \f$n \times n\f$ matrix.

    Returns true if calcuation succeeded. False otherwise.

    If \c X is positive definite, true is returned and \c F=chol(X)
    produces an upper triangular \c F. If also \c X is symmetric then \c F'*F = X.
    If \c X is not positive definite, false is returned.
  */
  bool chol(const cmat &X, cmat &F);

  /*!
    \brief Cholesky factorisation of complex hermitian and positive-definite matrix

    The Cholesky factorisation of a hermitian positive-definite matrix \f$\mathbf{X}\f$
    of size \f$n \times n\f$ is given by
    \f[
    \mathbf{X} = \mathbf{F}^H \mathbf{F}
    \f]
    where \f$\mathbf{F}\f$ is an upper trangular \f$n \times n\f$ matrix.
  */
  cmat chol(const cmat &X);

  //!@}

} // namespace itpp

#endif // #ifndef CHOLESKY_H
