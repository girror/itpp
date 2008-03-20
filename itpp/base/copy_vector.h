/*!
 * \file
 * \brief Vector copy functions for internal use
 * \author Tony Ottosson and Adam Piatyszek
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

#ifndef COPY_VECTOR_H
#define COPY_VECTOR_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#if defined (HAVE_BLAS)
#  include <itpp/base/blas.h>
#endif

#include <itpp/base/binary.h>
#include <cstring>


//! \cond

namespace itpp {


  /*
    Copy vector x to vector y. Both vectors are of size n
  */
  inline void copy_vector(const int n, const int *x, int *y) { memcpy(y, x, (unsigned int)n*sizeof(int)); }
  inline void copy_vector(const int n, const short *x, short *y) { memcpy(y, x, (unsigned int)n*sizeof(short)); }
  inline void copy_vector(const int n, const bin *x, bin *y) { memcpy(y, x, (unsigned int)n*sizeof(bin)); }
  inline void copy_vector(const int n, const float *x, float *y) { memcpy(y, x, (unsigned int)n*sizeof(float)); }
  inline void copy_vector(const int n, const std::complex<float> *x, std::complex<float> *y) { memcpy(y, x, (unsigned int)n*sizeof(std::complex<float>)); }

#if defined (HAVE_BLAS)
  inline void copy_vector(const int n, const double *x, double *y)
  {
    int incr = 1;
    blas::dcopy_(&n, x, &incr, y, &incr);
  }
  inline void copy_vector(const int n, const std::complex<double> *x,
			  std::complex<double> *y)
  {
    int incr = 1;
    blas::zcopy_(&n, x, &incr, y, &incr);
  }
#else
  inline void copy_vector(const int n, const double *x, double *y) { memcpy(y, x, (unsigned int)n*sizeof(double)); }
  inline void copy_vector(const int n, const std::complex<double> *x, std::complex<double> *y) { memcpy(y, x, (unsigned int)n*sizeof(std::complex<double>)); }
#endif

  template<class T> inline
  void copy_vector(const int n, const T *x, T *y)
  {
    for (int i=0; i<n; i++)
      y[i] = x[i];
  }




  /*
    Copy vector x to vector y. Both vectors are of size n
    vector x elements are stored linearly with element increament incx
    vector y elements are stored linearly with element increament incx
  */
#if defined (HAVE_BLAS)
  inline void copy_vector(const int n, const double *x, const int incx,
			  double *y, const int incy)
  {
    blas::dcopy_(&n, x, &incx, y, &incy);
  }
  inline void copy_vector(const int n, const std::complex<double> *x,
			  const int incx, std::complex<double> *y,
			  const int incy)
  {
    blas::zcopy_(&n, x, &incx, y, &incy);
  }
#endif

  template<class T> inline
  void copy_vector(const int n, const T *x, const int incx, T *y, const int incy)
  {
    for (int i=0;i<n; i++)
      y[i*incy] = x[i*incx];
  }


  /*
    Swap vector x and vector y. Both vectors are of size n
  */
  inline void swap_vector(const int n, int *x, int *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }
  inline void swap_vector(const int n, short *x, short *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }
  inline void swap_vector(const int n, bin *x, bin *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }
  inline void swap_vector(const int n, float *x, float *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }
  inline void swap_vector(const int n, std::complex<float> *x, std::complex<float> *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }

#if defined (HAVE_BLAS)
  inline void swap_vector(const int n, double *x, double *y)
  {
    int incr = 1;
    blas::dswap_(&n, x, &incr, y, &incr);
  }
  inline void swap_vector(const int n, std::complex<double> *x,
			  std::complex<double> *y)
  {
    int incr = 1;
    blas::zswap_(&n, x, &incr, y, &incr);
  }
#else
  inline void swap_vector(const int n, double *x, double *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }
  inline void swap_vector(const int n, std::complex<double> *x, std::complex<double> *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }
#endif

  template<class T> inline
  void swap_vector(const int n, T *x, T *y)
  {
    T tmp;
    for (int i=0; i<n; i++) {
      tmp = y[i];
      y[i] = x[i];
      x[i] = tmp;
    }
  }


  /*
    Swap vector x and vector y. Both vectors are of size n
    vector x elements are stored linearly with element increament incx
    vector y elements are stored linearly with element increament incx
  */
  inline void swap_vector(const int n, int *x, const int incx, int *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }
  inline void swap_vector(const int n, short *x, const int incx, short *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }
  inline void swap_vector(const int n, bin *x, const int incx, bin *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }
  inline void swap_vector(const int n, float *x, const int incx, float *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }
  inline void swap_vector(const int n, std::complex<float> *x, const int incx, std::complex<float> *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }

#if defined (HAVE_BLAS)
  inline void swap_vector(const int n, double *x, const int incx, double *y,
			  const int incy)
  {
    blas::dswap_(&n, x, &incx, y, &incy);
  }
  inline void swap_vector(const int n, std::complex<double> *x, const int incx,
			  std::complex<double> *y, const int incy)
  {
    blas::zswap_(&n, x, &incx, y, &incy);
  }
#else
  inline void swap_vector(const int n, double *x, const int incx, double *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }
  inline void swap_vector(const int n, std::complex<double> *x, const int incx, std::complex<double> *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }
#endif

  template<class T> inline
  void swap_vector(const int n, T *x, const int incx, T *y, const int incy)
  {
    T tmp;
    for (int i=0; i<n; i++) {
      tmp = y[i*incy];
      y[i*incy] = x[i*incx];
      x[i*incx] = tmp;
    }
  }


  /*
   * Realise scaling operation: x = alpha*x
   */
#if defined(HAVE_BLAS)
  inline void scal_vector(int n, double alpha, double *x)
  {
    int incr = 1;
    blas::dscal_(&n, &alpha, x, &incr);
  }
  inline void scal_vector(int n, std::complex<double> alpha,
			  std::complex<double> *x)
  {
    int incr = 1;
    blas::zscal_(&n, &alpha, x, &incr);
  }
#endif

  template<typename T> inline
  void scal_vector(int n, T alpha, T *x)
  {
    if (alpha != T(1)) {
      for (int i = 0; i < n; ++i) {
	x[i] *= alpha;
      }
    }
  }


  /*
   * Realise scaling operation: x = alpha*x
   * Elements of x are stored linearly with increament incx
   */
#if defined(HAVE_BLAS)
  inline void scal_vector(int n, double alpha, double *x, int incx)
  {
    blas::dscal_(&n, &alpha, x, &incx);
  }
  inline void scal_vector(int n, std::complex<double> alpha,
			  std::complex<double> *x, int incx)
  {
    blas::zscal_(&n, &alpha, x, &incx);
  }
#endif

  template<typename T> inline
  void scal_vector(int n, T alpha, T *x, int incx)
  {
    if (alpha != T(1)) {
      for (int i = 0; i < n; ++i) {
	x[i*incx] *= alpha;
      }
    }
  }


  /*
   * Realise the following equation on vectors: y = alpha*x + y
   */
#if defined(HAVE_BLAS)
  inline void axpy_vector(int n, double alpha, const double *x, double *y)
  {
    int incr = 1;
    blas::daxpy_(&n, &alpha, x, &incr, y, &incr);
  }
  inline void axpy_vector(int n, std::complex<double> alpha,
			  const std::complex<double> *x,
			  std::complex<double> *y)
  {
    int incr = 1;
    blas::zaxpy_(&n, &alpha, x, &incr, y, &incr);
  }
#endif

  template<typename T> inline
  void axpy_vector(int n, T alpha, const T *x, T *y)
  {
    if (alpha != T(1)) {
      for (int i = 0; i < n; ++i) {
	y[i] += alpha * x[i];
      }
    }
    else {
      for (int i = 0; i < n; ++i) {
	y[i] += x[i];
      }
    }
  }


  /*
   * Realise the following equation on vectors: y = alpha*x + y
   * Elements of x are stored linearly with increment incx
   * and elements of y are stored linearly with increment incx
   */
#if defined(HAVE_BLAS)
  inline void axpy_vector(int n, double alpha, const double *x, int incx,
			  double *y, int incy)
  {
    blas::daxpy_(&n, &alpha, x, &incx, y, &incy);
  }
  inline void axpy_vector(int n, std::complex<double> alpha,
			  const std::complex<double> *x, int incx,
			  std::complex<double> *y, int incy)
  {
    blas::zaxpy_(&n, &alpha, x, &incx, y, &incy);
  }
#endif

  template<typename T> inline
  void axpy_vector(int n, T alpha, const T *x, int incx, T *y, int incy)
  {
    if (alpha != T(1)) {
      for (int i = 0; i < n; ++i) {
	y[i*incy] += alpha * x[i*incx];
      }
    }
    else {
      for (int i = 0; i < n; ++i) {
	y[i*incy] += x[i*incx];
      }
    }
  }


} // namespace itpp

//! \endcond

#endif // #ifndef COPY_VECTOR_H
