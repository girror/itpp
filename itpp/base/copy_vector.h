/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2005 by Tony Ottosson, Thomas Eriksson, P�l Frenger,   *
 * Tobias Ringstr�m, and Jonas Samuelsson.                                   *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*!
  \file
  \brief Vector copy functions for internal use
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __copy_vector_h
#define __copy_vector_h

#include "itpp/itconfig.h"
#include "itpp/base/binary.h"

#if defined (HAVE_CBLAS) || defined(HAVE_MKL)
#include "itpp/base/cblas.h"
#endif

//using std::complex;
//using std::swap;

namespace itpp {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  /*
    Copy vector x to vector y. Both vectors are of size n
  */
  inline void copy_vector(const int n, const int *x, int *y) { memcpy(y, x, (unsigned int)n*sizeof(int)); }
  inline void copy_vector(const int n, const short *x, short *y) { memcpy(y, x, (unsigned int)n*sizeof(short)); }
  inline void copy_vector(const int n, const bin *x, bin *y) { memcpy(y, x, (unsigned int)n*sizeof(bin)); }
  inline void copy_vector(const int n, const float *x, float *y) { memcpy(y, x, (unsigned int)n*sizeof(float)); }
  inline void copy_vector(const int n, const std::complex<float> *x, std::complex<float> *y) { memcpy(y, x, (unsigned int)n*sizeof(std::complex<float>)); }

#if defined (HAVE_CBLAS) || defined(HAVE_MKL)
  inline void copy_vector(const int n, const double *x, double *y) { cblas_dcopy(n, x, 1, y, 1); }
  inline void copy_vector(const int n, const std::complex<double> *x, std::complex<double> *y) { cblas_zcopy(n, x, 1, y, 1); }
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
#if defined (HAVE_CBLAS) || defined(HAVE_MKL)
  inline void copy_vector(const int n, const double *x, const int incx, double *y, const int incy) { cblas_dcopy(n, x, incx, y, incy); }
  inline void copy_vector(const int n, const std::complex<double> *x, const int incx, std::complex<double> *y, const int incy) { cblas_zcopy(n, x, incx, y, incy); }
#endif

  template<class T> inline
  void copy_vector(const int n, const T *x, const int incx, T *y, const int incy)
  {
    for (int i=0;i<n; i++)
      y[i*incy] = x[i*incx];
  }


  /*
    Swap vector x to vector y. Both vectors are of size n
  */
  inline void swap_vector(const int n, int *x, int *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }
  inline void swap_vector(const int n, short *x, short *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }
  inline void swap_vector(const int n, bin *x, bin *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }
  inline void swap_vector(const int n, float *x, float *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }
  inline void swap_vector(const int n, std::complex<float> *x, std::complex<float> *y) { for (int i=0; i<n; i++) std::swap(x[i], y[i]); }

#if defined (HAVE_CBLAS) || defined(HAVE_MKL)
  inline void swap_vector(const int n, double *x, double *y) { cblas_dswap(n, x, 1, y, 1); }
  inline void swap_vector(const int n, std::complex<double> *x, std::complex<double> *y) { cblas_zswap(n, x, 1, y, 1); }
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
    Swap vector x to vector y. Both vectors are of size n
    vector x elements are stored linearly with element increament incx
    vector y elements are stored linearly with element increament incx
  */
  inline void swap_vector(const int n, int *x, const int incx, int *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }
  inline void swap_vector(const int n, short *x, const int incx, short *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }
  inline void swap_vector(const int n, bin *x, const int incx, bin *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }
  inline void swap_vector(const int n, float *x, const int incx, float *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }
  inline void swap_vector(const int n, std::complex<float> *x, const int incx, std::complex<float> *y, const int incy) { for (int i=0; i<n; i++) std::swap(x[i*incx], y[i*incy]); }

#if defined (HAVE_CBLAS) || defined(HAVE_MKL)
  inline void swap_vector(const int n, double *x, const int incx, double *y, const int incy) { cblas_dswap(n, x, incx, y, incy); }
  inline void swap_vector(const int n, std::complex<double> *x, const int incx, std::complex<double> *y, const int incy) { cblas_zswap(n, x, incx, y, incy); }
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


#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} //namespace itpp

#endif // __copy_vector_h
