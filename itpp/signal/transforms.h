/*!
 * \file
 * \brief Fourier, Hadamard, Walsh-Hadamard, and 2D Hadamard transforms -
 *        header file
 * \author Tony Ottosson, Thomas Eriksson, Simon Wood and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/matfunc.h>


namespace itpp
{

/*!
  \addtogroup fft
  \brief One dimensional fast fourier transform
  \author Tony Ottosson and Adam Piatyszek

  The functions \code X = fft(x) \endcode and \code x = ifft(X) \endcode are
  the fourier and inverse fourier transforms of size \a N defined as:
  \f[
  X(k) = \sum_{j=0}^{N-1} x(j) e^{-2\pi j k \cdot i / N}
  \f]

  \f[
  x(j) = \frac{1}{N} \sum_{k=0}^{N-1} X(k) e^{2\pi j k \cdot i / N}
  \f]

  \code Y = fft(X, N) \endcode performs zero-padding up to size N and then
  performs an N-size fft.

  The implementation is built upon one of the following libraries:
  - FFTW (version 3.0.0 or higher)
  - MKL (version 8.0.0 or higher)
  - ACML (version 2.5.3 or higher).

  \note FFTW-based implementation is the fastest for powers of two.
  Furthermore, the second time you call the routine with the same size,
  the calculation is much faster due to many things were calculated and
  stored the first time the routine was called.

  \note Achieving maximum runtime efficiency with the FFTW library on some
  computer architectures requires that data are stored in the memory with
  a special alignment (to 16-byte boundaries). The IT++ memory management
  functions and container classes do not generally allocate memory aligned
  this way, and as a result calling FFTW via the IT++ interface (i.e. the
  fft() function) may be slower than using the FFTW library directly.
  Therefore, FFTW users concerned about maximum possible performance may
  want to consider the possibility of calling the FFTW library and its
  memory management/allocation routines directly, bypassing the IT++
  storage classes and the fft() interface to FFTW.
*/

//!\addtogroup fft
//!@{

//! Fast Fourier Transform
void fft(const cvec &in, cvec &out);
//! Fast Fourier Transform
cvec fft(const cvec &in);
//! Fast Fourier Transform, with zero-padding up to size N
cvec fft(const cvec &in, const int N);
//! Inverse Fast Fourier Transform
void ifft(const cvec &in, cvec &out);
//! Inverse Fast Fourier Transform
cvec ifft(const cvec &in);
//! Inverse Fast Fourier Transform, with zero-padding up to size N
cvec ifft(const cvec &in, const int N);

//! Real Fast Fourier Transform
void fft_real(const vec& in, cvec &out);
//! Real Fast Fourier Transform
cvec fft_real(const vec& in);
//! Real Fast Fourier Transform, with zero-padding up to size N
cvec fft_real(const vec &in, const int N);
//! Inverse Real Fast Fourier Transform. Assumes even size.
void ifft_real(const cvec &in, vec &out);
//! Inverse Real Fast Fourier Transform. Assumes even size.
vec ifft_real(const cvec &in);
//! Inverse Real Fast Fourier Transform, with zero-padding up to size N
vec ifft_real(const cvec &in, const int N);
//!@}


/*!
  \addtogroup dct
  \brief One dimensional Dicrete Cosine Transform
  \author Tony Ottosson and Adam Piatyszek

  The functions \code X = dct(x) \endcode and \code x = idct(X) \endcode
  are the dicrete cosine and inverse discrete cosine transforms of size \a
  N defined as:
  \f[
  X(k) = w(k) \sum_{j=0}^{N-1} x(j) \cos \left(\frac{(2j+1)k \pi}{2N} \right)
  \f]

  \f[
  x(j) = \sum_{k=0}^{N-1} w(k) X(k) \cos \left(\frac{(2j+1)k \pi}{2N} \right)
  \f]
  where \f$w(k) = 1/sqrt{N}\f$ for \f$k=0\f$ and
  \f$w(k) = sqrt{2/N}\f$ for \f$k\geq 1\f$.

  The implementation is built upon one of the following libraries:
  - FFTW (version 3.0.0 or higher)
  - MKL (version 8.0.0 or higher)
  - ACML (version 2.5.3 or higher).

  \note FFTW-based implementation is the fastest for powers of two.
  Furthermore, the second time you call the routine with the same size,
  the calculation is much faster due to many things were calculated and
  stored the first time the routine was called.

  \note Achieving maximum runtime efficiency with the FFTW library on some
  computer architectures requires that data are stored in the memory with
  a special alignment (to 16-byte boundaries). The IT++ memory management
  functions and container classes do not generally allocate memory aligned
  this way, and as a result calling FFTW via the IT++ interface (i.e. the
  dct()/idct() function) may be slower than using the FFTW library
  directly. Therefore, FFTW users concerned about maximum possible
  performance may want to consider the possibility of calling the FFTW
  library and its memory management/allocation routines directly,
  bypassing the IT++ storage classes and the dct()/idct() interface to
  FFTW.
*/

//!\addtogroup dct
//!@{

//! Discrete Cosine Transform (DCT)
void dct(const vec &in, vec &out);
//! Discrete Cosine Transform (DCT)
vec dct(const vec &in);
//! Inverse Discrete Cosine Transform (IDCT)
void idct(const vec &in, vec &out);
//! Inverse Discrete Cosine Transform (IDCT)
vec idct(const vec &in);
//!@}


//!\addtogroup fht
//!@{

//! Fast Hadamard Transform
template <class T> Vec<T> dht(const Vec<T> &v);
//! Fast Hadamard Transform
template <class T> void dht(const Vec<T> &vin, Vec<T> &vout);
//! Fast Hadamard Transform - memory efficient. Stores the result in \c v
template <class T> void self_dht(Vec<T> &v);

//! Fast Walsh Hadamard Transform
template <class T> Vec<T> dwht(const Vec<T> &v);
//! Fast Walsh Hadamard Transform
template <class T> void dwht(const Vec<T> &vin, Vec<T> &vout);
//! Fast Walsh Hadamard Transform - memory efficient (result in \c v)
template <class T> void self_dwht(Vec<T> &v);

//! Fast 2D Hadamard Transform
template <class T> Mat<T> dht2(const Mat<T> &m);
//! Fast 2D Walsh Hadamard Transform
template <class T> Mat<T> dwht2(const Mat<T> &m);
//!@}

template <class T>
Vec<T> dht(const Vec<T> &v)
{
  Vec<T> ret(v.size());
  dht(v, ret);
  return ret;
}

//! Bit reverse
template <class T>
void bitrv(Vec<T> &out)
{
  int N = out.size();
  int j = 0;
  int N1 = N - 1;
  for (int i = 0; i < N1; ++i) {
    if (i < j) {
      T temp = out[j];
      out[j] = out[i];
      out[i] = temp;
    }
    int K = N / 2;
    while (K <= j) {
      j -= K;
      K /= 2;
    }
    j += K;
  }
}

template <class T>
void dht(const Vec<T> &vin, Vec<T> &vout)
{
  int N = vin.size();
  int m = levels2bits(N);
  it_assert_debug((1 << m) == N, "dht(): The vector size must be a power of two");

  vout.set_size(N);

  // This step is separated because it copies vin to vout
  for (int ib = 0; ib < N; ib += 2) {
    vout(ib) = vin(ib) + vin(ib + 1);
    vout(ib + 1) = vin(ib) - vin(ib + 1);
  }
  N /= 2;

  int l = 2;
  for (int i = 1; i < m; ++i) {
    N /= 2;
    int ib = 0;
    for (int k = 0; k < N; ++k) {
      for (int j = 0; j < l; ++j) {
        T t = vout(ib + j);
        vout(ib + j) += vout(ib + j + l);
        vout(ib + j + l) = t - vout(ib + j + l);
      }
      ib += 2 * l;
    }
    l *= 2;
  }

  vout /= static_cast<T>(std::sqrt(static_cast<double>(vin.size())));
}

template <class T>
void self_dht(Vec<T> &v)
{
  int N = v.size();
  int m = levels2bits(N);
  it_assert_debug((1 << m) == N, "self_dht(): The vector size must be a power "
                  "of two");

  int l = 1;
  for (int i = 0; i < m; ++i) {
    N /= 2;
    int ib = 0;
    for (int k = 0; k < N; ++k) {
      for (int j = 0; j < l; ++j) {
        T t = v(ib + j);
        v(ib + j) += v(ib + j + l);
        v(ib + j + l) = t - v(ib + j + l);
      }
      ib += 2 * l;
    }
    l *= 2;
  }

  v /= static_cast<T>(std::sqrt(static_cast<double>(v.size())));
}

template <class T>
Vec<T> dwht(const Vec<T> &v)
{
  Vec<T> ret(v.size());
  dwht(v, ret);
  return ret;
}

template <class T>
void dwht(const Vec<T> &vin, Vec<T> &vout)
{
  dht(vin, vout);
  bitrv(vout);
}


template <class T>
void self_dwht(Vec<T> &v)
{
  self_dht(v);
  bitrv(v);
}

template <class T>
Mat<T> dht2(const Mat<T> &m)
{
  Mat<T> ret(m.rows(), m.cols());
  Vec<T> v;

  for (int i = 0; i < m.rows(); ++i) {
    v = m.get_row(i);
    self_dht(v);
    ret.set_row(i, v);
  }
  for (int i = 0; i < m.cols(); ++i) {
    v = ret.get_col(i);
    self_dht(v);
    ret.set_col(i, v);
  }

  return transpose(ret);
}

template <class T>
Mat<T> dwht2(const Mat<T> &m)
{
  Mat<T> ret(m.rows(), m.cols());
  Vec<T> v;

  for (int i = 0; i < m.rows(); ++i) {
    v = m.get_row(i);
    self_dwht(v);
    ret.set_row(i, v);
  }
  for (int i = 0; i < m.cols(); ++i) {
    v = ret.get_col(i);
    self_dwht(v);
    ret.set_col(i, v);
  }

  return transpose(ret);
}

//! \cond

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

#ifdef HAVE_EXTERN_TEMPLATE

extern template vec dht(const vec &v);
extern template cvec dht(const cvec &v);
extern template void dht(const vec &vin, vec &vout);
extern template void dht(const cvec &vin, cvec &vout);

extern template void self_dht(vec &v);
extern template void self_dht(cvec &v);

extern template vec dwht(const vec &v);
extern template cvec dwht(const cvec &v);
extern template void dwht(const vec &vin, vec &vout);
extern template void dwht(const cvec &vin, cvec &vout);

extern template void self_dwht(vec &v);
extern template void self_dwht(cvec &v);

extern template mat  dht2(const mat &m);
extern template cmat dht2(const cmat &m);

extern template mat  dwht2(const mat &m);
extern template cmat dwht2(const cmat &m);

#endif // HAVE_EXTERN_TEMPLATE

//! \endcond

} // namespace itpp

#endif // #ifndef TRANSFORMS_H
