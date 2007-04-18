/*!
 * \file
 * \brief Definitions of scalar functions
 * \author Tony Ottosson, Pal Frenger and Adam Piatyszek
 * 
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#ifndef SCALFUNC_H
#define SCALFUNC_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <itpp/base/itassert.h>
#include <complex>
#include <limits>


  /*! \addtogroup miscfunc */
  //!@{

#ifndef HAVE_ASINH
  //! Arcus sinhyp
  double asinh(double x);
#endif

#ifndef HAVE_ACOSH
  //! Arcus coshyp
  double acosh(double x);
#endif

#ifndef HAVE_ATANH
  //! Arcus tanhyp
  double atanh(double x);
#endif
  //!@}


#ifdef _MSC_VER
// These functions are part of C99. But apparently not in Visual C++.
/*!
  \brief Error function
  \ingroup errorfunc
*/
double erf(double x);
/*!
  \brief Complementary error function
  \ingroup errorfunc
*/
double erfc(double x);
#endif

#ifndef HAVE_TGAMMA
/*!
  \brief True gamma function
  \ingroup miscfunc
*/
double tgamma(double x);
#endif

#if !defined(HAVE_LGAMMA) || (HAVE_DECL_SIGNGAM != 1)
/*!
  \brief Lograrithm of an absolute gamma function
  \ingroup miscfunc
*/
double lgamma(double x);
extern int signgam;
#endif

#ifndef HAVE_CBRT
/*
  \brief Cubic root
  \ingroup miscfunc
*/
double cbrt(double x);
#endif


// Fix log2 for some platforms, that have it defined as a macro
#if defined (log2)
#undef log2
#endif

namespace itpp {

  //! Constant pi
  const double pi = 3.14159265358979323846;

  //! Constant 2*Pi
  const double m_2pi = 2 * pi;

  //! Constant eps
  const double eps = std::numeric_limits<double>::epsilon();

  //! Constant definition to speed up trunc_log and trunc_exp functions
  const double log_double_max = std::log(std::numeric_limits<double>::max());
  //! Constant definition to speed up trunc_log and trunc_exp functions
  const double log_double_min = std::log(std::numeric_limits<double>::min());

  /*!
    \brief Truncated natural logarithm function

    This truncated function provides a solution in the cases when the
    logarithm argument is less or equal to zero or infinity. The function
    checks for such extreme values and use some kind of truncation
    (saturation) before calculating the logarithm.

    The truncated logarithm function can be used for calculation of
    log-likelihood in soft demodulators, when numerical instability problem
    might occur.
  */
  inline double trunc_log(double x)
  {
    if (std::numeric_limits<double>::is_iec559) {
      if (x == std::numeric_limits<double>::infinity())
	return log_double_max;
      if (x <= 0)
	return log_double_min;
    }
    return std::log(x);
  }

  /*!
    \brief Truncated exponential function

    This truncated function provides a solution in the case when the
    exponent function results in infinity. The function checks for an
    extreme value and use truncation (saturation) before calculating 
    the result.

    The truncated exponential function can be used  when numerical
    instability problem occurs for a standard exp function.
  */
  inline double trunc_exp(double x)
  {
    if (std::numeric_limits<double>::is_iec559
	&& (x >= log_double_max))
      return std::numeric_limits<double>::max();
    return std::exp(x);
  }

  /*!
    \brief Gamma function
    \ingroup miscfunc
  */
  double gamma(double x);

  /*!
    \brief Q-function
    \ingroup errorfunc
  */
  double Qfunc(double x);

  /*!
    \brief Inverse of error function
    \ingroup errorfunc
  */
  double erfinv(double x);

  /*!
   * \brief Error function for complex argument
	 * \ingroup errorfunc
	 * \author Adam Piatyszek
   *
   * This function calculates a well known error function \c erf(z)
   * for complex \c z. The implementation is based on unofficial
   * implementation for Octave. Here is a part of the author's note
   * from original sources:
	 *
	 * Put together by John Smith john at arrows dot demon dot co dot uk, 
	 * using ideas by others.
	 *
	 * Calculate \c erf(z) for complex \c z.
	 * Three methods are implemented; which one is used depends on z.
	 *
	 * The code includes some hard coded constants that are intended to
	 * give about 14 decimal places of accuracy. This is appropriate for
	 * 64-bit floating point numbers. 
	 */
  std::complex<double> erf(const std::complex<double>& z);

  /*!
    \brief Base-2 logarithm
    \ingroup logexpfunc
  */
  inline double log2(double x) { return (std::log(x)/0.693147180559945309417); }

  /*!
    \brief Base-b logarithm
    \ingroup logexpfunc
  */
  inline double logb(double b, double x) { return std::log(x)/std::log(b); }

  /*!
    \brief Sinc function. sinc(x) = sin(pi*x)/pi*x
    \ingroup miscfunc
  */
  inline double sinc(double x) { return ( (x==0) ? 1 : (sin(itpp::pi*x)/itpp::pi/x) ); }

  /*! \addtogroup miscfunc */
  //!@{

#ifdef _MSC_VER
  //! Round to nearest integer
  inline double rint(double x) { return floor(x+0.5); }
#endif

  //! Round to nearest integer
  inline int round_i(double x) { return int(rint(x)); }
  //! The nearest larger integer
  inline int ceil_i(double x) { return int(ceil(x)); }
  //! The nearest smaller integer
  inline int floor_i(double x) { return int(floor(x)); }

  //! Round to nearest integer, return result in double
  inline double round(double x) { return rint(x); }

  //! Return true if x is an integer
  inline bool is_int(double x) { double dummy; return (modf(x, &dummy) == 0.0); }

  //! Return true if x is an even integer
  inline bool is_even(int x) { return ((x&1) == 0); }
  //!@}

  /*! \addtogroup logexpfunc */
  //!@{

  //! Calculate the number of bits needed to represent an inteager \c n
  inline int int2bits(int n)
  {
    it_assert(n >= 0, "int2bits(): Improper argument value");

    if (n == 0)
      return 1;

    int b = 0;
    while (n) { 
      n >>= 1; 
      ++b; 
    }
    return b;
  }
  
  //! Calculate the number of bits needed to represent \c n different values (levels).
  inline int levels2bits(int n)
  {
    it_assert(n > 0,"levels2bits(): Improper argument value");
    return int2bits(--n);
  }

  //! Deprecated function. Might be removed in future releases. Please use int2bits() or levels2bits() instead.
  inline int needed_bits(int n)
  {
    it_warning("needed_bits(): This function is depreceted. Might be removed in future releases. Depending on your needs, please use int2bits() or levels2bits() instead.");
    return int2bits(n);
  }

  //! Deprecated function. Might be removed in future releases. Please use int2bits() or levels2bits() instead.
  inline int needed_bits(double n) 
  { 
    it_warning("needed_bits(): This function is depreceted. Might be removed in future releases. Depending on your needs, please use int2bits() or levels2bits() instead.");
    it_assert(n > 0, "needed_bits(): Improper argument value"); 
    return ceil_i(log2(n));
  }

  //! Integer 2^x
  inline int pow2i(int x) { return ((x < 0) ? 0 : (1 << x)); }
  //! Calculate two to the power of x (2^x)
  inline int pow2(int x) { return pow2i(x); }

  //! Calculate two to the power of x (2^x)
  inline double pow2(double x) { return pow(2.0, x); }
  //! Calculate ten to the power of x (10^x)
  inline double pow10(double x) { return pow(10.0, x); }

  //! Decibel of x (10*log10(x))
  inline double dB(double x) { return 10.0 * log10(x); }
  //! Inverse of decibel
  inline double inv_dB(double x) { return pow(10.0, x/10.0); }
  //!@}

  /*! \addtogroup miscfunc */
  //!@{

  //! Convert to Gray Code
  inline int gray_code(int x) { return x^(x>>1); }

  //! Compute the binomial coefficient "n over k" as a float.
  double binom(int n, int k);

  //! Compute the binomial coefficient "n over k" as an integer.
  int binom_i(int n, int k);

  //! Compute the base 10 logarithm of the binomial coefficient "n over k".
  double log_binom(int n, int k);

  //! Convert radians to degrees
  inline double rad_to_deg(double x) { return 180.0 / itpp::pi * x; }
  //! Convert degrees to radians
  inline double deg_to_rad(double x) { return itpp::pi / 180.0 * x; }
  //!@}


  /*! \addtogroup miscfunc */
  //!@{

  //! Square of x
  inline double sqr(double x) { return x*x; }
  //! Square of complex-valued x, ||x||^2
  inline double sqr(std::complex<double> x) { return (x.real()*x.real()+x.imag()*x.imag()); }
  //! The reminder of the division x/y
  inline double rem(double x, double y) {return fmod(x,y);}
  //! The sign of x
  inline double sign(double x) { return x==0.0 ? 0.0 : (x<0.0 ? -1.0 : 1.0); }
  //! The sign of x
  inline double sgn(double x) { return x==0.0 ? 0.0 : (x<0.0 ? -1.0 : 1.0); }

  //! Absolute value
  inline signed char abs(signed char x) { return x>0 ? x : -x; }
  //! Absolute value
  inline short abs(short x) { return x>0 ? x : -x; }
  //! Absolute value
  inline int abs(int x) { return x>0 ? x : -x; }

  //double sigmoid(double x) { return 1.0/(1.0+exp(-x)); }

  //! Calculates factorial coefficient for index <= 170.
  double fact(int index);

  //! Calculates the modulus, i.e. the signed reminder after division
  long mod(long k, long n);

  /*!
    \brief returns the greatest common divisor (GCD) \a g of the elements \a a and \a b.

    \a a and \a b must be non-negative integers. \a gdc(0,0) is 0 by convention; all other
    GCDs are positive integers.
  */
  long gcd(long a, long b);

  //! Round \a x to zero if \a abs(x) is smaller than \a threshold
  inline double round_to_zero(double x, double threshold = 1e-14) {
    return (std::abs(x) < threshold) ? 0.0 : x;
  }
  
  //! Round each part of \a x smaller than \a threshold to zero  
  std::complex<double> round_to_zero(const std::complex<double>& x, 
				     double threshold = 1e-14);

  //!@}

} // namespace itpp

#endif // #ifndef SCALFUNC_H
