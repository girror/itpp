/*!
 * \file
 * \brief Trigonometric and hyperbolic functions - header file
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2009  (see AUTHORS file for a list of contributors)
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

#ifndef TRIG_HYP_H
#define TRIG_HYP_H

#include <itpp/base/itcompat.h>
#include <itpp/base/help_functions.h>


namespace itpp
{

//!\addtogroup trifunc
//!@{

//! Sinc function: sinc(x) = sin(pi*x)/pi*x
inline double sinc(double x)
{
  if (x == 0) {
    return 1.0;
  }
  else {
    double pix = itpp::pi * x;
    return sin(pix) / pix;
  }
}

//! Sine function
inline vec sin(const vec &x) { return apply_function<double>(std::sin, x); }
//! Sine function
inline mat sin(const mat &x) { return apply_function<double>(std::sin, x); }
//! Cosine function
inline vec cos(const vec &x) { return apply_function<double>(std::cos, x); }
//! Cosine function
inline mat cos(const mat &x) { return apply_function<double>(std::cos, x); }
//! Tan function
inline vec tan(const vec &x) { return apply_function<double>(std::tan, x); }
//! Tan function
inline mat tan(const mat &x) { return apply_function<double>(std::tan, x); }
//! Inverse sine function
inline vec asin(const vec &x) { return apply_function<double>(std::asin, x); }
//! Inverse sine function
inline mat asin(const mat &x) { return apply_function<double>(std::asin, x); }
//! Inverse cosine function
inline vec acos(const vec &x) { return apply_function<double>(std::acos, x); }
//! Inverse cosine function
inline mat acos(const mat &x) { return apply_function<double>(std::acos, x); }
//! Inverse tan function
inline vec atan(const vec &x) { return apply_function<double>(std::atan, x); }
//! Inverse tan function
inline mat atan(const mat &x) { return apply_function<double>(std::atan, x); }
//! Sinc function, sin(pi*x)/(pi*x)
inline vec sinc(const vec &x) { return apply_function<double>(sinc, x); }
//! Sinc function, sin(pi*x)/(pi*x)
inline mat sinc(const mat &x) { return apply_function<double>(sinc, x); }

//!@}


//!\addtogroup hypfunc
//!@{

//! Sine hyperbolic function
inline vec sinh(const vec &x) { return apply_function<double>(std::sinh, x); }
//! Sine hyperbolic function
inline mat sinh(const mat &x) { return apply_function<double>(std::sinh, x); }
//! Cosine hyperbolic function
inline vec cosh(const vec &x) { return apply_function<double>(std::cosh, x); }
//! Cosine hyperbolic function
inline mat cosh(const mat &x) { return apply_function<double>(std::cosh, x); }
//! Tan hyperbolic function
inline vec tanh(const vec &x) { return apply_function<double>(std::tanh, x); }
//! Tan hyperbolic function
inline mat tanh(const mat &x) { return apply_function<double>(std::tanh, x); }
//! Inverse sine hyperbolic function
inline vec asinh(const vec &x) { return apply_function<double>(::asinh, x); }
//! Inverse sine hyperbolic function
inline mat asinh(const mat &x) { return apply_function<double>(::asinh, x); }
//! Inverse cosine hyperbolic function
inline vec acosh(const vec &x) { return apply_function<double>(::acosh, x); }
//! Inverse cosine hyperbolic function
inline mat acosh(const mat &x) { return apply_function<double>(::acosh, x); }
//! Inverse tan hyperbolic function
inline vec atanh(const vec &x) { return apply_function<double>(::atanh, x); }
//! Inverse tan hyperbolic function
inline mat atanh(const mat &x) { return apply_function<double>(::atanh, x); }

//!@}

} // namespace itpp

#endif // #ifndef TRIG_HYP_H
