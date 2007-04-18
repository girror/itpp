/*!
 * \file 
 * \brief Implementation of Bessel functions
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

#include <itpp/base/bessel.h>
#include <itpp/base/bessel/bessel_internal.h>
#include <cmath>


namespace itpp { 

  // Bessel function of order nu
  double besselj(int nu, double x) 
  { 
#ifdef _MSC_VER
    return _jn(nu, x); 
#else
    return jn(nu, x); 
#endif
  }

  vec besselj(int nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
#ifdef _MSC_VER
      out(i) = _jn(nu, x(i));
#else
      out(i) = jn(nu, x(i));
#endif

    return out;
  }

  // Bessel function of order nu. nu is real.
  double besselj(double nu, double x) { return jv(nu, x); }

  vec besselj(double nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
      out(i) = jv(nu, x(i));

    return out;
  }

  // Bessel function of second kind of order nu
  double bessely(int nu, double x) 
  { 
#ifdef _MSC_VER
    return _yn(nu, x); 
#else
    return yn(nu, x); 
#endif
  }

  vec bessely(int nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
#ifdef _MSC_VER
      out(i) = _yn(nu, x(i));
#else
      out(i) = yn(nu, x(i));
#endif

    return out;
  }
  // Bessel function of second kind of order nu
  double bessely(double nu, double x) { return yv(nu, x); }

  vec bessely(double nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
      out(i) = yv(nu, x(i));

    return out;
  }

  // Modified Bessel function of order nu
  double besseli(double nu, double x) { return iv(nu, x); }

  vec besseli(double nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
      out(i) = iv(nu, x(i));

    return out;
  }

  // Modified Bessel function of second kind of order n
  double besselk(int n, double x) { return kn(n, x); }

  vec besselk(int nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
      out(i) = kn(nu, x(i));

    return out;
  }

} // namespace itpp
