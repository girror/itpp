/*!
 * \file
 * \brief Definition of numerical integration
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
 *
 * -------------------------------------------------------------------------
 */

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <itpp/itconfig.h>
#include <itpp/base/vec.h>
#include <limits>


namespace itpp {

  /*! 
    \addtogroup integration
    \brief Numerical integration routines
  */
  
  //@{

  /*!
    1-dimensional numerical Simpson quadrature integration

    Calculate the 1-dimensional integral
    \f[
    \int_a^b f(x) dx
    \f]

    Uses an adaptive Simpson quadrature method. See [Gander] for more
    details. The integrand is specified as a function \code double
    f(double) \endcode.

    Example:
    \code
    #include "itpp/itbase.h"

    double f(const double x)
    {
      return x*log(x);
    }

    int main()
    {
      double res = quad( f, 1.5, 3.5);
      cout << "res = " << res << endl;

      return 0;
    }
    \endcode

    References:
    
    [Gander] Gander, W. and W. Gautschi, "Adaptive Quadrature -
    Revisited", BIT, Vol. 40, 2000, pp. 84-101.  
		This document is also available at http://www.inf.ethz.ch/personal/gander.
  */
  double quad(double (*f)(double), const double a, const double b, const double tol=std::numeric_limits<double>::epsilon());

  /*!
    1-dimensional numerical adaptive Lobatto quadrature integration
    
    Calculate the 1-dimensional integral
    \f[
    \int_a^b f(x) dx
    \f]

    Uses an adaptive Lobatto quadrature method. See [Gander] for more
    details. The integrand is specified as a function \code double
    f(double) \endcode.

    Example:
    \code
    #include "itpp/itbase.h"

    double f(const double x)
    {
      return x*log(x);
    }

    int main()
    {
      double res = quadl( f, 1.5, 3.5);
      cout << "res = " << res << endl;

      return 0;
    }
    \endcode

    References:
    
    [Gander] Gander, W. and W. Gautschi, "Adaptive Quadrature -
    Revisited", BIT, Vol. 40, 2000, pp. 84-101.
  	This document is also available at http:// www.inf.ethz.ch/personal/gander.
  */
  double quadl(double (*f)(double), const double a, const double b, const double tol=std::numeric_limits<double>::epsilon());

  //@}
  
} // namespace itpp

#endif // #ifndef INTEGRATION_H
