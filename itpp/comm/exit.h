/*!
 * \file
 * \brief Definitions for EXtrinsic Information Transfer (EXIT) chart class
 * \author Bogdan Cristea
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

#ifndef EXIT_H
#define EXIT_H

#include <itpp/itbase.h>
#include <itpp/comm/modulator.h> //BPSK class for a priori information generation
#include <itpp/stat/histogram.h> //histogram class for mutual information computation
//#include <itpp/base/itcompat.h>

namespace itpp
{

/*!
  \ingroup misccommfunc
  \brief EXtrinsic Information Transfer (%EXIT) chart

  Computes the A priori Mutual Information assuming a Gaussian distribution of the a priori information and
  the Extrinsic Mutual Information between the emitted bits and their extrinsic information

  Description:
  - the a priori mutual information is computed using relation (14)
  - the extrinsic mutual information is computed by estimating first the conditional Probability Density Functions (PDF),
  given the emitted bits, and then numerically integrating according to relation (19)

  Reference:
  Stephan ten Brink, ''Convergence behavior of iteratively decoded parallel concatenated codes,``
  IEEE Transactions on Communications, oct. 2001
*/
class EXIT
{
public:
    //! Computes the a priori mutual information
    /*! It is assumed that the a priori information has a Gaussian distribution
     */
    double apriori_mutual_info(const double &in_sigma2A, //!< variance of the a priori information
                               const double &lim=100 //!< [-lim,+lim] is the integration interval (theoretically it should be [-inf,+inf])
                              )
    {
        sigma2A = in_sigma2A;
        return double(1)-itpp::quad(&gaussian_fct, -lim, lim);
    };
    //! Generates a priori information assuming a Gaussian distribution of the a priori information
    /*! The variance of the a priori information must be already initialized through EXIT::apriori_mutual_info function.
     * The information generated in this way is used sometimes as intrinsic information at the SISO module input.
     */
    itpp::vec generate_apriori_info(const itpp::bvec &bits)
    {
        itpp::BPSK bpsk;
        return (-sigma2A/2)*bpsk.modulate_bits(bits)+std::sqrt(sigma2A)*itpp::randn(bits.length());
    };
    //! Computes the extrinsic mutual information
    /*! The conditional Probability Density Function (PDF) of the extrinsic information is estimated using the histogram of the
     * extrinsic information and the knowledge of the emitted bits corresponding to the extrinsic information.
     */
    double extrinsic_mutual_info(const itpp::vec &obs, //!< extrinsic information obtained from the SISO module output
                                 const itpp::bvec &cond, //!< emitted bits corresponding to the extrinsic information
                                 const int &N=100 //!< number of subintervals used to compute the histogram
                                );
private:
    static double sigma2A;
    friend double itpp::quad(double (*f)(double), double a, double b, double tol);
    static double gaussian_fct(const double x)
    {
        return (1.0/std::sqrt(sigma2A*itpp::m_2pi))*std::exp(-itpp::sqr(x-(sigma2A/2.0))/(2.0*sigma2A))*::log2(1+std::exp(-x));
    };
};

}
#endif /* EXIT_H_ */
