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
  \brief 
  \author Thomas Eriksson

  $Revision$

  $Date$
*/

#ifndef __gmm_h
#define __gmm_h

#include <string>
#include <itpp/itconfig.h>
#include <itpp/base/itassert.h>
#include <iostream>
#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/elmatfunc.h>
#include <itpp/base/array.h>
#include <cstring>

namespace itpp {


  /*! 
    \ingroup sourcecoding
    \brief Gaussian Mixture Model Class
    \author Thomas Eriksson
   */
  class GMM {
  public:
    GMM();
    GMM(int nomix, int dim);
    GMM(std::string filename);
    void init_from_vq(const vec &codebook, int dim);
    //	void init(const vec &w_in, const vec &m_in, const vec &sigma_in);
    void init(const vec &w_in, const mat &m_in, const mat &sigma_in);
    void load(std::string filename);
    void save(std::string filename);
    void set_weight(const vec &weights, bool compflag=true);
    void set_weight(int i, double weight, bool compflag=true);
    void set_mean(const mat &m_in);
    void set_mean(const vec &means, bool compflag=true);
    void set_mean(int i, const vec &means, bool compflag=true);
    void set_covariance(const mat &sigma_in);
    void set_covariance(const vec &covariances, bool compflag=true);
    void set_covariance(int i, const vec &covariances, bool compflag=true);
    int get_no_mixtures();
    int get_dimension();
    vec get_weight();
    double get_weight(int i);
    vec get_mean();
    vec get_mean(int i);
    vec get_covariance();
    vec get_covariance(int i);
    void marginalize(int d_new);
    void join(const GMM &newgmm);
    void clear();
    double likelihood(const vec &x);
    double likelihood_aposteriori(const vec &x, int mixture);
    vec likelihood_aposteriori(const vec &x);
    vec draw_sample();
  protected:
    vec			m,sigma,w;
    int			M,d;
  private:
    void		compute_internals();
    vec			normweight,normexp;
  };

  inline void GMM::set_weight(const vec &weights, bool compflag) {w=weights; if (compflag) compute_internals(); }
  inline void GMM::set_weight(int i, double weight, bool compflag) {w(i)=weight; if (compflag) compute_internals(); }
  inline void GMM::set_mean(const vec &means, bool compflag) {m=means; if (compflag) compute_internals(); }
  inline void GMM::set_covariance(const vec &covariances, bool compflag) {sigma=covariances; if (compflag) compute_internals(); }
  inline int GMM::get_no_mixtures() {return M;}
  inline int GMM::get_dimension() {return d;}
  inline vec GMM::get_weight() {return w;}
  inline double GMM::get_weight(int i) {return w(i);}
  inline vec GMM::get_mean() {return m;}
  inline vec GMM::get_mean(int i) {return m.mid(i*d,d);}
  inline vec GMM::get_covariance() {return sigma;}
  inline vec GMM::get_covariance(int i) {return sigma.mid(i*d,d);}

  GMM gmmtrain(Array<vec> &TrainingData, int M, int NOITER=30, bool VERBOSE=true);

} // namespace itpp

#endif // __gmm_h
