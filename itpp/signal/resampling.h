/*!
 * \file
 * \brief Resampling functions - header file
 * \author Tony Ottosson and Adam Piatyszek
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

#ifndef RESAMPLING_H
#define RESAMPLING_H

#include <itpp/base/mat.h>


namespace itpp {

  /*!
   * \addtogroup resampling
   * @{
   */

  //! Repeat each element in the vector \a norepeats times in sequence
  template<class T>
  Vec<T> repeat(const Vec<T> &v, int norepeats)
  {
    Vec<T> temp(v.length()*norepeats);

    for(int i=0; i<v.length(); i++) {
      for(int j=0;j<norepeats;j++)
	temp(i*norepeats+j)=v(i);
    }
    return temp;
  }

  //! Repeats each column \a norepeats times in sequence
  template<class T>
  Mat<T> repeat(const Mat<T> &m, int norepeats)
  {
    Mat<T> temp(m.rows(), m.cols()*norepeats);

    for (int j=0; j<m.cols(); j++) {
      for (int i=0;i<norepeats;i++) {
	temp.set_col(j*norepeats+i, m.get_col(j));
      }
    }
    return temp;
  }

  //! Upsample a vector by inserting \a (usf-1) zeros after each sample
  template<class T>
  void upsample(const Vec<T> &v, int usf, Vec<T> &u)
  {
    it_assert_debug(usf >= 1, "upsample: upsampling factor must be equal or greater than one" );
    u.set_size(v.length()*usf);
    u.clear();
    for(int i=0;i<v.length();i++)
      u(i*usf)=v(i);
  }


  //! Upsample a vector by inserting \a (usf-1) zeros after each sample
  template<class T>
  Vec<T> upsample(const Vec<T> &v, int usf)
  {
    Vec<T> u;
    upsample(v,usf,u);
    return u;
  }

  //! Upsample each column by inserting \a (usf-1) zeros after each column
  template<class T>
  void upsample(const Mat<T> &v, int usf, Mat<T> &u)
  {
    it_assert_debug(usf >= 1, "upsample: upsampling factor must be equal or greater than one" );
    u.set_size(v.rows(),v.cols()*usf);
    u.clear();
    for (int j=0;j<v.cols();j++)
      u.set_col(j*usf,v.get_col(j));
  }

  //! Upsample each column by inserting \a (usf-1) zeros after each column
  template<class T>
  Mat<T> upsample(const Mat<T> &v, int usf)
  {
    Mat<T> u;
    upsample(v,usf,u);
    return u;
  }

  //! Upsample each column by a factor of \a (usf-1) by linear interpolation
  template<class T>
  void lininterp(const Mat<T> &m, int usf, Mat<T> &u)
  {
    it_assert_debug(usf >= 1, "lininterp: upsampling factor must be equal or greater than one" );
    int L = (m.cols()-1)*usf+1;
    u.set_size(m.rows(),L);
    for (int i = 0; i < m.rows(); i++){
      for (int j = 0; j < L-1; j++)
	u(i,j) = (m(i,j/usf) + (j % usf)/((double)usf)*(m(i,(j+usf)/usf)-m(i,j/usf)));
      u(i,L-1) = m(i,m.cols()-1);
    }
  }

  /*!
   * \brief Upsample each column of matrix \a m to achieve \a f_ups
   * frequency using linear interpolation
   * \author Adam Piatyszek
   *
   * This function performs upsampling of matrix \a m to achieve
   * \a nrof_samples samples at \a f_ups frequency starting from the sample
   * at \a t_start time. The frequency of input samples stored in the matrix
   * \a m is defined by the \a f_base parameter.
   */
  template<class T>
  Mat<T> lininterp(const Mat<T> &m, double f_base, double f_ups,
		   int nrof_samples, double t_start = 0)
  {
    double t_base = 1 / f_base;
    double t_ups = 1 / f_ups;
    int rows = m.rows();
    int cols = m.cols();
    it_assert_debug(f_ups > f_base, "lininterp(): upsampled frequency must be greater than base frequency" );
    it_assert_debug((t_start >= 0) && (t_start < cols * t_base), "lininterp(): incorrect start time offset");
    it_assert_debug((nrof_samples * t_ups + t_start) <= (cols * t_base), "lininterp(): too many samples required or input data to short");
    Mat<T> u(rows, nrof_samples);
    double curr_time = t_start;

    int i = 0;
    int k = 0;
    while (i < cols - 1) {
      while ((curr_time < (i + 1) * t_base) && (k < nrof_samples)) {
	for (int j = 0; j < rows; j++) {
	  u(j, k) = (m(j, i) * ((i + 1) * t_base - curr_time)
		     - m(j, i + 1) * (i * t_base - curr_time)) / t_base;
	}
	k++;
	curr_time += t_ups;
      }
      i++;
    }
    return u;
  }

  //! Upsample each column by a factor of  \a (usf-1) by linear interpolation
  template<class T>
  Mat<T> lininterp(const Mat<T> &m, int usf)
  {
    Mat<T> u;
    lininterp(m,usf,u);
    return u;
  }

  //! Upsample by a factor of  \a (usf-1) by linear interpolation
  template<class T>
  void lininterp(const Vec<T> &v, int usf, Vec<T> &u)
  {
    it_assert_debug(usf >= 1, "lininterp(): upsampling factor must be equal or greater than one" );
    int L = (v.length()-1)*usf+1;
    u.set_size(L);
    for (int j = 0; j < L-1; j++) {
      u(j) = (v(j/usf) + (j % usf)/((double)usf)*(v((j+usf)/usf)-v(j/usf)));
    }
    u(L-1) = v(v.length()-1);
  }

  //! Upsample by a factor of  \a (usf-1) by linear interpolation
  template<class T>
  Vec<T> lininterp(const Vec<T> &v, int usf)
  {
    Vec<T> u;
    lininterp(v,usf,u);
    return u;
  }

  /*!
   * \brief Upsample vector \a v to achieve \a f_ups frequency using linear
   * interpolation
   * \author Adam Piatyszek
   *
   * This function performs upsampling of vector \a v to achieve
   * \a nrof_samples samples at \a f_ups frequency starting from the sample
   * at \a t_start time. The frequency of input samples stored in the vector
   * \a v is defined by the \a f_base parameter.
   */
  template<class T>
  Vec<T> lininterp(const Vec<T> &v, double f_base, double f_ups,
		   int nrof_samples, double t_start = 0)
  {
    double t_base = 1 / f_base;
    double t_ups = 1 / f_ups;
    int len = v.length();
    it_assert_debug(f_ups > f_base, "lininterp(): upsampled frequency must be greater than base frequency" );
    it_assert_debug((t_start >= 0) && (t_start < len * t_base), "lininterp(): incorrect start time offset");
    it_assert_debug((nrof_samples * t_ups + t_start) <= (len * t_base), "lininterp(): too many samples required or input data to short");
    Vec<T> u(nrof_samples);
    double curr_time = t_start;

    int i = 0;
    int k = 0;
    while (i < len - 1) {
      while ((curr_time < (i + 1) * t_base) && (k < nrof_samples)) {
	u(k) = (v(i) * ((i + 1) * t_base - curr_time)
		- v(i + 1) * (i * t_base - curr_time)) / t_base;
	k++;
	curr_time += t_ups;
      }
      i++;
    }
    return u;
  }

  /*!
   * @}
   */

#ifndef _MSC_VER

  // ----------------------------------------------------------------------
  // Instantiations
  // ----------------------------------------------------------------------

  //! Extern Template instantiation of repeat
  extern template vec repeat(const vec &v, int norepeats);
  //! Extern Template instantiation of repeat
  extern template cvec repeat(const cvec &v, int norepeats);
  //! Extern Template instantiation of repeat
  extern template svec repeat(const svec &v, int norepeats);
  //! Extern Template instantiation of repeat
  extern template ivec repeat(const ivec &v, int norepeats);
  //! Extern Template instantiation of repeat
  extern template bvec repeat(const bvec &v, int norepeats);

  //! Extern Template instantiation of repeat
  extern template mat repeat(const mat &m, int norepeats);
  //! Extern Template instantiation of repeat
  extern template cmat repeat(const cmat &m, int norepeats);
  //! Extern Template instantiation of repeat
  extern template smat repeat(const smat &m, int norepeats);
  //! Extern Template instantiation of repeat
  extern template imat repeat(const imat &m, int norepeats);
  //! Extern Template instantiation of repeat
  extern template bmat repeat(const bmat &m, int norepeats);

  //! Extern Template instantiation of upsample
  extern template vec upsample(const vec &v, int usf);
  //! Extern Template instantiation of upsample
  extern template cvec upsample(const cvec &v, int usf);
  //! Extern Template instantiation of upsample
  extern template svec upsample(const svec &v, int usf);
  //! Extern Template instantiation of upsample
  extern template ivec upsample(const ivec &v, int usf);
  //! Extern Template instantiation of upsample
  extern template bvec upsample(const bvec &v, int usf);

  //! Extern Template instantiation of upsample
  extern template mat upsample(const mat &v, int usf);
  //! Extern Template instantiation of upsample
  extern template cmat upsample(const cmat &v, int usf);
  //! Extern Template instantiation of upsample
  extern template smat upsample(const smat &v, int usf);
  //! Extern Template instantiation of upsample
  extern template imat upsample(const imat &v, int usf);
  //! Extern Template instantiation of upsample
  extern template bmat upsample(const bmat &v, int usf);

  //! Extern Template instantiation of upsample
  extern template void upsample(const vec &v, int usf,  vec &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const cvec &v, int usf,  cvec &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const svec &v, int usf,  svec &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const ivec &v, int usf,  ivec &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const bvec &v, int usf,  bvec &u);

  //! Extern Template instantiation of upsample
  extern template void upsample(const mat &v, int usf,  mat &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const cmat &v, int usf,  cmat &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const smat &v, int usf,  smat &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const imat &v, int usf,  imat &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const bmat &v, int usf,  bmat &u);

  //! Extern Template instantiation of liniterp
  extern template vec lininterp(const vec &v, int usf);
  //! Extern Template instantiation of liniterp
  extern template cvec lininterp(const cvec &v, int usf);

  //! Extern Template instantiation of liniterp
  extern template mat lininterp(const mat &v, int usf);
  //! Extern Template instantiation of liniterp
  extern template cmat lininterp(const cmat &v, int usf);

  //! Extern Template instantiation of liniterp
  extern template void lininterp(const vec &v, int usf,  vec &u);
  //! Extern Template instantiation of liniterp
  extern template void lininterp(const cvec &v, int usf,  cvec &u);

  //! Extern Template instantiation of liniterp
  extern template void lininterp(const mat &v, int usf,  mat &u);
  //! Extern Template instantiation of liniterp
  extern template void lininterp(const cmat &v, int usf,  cmat &u);

  //! Extern Template instantiation of liniterp
  extern template mat lininterp(const mat &m, double f_base, double f_ups, int nrof_samples, double t_start);
  //! Extern Template instantiation of liniterp
  extern template cmat lininterp(const cmat &m, double f_base, double f_ups, int nrof_samples, double t_start);

  //! Extern Template instantiation of liniterp
  extern template vec lininterp(const vec &v, double f_base, double f_ups, int nrof_samples, double t_start);
  //! Extern Template instantiation of liniterp
  extern template cvec lininterp(const cvec &v, double f_base, double f_ups, int nrof_samples, double t_start);

#endif

} // namespace itpp

#endif // #ifndef RESAMPLING_H

