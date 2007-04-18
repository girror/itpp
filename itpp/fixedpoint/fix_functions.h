/*!
 * \file 
 * \brief Definitions of a set of functions for Fix, Fixed, CFix and 
 * CFixed classes
 * \author Johan Bergman
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

#ifndef FIX_FUNCTIONS_H
#define FIX_FUNCTIONS_H

#include <itpp/fixedpoint/cfix.h>
#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/array.h>


using std::complex;

namespace itpp {

  //! \addtogroup fixtypes
  //!@{

  //! Return true only if argument is of type Fix or CFix (or an Array/Vec/Mat of Fix or CFix)
  template<class T> inline bool is_fix(const T &x) {return false;}
  //! Return true only if argument is of type Fix or CFix (or an Array/Vec/Mat of Fix or CFix)
  template<> inline bool is_fix(const Fix &x) {return true;}
  //! Return true only if argument is of type Fix or CFix (or an Array/Vec/Mat of Fix or CFix)
  template<> inline bool is_fix(const fixvec &x) {return true;}
  //! Return true only if argument is of type Fix or CFix (or an Array/Vec/Mat of Fix or CFix)
  template<> inline bool is_fix(const fixmat &x) {return true;}
  //! Return true only if argument is of type Fix or CFix (or an Array/Vec/Mat of Fix or CFix)
  template<> inline bool is_fix(const CFix &x) {return true;}
  //! Return true only if argument is of type Fix or CFix (or an Array/Vec/Mat of Fix or CFix)
  template<> inline bool is_fix(const cfixvec &x) {return true;}
  //! Return true only if argument is of type Fix or CFix (or an Array/Vec/Mat of Fix or CFix)
  template<> inline bool is_fix(const cfixmat &x) {return true;}
  //! Return true only if argument is of type Fix or CFix (or an Array/Vec/Mat of Fix or CFix)
  template<class T> inline bool is_fix(const Array<T> &x) {return is_fix(T());}

  //! Set <tt>y = x * pow2(n)</tt> using the quantization mode of \c y
  inline void set_fix(Fix &y, double x, int n) {y.set(x, n);}
  //! Set <tt>y = x * pow2(n)</tt> using the specified quantization mode \c q
  inline void set_fix(Fix &y, double x, int n, q_mode q) {y.set(x, n, q);}
  //! Set <tt>y = x * pow2(n)</tt> using the quantization mode of \c y
  inline void set_fix(fixvec &y, const vec &x, int n)
  {
    y.set_size(x.length());
    for (int i=0; i<y.size(); i++) y(i).set(x(i), n);
  }
  //! Set <tt>y = x * pow2(n)</tt> using the specified quantization mode \c q
  inline void set_fix(fixvec &y, const vec &x, int n, q_mode q)
  {
    y.set_size(x.length());
    for (int i=0; i<y.size(); i++) y(i).set(x(i), n, q);
  }
  //! Set <tt>y = x * pow2(n)</tt> using the quantization mode of \c y
  inline void set_fix(fixmat &y, const mat &x, int n)
  {
    y.set_size(x.rows(), x.cols());
    for (int i=0; i<y.size(); i++) y(i).set(x(i), n);
  }
  //! Set <tt>y = x * pow2(n)</tt> using the specified quantization mode \c q
  inline void set_fix(fixmat &y, const mat &x, int n, q_mode q)
  {
    y.set_size(x.rows(), x.cols());
    for (int i=0; i<y.size(); i++) y(i).set(x(i), n, q);
  }
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(double &y, double x, int n) {y = x;}
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(double &y, double x, int n, q_mode q) {y = x;}
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(vec &y, const vec &x, int n) {y = x;}
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(vec &y, const vec &x, int n, q_mode q) {y = x;}
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(mat &y, const mat &x, int n) {y = x;}
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(mat &y, const mat &x, int n, q_mode q) {y = x;}

  //! Set <tt>y = x * pow2(n)</tt> using the quantization mode of \c y
  inline void set_fix(CFix &y, complex<double> x, int n) {y.set(x, n);}
  //! Set <tt>y = (real + i*imag) * pow2(n)</tt> using the quantization mode of \c y
  inline void set_fix(CFix &y, double real, double imag, int n) {y.set(real, imag, n);}
  //! Set <tt>y = x * pow2(n)</tt> using the specified quantization mode \c q
  inline void set_fix(CFix &y, complex<double> x, int n, q_mode q) {y.set(x, n, q);}
  //! Set <tt>y = (real + i*imag) * pow2(n)</tt> using the specified quantization mode \c q
  inline void set_fix(CFix &y, double real, double imag, int n, q_mode q) {y.set(real, imag, n, q);}
  //! Set <tt>y = x * pow2(n)</tt> using the quantization mode of \c y
  inline void set_fix(cfixvec &y, const cvec &x, int n)
  {
    y.set_size(x.length());
    for (int i=0; i<y.size(); i++) y(i).set(x(i), n);
  }
  //! Set <tt>y = (real + i*imag) * pow2(n)</tt> using the quantization mode of \c y
  inline void set_fix(cfixvec &y, const vec &real, const vec &imag, int n)
  {
    it_assert1(real.length()==imag.length(), "set_fix: real and imag should have the same size");
    y.set_size(real.length());
    for (int i=0; i<y.size(); i++) y(i).set(real(i), imag(i), n);
  }
  //! Set <tt>y = x * pow2(n)</tt> using the specified quantization mode \c q
  inline void set_fix(cfixvec &y, const cvec &x, int n, q_mode q)
  {
    y.set_size(x.length());
    for (int i=0; i<y.size(); i++) y(i).set(x(i), n, q);
  }
  //! Set <tt>y = (real + i*imag) * pow2(n)</tt> using the specified quantization mode \c q
  inline void set_fix(cfixvec &y, const vec &real, const vec &imag, int n, q_mode q)
  {
    it_assert1(real.length()==imag.length(), "set_fix: real and imag should have the same size");
    y.set_size(real.length());
    for (int i=0; i<y.size(); i++) y(i).set(real(i), imag(i), n, q);
  }
  //! Set <tt>y = x * pow2(n)</tt> using the quantization mode of \c y
  inline void set_fix(cfixmat &y, const cmat &x, int n)
  {
    y.set_size(x.rows(), x.cols());
    for (int i=0; i<y.size(); i++) y(i).set(x(i), n);
  }
  //! Set <tt>y = (real + i*imag) * pow2(n)</tt> using the quantization mode of \c y
  inline void set_fix(cfixmat &y, const mat &real, const mat &imag, int n)
  {
    it_assert1(real.rows()==imag.rows() && real.cols()==imag.cols(), "set_fix: real and imag should have the same size");
    y.set_size(real.rows(), real.cols());
    for (int i=0; i<y.size(); i++) y(i).set(real(i), imag(i), n);
  }
  //! Set <tt>y = x * pow2(n)</tt> using the specified quantization mode \c q
  inline void set_fix(cfixmat &y, const cmat &x, int n, q_mode q)
  {
    y.set_size(x.rows(), x.cols());
    for (int i=0; i<y.size(); i++) y(i).set(x(i), n, q);
  }
  //! Set <tt>y = (real + i*imag) * pow2(n)</tt> using the specified quantization mode \c q
  inline void set_fix(cfixmat &y, const mat &real, const mat &imag, int n, q_mode q)
  {
    it_assert1(real.rows()==imag.rows() && real.cols()==imag.cols(), "set_fix: real and imag should have the same size");
    y.set_size(real.rows(), real.cols());
    for (int i=0; i<y.size(); i++) y(i).set(real(i), imag(i), n, q);
  }
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(complex<double> &y, const complex<double> &x, int n) {y = x;}
  //! Set <tt>y = real + i*imag</tt>. Useful in templated code
  inline void set_fix(complex<double> &y, double real, double imag, int n) {y = complex<double>(real, imag);}
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(complex<double> &y, const complex<double> &x, int n, q_mode q) {y = x;}
  //! Set <tt>y = real + i*imag</tt>. Useful in templated code
  inline void set_fix(complex<double> &y, double real, double imag, int n, q_mode q) {y = complex<double>(real, imag);}
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(cvec &y, const cvec &x, int n) {y = x;}
  //! Set <tt>y = real + i*imag</tt>. Useful in templated code
  inline void set_fix(cvec &y, const vec &real, const vec &imag, int n) {y = to_cvec(real, imag);}
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(cvec &y, const cvec &x, int n, q_mode q) {y = x;}
  //! Set <tt>y = real + i*imag</tt>. Useful in templated code
  inline void set_fix(cvec &y, const vec &real, const vec &imag, int n, q_mode q) {y = to_cvec(real, imag);}
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(cmat &y, const cmat &x, int n) {y = x;}
  //! Set <tt>y = real + i*imag</tt>. Useful in templated code
  inline void set_fix(cmat &y, const mat &real, const mat &imag, int n) {y = to_cmat(real, imag);}
  //! Set <tt>y = x</tt>. Useful in templated code
  inline void set_fix(cmat &y, const cmat &x, int n, q_mode q) {y = x;}
  //! Set <tt>y = real + i*imag</tt>. Useful in templated code
  inline void set_fix(cmat &y, const mat &real, const mat &imag, int n, q_mode q) {y = to_cmat(real, imag);}

  //! Call set_fix for each Array element
  template<class T1, class T2> inline void set_fix(Array<T1> &y, const Array<T2> &x, int n)
  {
    y.set_size(x.size());
    for (int i=0; i<y.size(); i++) set_fix(y(i), x(i), n);
  }
  //! Call set_fix for each Array element
  template<class T1, class T2> inline void set_fix(Array<T1> &y, const Array<T2> &real, const Array<T2> &imag, int n)
  {
    it_assert1(real.size()==imag.size(), "set_fix: real and imag should have the same size");
    y.set_size(real.size());
    for (int i=0; i<y.size(); i++) set_fix(y(i), real(i), imag(i), n);
  }
  //! Call set_fix for each Array element
  template<class T1, class T2> inline void set_fix(Array<T1> &y, const Array<T2> &x, int n, q_mode q)
  {
    y.set_size(x.size());
    for (int i=0; i<y.size(); i++) set_fix(y(i), x(i), n, q);
  }
  //! Call set_fix for each Array element
  template<class T1, class T2> inline void set_fix(Array<T1> &y, const Array<T2> &real, const Array<T2> &imag, int n, q_mode q)
  {
    it_assert1(real.size()==imag.size(), "set_fix: real and imag should have the same size");
    y.set_size(real.size());
    for (int i=0; i<y.size(); i++) set_fix(y(i), real(i), imag(i), n, q);
  }

  //! Left shift \c n bits
  inline void lshift_fix(Fix &y, int n) {y.lshift(n);}
  //! Right shift \c n bits using the quantization mode of \c y
  inline void rshift_fix(Fix &y, int n) {y.rshift(n);}
  //! Right shift \c n bits using the specified quantization mode \c q
  inline void rshift_fix(Fix &y, int n, q_mode q) {y.rshift(n, q);}
  //! Left shift \c n bits
  inline void lshift_fix(fixvec &y, int n)
  {for (int i=0; i<y.size(); i++) y(i).lshift(n);}
  //! Right shift \c n bits using the quantization mode of \c y
  inline void rshift_fix(fixvec &y, int n)
  {for (int i=0; i<y.size(); i++) y(i).rshift(n);}
  //! Right shift \c n bits using the specified quantization mode \c q
  inline void rshift_fix(fixvec &y, int n, q_mode q)
  {for (int i=0; i<y.size(); i++) y(i).rshift(n, q);}
  //! Left shift \c n bits
  inline void lshift_fix(fixmat &y, int n)
  {for (int i=0; i<y.size(); i++) y(i).lshift(n);}
  //! Right shift \c n bits using the quantization mode of \c y
  inline void rshift_fix(fixmat &y, int n)
  {for (int i=0; i<y.size(); i++) y(i).rshift(n);}
  //! Right shift \c n bits using the specified quantization mode \c q
  inline void rshift_fix(fixmat &y, int n, q_mode q)
  {for (int i=0; i<y.size(); i++) y(i).rshift(n, q);}
  //! Dummy function useful in templated code
  inline void lshift_fix(double &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(double &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(double &y, int n, q_mode q) {}
  //! Dummy function useful in templated code
  inline void lshift_fix(vec &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(vec &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(vec &y, int n, q_mode q) {}
  //! Dummy function useful in templated code
  inline void lshift_fix(mat &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(mat &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(mat &y, int n, q_mode q) {}
  //! Left shift \c n bits
  inline void lshift_fix(CFix &y, int n) {y.lshift(n);}
  //! Right shift \c n bits using the quantization mode of \c y
  inline void rshift_fix(CFix &y, int n) {y.rshift(n);}
  //! Right shift \c n bits using the specified quantization mode \c q
  inline void rshift_fix(CFix &y, int n, q_mode q) {y.rshift(n, q);}
  //! Left shift \c n bits
  inline void lshift_fix(cfixvec &y, int n)
  {for (int i=0; i<y.size(); i++) y(i).lshift(n);}
  //! Right shift \c n bits using the quantization mode of \c y
  inline void rshift_fix(cfixvec &y, int n)
  {for (int i=0; i<y.size(); i++) y(i).rshift(n);}
  //! Right shift \c n bits using the specified quantization mode \c q
  inline void rshift_fix(cfixvec &y, int n, q_mode q)
  {for (int i=0; i<y.size(); i++) y(i).rshift(n, q);}
  //! Left shift \c n bits
  inline void lshift_fix(cfixmat &y, int n)
  {for (int i=0; i<y.size(); i++) y(i).lshift(n);}
  //! Right shift \c n bits using the quantization mode of \c y
  inline void rshift_fix(cfixmat &y, int n)
  {for (int i=0; i<y.size(); i++) y(i).rshift(n);}
  //! Right shift \c n bits using the specified quantization mode \c q
  inline void rshift_fix(cfixmat &y, int n, q_mode q)
  {for (int i=0; i<y.size(); i++) y(i).rshift(n, q);}
  //! Dummy function useful in templated code
  inline void lshift_fix(complex<double> &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(complex<double> &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(complex<double> &y, int n, q_mode q) {}
  //! Dummy function useful in templated code
  inline void lshift_fix(cvec &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(cvec &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(cvec &y, int n, q_mode q) {}
  //! Dummy function useful in templated code
  inline void lshift_fix(cmat &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(cmat &y, int n) {}
  //! Dummy function useful in templated code
  inline void rshift_fix(cmat &y, int n, q_mode q) {}
  //! Call lshift_fix for each Array element
  template<class T> inline void lshift_fix(Array<T> &y, int n)
  {for (int i=0; i<y.size(); i++) lshift_fix(y(i), n);}
  //! Call rshift_fix for each Array element
  template<class T> inline void rshift_fix(Array<T> &y, int n)
  {for (int i=0; i<y.size(); i++) rshift_fix(y(i), n);}
  //! Call rshift_fix for each Array element
  template<class T> inline void rshift_fix(Array<T> &y, int n, q_mode q)
  {for (int i=0; i<y.size(); i++) rshift_fix(y(i), n, q);}

  //! If x is a fixed-point variable, assert that x has the specified shift value, otherwise do nothing
  inline void assert_fixshift(const double x, const int shift) {}
  //! If x is a fixed-point variable, assert that x has the specified shift value, otherwise do nothing
  inline void assert_fixshift(const complex<double> &x, const int shift) {}
  //! If x is a fixed-point variable, assert that x has the specified shift value, otherwise do nothing
  inline void assert_fixshift(const Fix &x, const int shift)
  {it_assert1(x.get_shift()==shift, "Shift should be " + to_str(shift) + " but it is " + to_str(x.get_shift()) + ".");}
  //! If x is a fixed-point variable, assert that x has the specified shift value, otherwise do nothing
  inline void assert_fixshift(const CFix &x, const int shift)
  {it_assert1(x.get_shift()==shift, "Shift should be " + to_str(shift) + " but it is " + to_str(x.get_shift()) + ".");}

  //! Converts a fixvec to vec
  vec to_vec(const fixvec &v);
  //! Converts a cfixvec to cvec
  cvec to_cvec(const cfixvec &v);
  //! Converts a fixmat to mat
  mat to_mat(const fixmat &m);
  //! Converts a cfixmat to cmat
  cmat to_cmat(const cfixmat &m);

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  //! Help class used by the conversion function to<T>(const Array<...> &x). To be merged with Convert?
  template<class T, class U>
  class ConvertU2T {
  public:
    typedef T result;
  };
  //! Template specialization for Array<T>
  template<class T, class U>
  class ConvertU2T<T, Array<U> > {
  public:
    typedef Array<typename ConvertU2T<T, U>::result> result;  // Recursive
  };
  //! Template specialization for Vec<T>
  template<class T, class U>
  class ConvertU2T<T, Vec<U> > {
  public:
    typedef Vec<T> result;
  };
  //! Template specialization for Mat<T>
  template<class T, class U>
  class ConvertU2T<T, Mat<U> > {
  public:
    typedef Mat<T> result;
  };

#endif

  //! Convert double to T
  template<class T> inline T to(const double x) {return T(x);}
  //! Convert Fix to T
  template<class T> inline T to(const Fix &x) {return T(x);}
  //! Convert complex<double> to T
  template<class T> inline T to(const complex<double> &x) {return T(x);}
  //! Convert CFix to T
  template<class T> inline T to(const CFix &x) {return T(x);}
  //! Convert double (real and imaginary parts) to T
  template<class T> inline T to(const double real, const double imag) {return T(real, imag);}
  //! Convert Fix (real and imaginary parts) to T
  template<class T> inline T to(const Fix &real, const Fix &imag) {return T(real, imag);}

  //! Convert Vec<U> to Vec<T>
  template<class T, class U> Vec<T> to(const Vec<U> &x)
  {
    Vec<T> y(x.length());
    for (int i=0; i<x.length(); i++) {
      y(i) = T(x(i));
    }
    return y;
  }
  //! Convert vec to vec
  template<> inline vec to<double>(const vec &x) {return x;}
  //! Convert cvec to cvec
  template<> inline cvec to<complex<double> >(const cvec &x) {return x;}
  //! Convert fixvec to fixvec
  template<> inline fixvec to<Fix>(const fixvec &x) {return x;}
  //! Convert cfixvec to cfixvec
  template<> inline cfixvec to<CFix>(const cfixvec &x) {return x;}

  //! Convert Vec<U> (real and imaginary parts) to Vec<T>
  template<class T, class U> Vec<T> to(const Vec<U> &real, const Vec<U> &imag)
  {
    it_assert1(real.length()==imag.length(), "to: real and imag should have the same size");
    Vec<T> y(real.length());
    for (int i=0; i<real.length(); i++) {
      y(i) = T(real(i), imag(i));
    }
    return y;
  }

  //! Convert Mat<U> to Mat<T>
  template<class T, class U> Mat<T> to(const Mat<U> &x)
  {
    Mat<T> y(x.rows(), x.cols());
    for (int i=0; i<x.rows(); i++) {
      for (int j=0; j<x.cols(); j++) {
        y(i,j) = T(x(i,j));
      }
    }
    return y;
  }
  //! Convert mat to mat
  template<> inline mat to<double>(const mat &x) {return x;}
  //! Convert cmat to cmat
  template<> inline cmat to<complex<double> >(const cmat &x) {return x;}
  //! Convert fixmat to fixmat
  template<> inline fixmat to<Fix>(const fixmat &x) {return x;}
  //! Convert cfixmat to cfixmat
  template<> inline cfixmat to<CFix>(const cfixmat &x) {return x;}

  //! Convert Mat<U> (real and imaginary parts) to Mat<T>
  template<class T, class U> Mat<T> to(const Mat<U> &real, const Mat<U> &imag)
  {
    it_assert1(real.rows()==imag.rows() && real.cols()==imag.cols(), "to: real and imag should have the same size");
    Mat<T> y(real.rows(), real.cols());
    for (int i=0; i<real.rows(); i++) {
      for (int j=0; j<real.cols(); j++) {
        y(i,j) = T(real(i,j), imag(i,j));
      }
    }
    return y;
  }

  //! Convert Array<U>, where U can be an Array/Vec/Mat, to a corresponding Array with T elements
  template<class T, class U>
  Array<typename ConvertU2T<T, U>::result> to(const Array<U> &x)
  {
    Array<typename ConvertU2T<T, U>::result> y(x.size());
    for (int i=0; i<x.size(); i++) {
      y(i) = to<T>(x(i));
    }
    return y;
  }

  //! Convert Array<U> (real and imaginary parts), where U can be an Array/Vec/Mat, to a corresponding Array with T elements
  template<class T, class U>
  Array<typename ConvertU2T<T, U>::result> to(const Array<U> &real, const Array<U> &imag)
  {
    it_assert1(real.size()==imag.size(), "to: real and imag should have the same size");
    Array<typename ConvertU2T<T, U>::result> y(real.size());
    for (int i=0; i<real.size(); i++) {
      y(i) = to<T>(real(i), imag(i));
    }
    return y;
  }

  //! Convert Fix to double by multiplying the bit representation with pow2(-shift)
  inline double unfix(const Fix &x) {return x.unfix();}
  //! Convert CFix to complex<double> by multiplying the bit representation with pow2(-shift)
  inline complex<double> unfix(const CFix &x) {return x.unfix();}
  //! Convert fixvec to vec by multiplying the bit representations with pow2(-shift)
  inline vec unfix(const fixvec &x) {return to_vec(x);}
  //! Convert cfixvec to cvec by multiplying the bit representations with pow2(-shift)
  inline cvec unfix(const cfixvec &x) {return to_cvec(x);}
  //! Convert fixmat to mat by multiplying the bit representations with pow2(-shift)
  inline mat unfix(const fixmat &x) {return to_mat(x);}
  //! Convert cfixmat to cmat by multiplying the bit representations with pow2(-shift)
  inline cmat unfix(const cfixmat &x) {return to_cmat(x);}

  //! Convert double to double i.e. do nothing
  inline double unfix(const double &x) {return x;}
  //! Convert complex<double> to complex<double> i.e. do nothing
  inline complex<double> unfix(const complex<double> &x) {return x;}
  //! Convert vec to vec i.e. do nothing
  inline vec unfix(const vec &x) {return x;}
  //! Convert cvec to cvec i.e. do nothing
  inline cvec unfix(const cvec &x) {return x;}
  //! Convert mat to mat i.e. do nothing
  inline mat unfix(const mat &x) {return x;}
  //! Convert cmat to cmat i.e. do nothing
  inline cmat unfix(const cmat &x) {return x;}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  //! Help class used by the conversion function unfix(const Array<T> &x)
  template<class T>
  class Convert {
  public:
    typedef double to_double;
  };
  //! Template specialization for CFix
  template<>
  class Convert<CFix> {
  public:
    typedef complex<double> to_double;
  };
  //! Template specialization for complex<T>
  template<class T>
  class Convert<complex<T> > {
  public:
    typedef complex<double> to_double;
  };
  //! Template specialization for Array<T>
  template<class T>
  class Convert<Array<T> > {
  public:
    typedef Array<typename Convert<T>::to_double> to_double;  // Recursive
  };
  //! Template specialization for Vec<T>
  template<class T>
  class Convert<Vec<T> > {
  public:
    typedef Vec<typename Convert<T>::to_double> to_double;  // Recursive
  };
  //! Template specialization for Mat<T>
  template<class T>
  class Convert<Mat<T> > {
  public:
    typedef Mat<typename Convert<T>::to_double> to_double;  // Recursive
  };

#endif //DOXYGEN_SHOULD_SKIP_THIS

  //! Convert floating- or fixed-point Array to floating-point Array
  template<class T>
  Array<typename Convert<T>::to_double> unfix(const Array<T> &x)
  {
    Array<typename Convert<T>::to_double> y(x.size());
    for (int i=0; i<x.size(); i++) {
      y(i) = unfix(x(i));
    }
    return y;
  }

  //! Absolute value
  Fix abs(const Fix &x);
  //! Real part of complex value
  Fix real(const CFix &x);
  //! Imaginary part of complex value
  Fix imag(const CFix &x);
  //! Conjugate of complex value
  CFix conj(const CFix &x);

  //!@}

} // namespace itpp

#endif // #ifndef FIX_FUNCTIONS_H
