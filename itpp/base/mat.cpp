/*!
 * \file
 * \brief Matrix Class Implementation
 * \author Tony Ottosson and Tobias Ringstrom
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

#include <itpp/base/mat.h>


namespace itpp {

  template<>
  bool Mat<std::complex<double> >::set(const char *values)
  {
    std::istringstream buffer(values);
    int rows=0, maxrows=10, cols=0, nocols=0, maxcols=10;

    alloc(maxrows, maxcols);

    while (buffer.peek()!=EOF) {
      rows++;
      if (rows > maxrows) {
	maxrows=maxrows*2;
	set_size(maxrows, maxcols, true);
      }

      cols=0;
      while ( (buffer.peek() != ';') && (buffer.peek() != EOF) ) {
	if (buffer.peek()==',') {
	  buffer.get();
	} else {
	  cols++;
	  if (cols > nocols) {
	    nocols=cols;
	    if (cols > maxcols) {
	      maxcols=maxcols*2;
	      set_size(maxrows, maxcols, true);
	    }
	  }
	  buffer >> this->operator()(rows-1,cols-1);
	  while (buffer.peek()==' ') { buffer.get(); }
	}
      }

      if (!buffer.eof())
	buffer.get();
    }
    set_size(rows, nocols, true);

    return true;
  }


  template<>
  bool Mat<bin>::set(const char *values)
  {
    std::istringstream buffer(values);
    int rows=0, maxrows=10, cols=0, nocols=0, maxcols=10;
    short intemp;

    alloc(maxrows, maxcols);

    while (buffer.peek()!=EOF) {
      rows++;
      if (rows > maxrows) {
	maxrows *= 2;
	set_size(maxrows, maxcols, true);
      }
      cols=0;
      while ( (buffer.peek() != ';') && (buffer.peek() != EOF) ) {
	if (buffer.peek()==',') {
	  buffer.get();
	} else {
	  cols++;
	  if (cols > nocols) {
	    nocols=cols;
	    if (cols > maxcols) {
	      maxcols=maxcols*2;
	      set_size(maxrows, maxcols, true);
	    }
	  }
	  buffer >> intemp;
	  this->operator()(rows-1,cols-1)=intemp;
	  while (buffer.peek()==' ') { buffer.get(); }
	}
      }
      if (!buffer.eof())
	buffer.get();
    }
    set_size(rows, nocols, true);

    return true;
  }

  template<>
  const cmat Mat<std::complex<double> >::hermitian_transpose() const
  {
    cmat temp(no_cols, no_rows);
    for (int i=0; i<no_rows; i++)
      for (int j=0; j<no_cols; j++)
	temp(j,i) = std::conj(operator()(i,j));

    return temp;
  }


  // -------- Multiplication operator -------------

#if defined(HAVE_CBLAS)
  template<>
  mat& Mat<double>::operator*=(const mat &m)
  {
    it_assert1(no_cols == m.no_rows,"mat::operator*=: wrong sizes");
    mat r(no_rows, m.no_cols); // unnecessary memory??

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, no_rows, m.no_cols, no_cols, 1.0,
		data, no_rows, m.data, m.no_rows, 0.0, r.data, r.no_rows);

    operator=(r); // time consuming 
    return *this;
  }

  template<>
  cmat& Mat<std::complex<double> >::operator*=(const cmat &m)
  {
    it_assert1(no_cols == m.no_rows,"cmat::operator*=: wrong sizes");
    std::complex<double> alpha = std::complex<double>(1.0), beta = std::complex<double>(0.0);
    cmat r(no_rows, m.no_cols); // unnecessary memory??

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, no_rows, m.no_cols, no_cols, (const void *)&alpha,
		data, no_rows, m.data, m.no_rows, (const void*)&beta, &r.data, r.no_rows);

    operator=(r); // time consuming
    return *this;
  }

  template<>
  const mat operator*(const mat &m1, const mat &m2)
  {
    it_assert1(m1.no_cols == m2.no_rows,"mat::operator*: wrong sizes");
    mat r(m1.no_rows, m2.no_cols);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m1.no_rows, m2.no_cols, m1.no_cols, 1.0,
		m1.data, m1.no_rows, m2.data, m2.no_rows, 0.0, r.data, r.no_rows);

    return r;
  }

  template<>
  const cmat operator*(const cmat &m1, const cmat &m2)
  {
    it_assert1(m1.no_cols == m2.no_rows,"cmat::operator*: wrong sizes");
    std::complex<double> alpha = std::complex<double>(1.0), beta = std::complex<double>(0.0);
    cmat r(m1.no_rows, m2.no_cols);

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m1.no_rows, m2.no_cols, m1.no_cols, (const void *)&alpha,
		m1.data, m1.no_rows, m2.data, m2.no_rows, (const void*)&beta, r.data, r.no_rows);

    return r;
  }

  template<>
  const Vec<double> operator*(const mat &m, const Vec<double> &v)
  {
    it_assert1(m.no_cols == v.size(),"mat::operator*: wrong sizes");
    vec r(m.no_rows);

    cblas_dgemv(CblasColMajor, CblasNoTrans, m.no_rows, m.no_cols, 1.0, m.data, m.no_rows, v._data(), 1,
		0.0, r._data(), 1);

    return r;
  }

  template<>
  const Vec<std::complex<double> > operator*(const cmat &m, const Vec<std::complex<double> > &v)
  {
    it_assert1(m.no_cols == v.size(),"cmat::operator*: wrong sizes");
    std::complex<double> alpha = std::complex<double>(1.0), beta = std::complex<double>(0.0);
    cvec r(m.no_rows);

    cblas_zgemv(CblasColMajor, CblasNoTrans, m.no_rows, m.no_cols, (const void *)&alpha, m.data, m.no_rows,
		v._data(), 1, (const void*)&beta, r._data(), 1);

    return r;
  }

#endif // HAVE_CBLAS


  // ---------------------- Instantiations --------------------------------

  //------- class instantiations --------

  template class Mat<double>;
  template class Mat<std::complex<double> >;
  template class Mat<int>;
  template class Mat<short int>;
  template class Mat<bin>;


  //-------------------- Operator instantiations --------------------

  //-------- Addition operators ---------------

  template const mat operator+(const mat &m1, const mat &m2);
  template const cmat operator+(const cmat &m1, const cmat &m2);
  template const imat operator+(const imat &m1, const imat &m2);
  template const smat operator+(const smat &m1, const smat &m2);
  template const bmat operator+(const bmat &m1, const bmat &m2);

  template const mat operator+(const mat &m, double t);
  template const cmat operator+(const cmat &m, std::complex<double> t);
  template const imat operator+(const imat &m, int t);
  template const smat operator+(const smat &m, short t);
  template const bmat operator+(const bmat &m, bin t);

  template const mat operator+(double t, const mat &m);
  template const cmat operator+(std::complex<double> t, const cmat &m);
  template const imat operator+(int t, const imat &m);
  template const smat operator+(short t, const smat &m);
  template const bmat operator+(bin t, const bmat &m);

  //-------- Subraction operators ---------------

  template const mat operator-(const mat &m1, const mat &m2);
  template const cmat operator-(const cmat &m1, const cmat &m2);
  template const imat operator-(const imat &m1, const imat &m2);
  template const smat operator-(const smat &m1, const smat &m2);
  template const bmat operator-(const bmat &m1, const bmat &m2);

  template const mat operator-(const mat &m, double t);
  template const cmat operator-(const cmat &m, std::complex<double> t);
  template const imat operator-(const imat &m, int t);
  template const smat operator-(const smat &m, short t);
  template const bmat operator-(const bmat &m, bin t);

  template const mat operator-(double t, const mat &m);
  template const cmat operator-(std::complex<double> t, const cmat &m);
  template const imat operator-(int t, const imat &m);
  template const smat operator-(short t, const smat &m);
  template const bmat operator-(bin t, const bmat &m);

  //--------- Unary minus ---------------

  template const mat operator-(const mat &m);
  template const cmat operator-(const cmat &m);
  template const imat operator-(const imat &m);
  template const smat operator-(const smat &m);
  template const bmat operator-(const bmat &m);

  //-------- Multiplication operators ---------------

#if !defined(HAVE_CBLAS)
  template const mat operator*(const mat &m1, const mat &m2);
  template const cmat operator*(const cmat &m1, const cmat &m2);
#endif
  template const imat operator*(const imat &m1, const imat &m2);
  template const smat operator*(const smat &m1, const smat &m2);
  template const bmat operator*(const bmat &m1, const bmat &m2);

#if !defined(HAVE_CBLAS)
  template const vec operator*(const mat &m, const vec &v);
  template const cvec operator*(const cmat &m, const cvec &v);
#endif
  template const ivec operator*(const imat &m, const ivec &v);
  template const svec operator*(const smat &m, const svec &v);
  template const bvec operator*(const bmat &m, const bvec &v);

  template const mat operator*(const vec &v, const mat &m);
  template const cmat operator*(const cvec &v, const cmat &m);
  template const imat operator*(const ivec &v, const imat &m);
  template const smat operator*(const svec &v, const smat &m);
  template const bmat operator*(const bvec &v, const bmat &m);

  template const mat operator*(const mat &m, double t);
  template const cmat operator*(const cmat &m, std::complex<double> t);
  template const imat operator*(const imat &m, int t);
  template const smat operator*(const smat &m, short t);
  template const bmat operator*(const bmat &m, bin t);

  template const mat operator*(double t, const mat &m);
  template const cmat operator*(std::complex<double> t, const cmat &m);
  template const imat operator*(int t, const imat &m);
  template const smat operator*(short t, const smat &m);
  template const bmat operator*(bin t, const bmat &m);

  // ------------ Elementwise multiplication -----------

  template const mat elem_mult(const mat &m1, const mat &m2);
  template const cmat elem_mult(const cmat &m1, const cmat &m2);
  template const imat elem_mult(const imat &m1, const imat &m2);
  template const smat elem_mult(const smat &m1, const smat &m2);
  template const bmat elem_mult(const bmat &m1, const bmat &m2);

  // ------------ Division operator -----------

  template const mat operator/(const mat &m, double t);
  template const cmat operator/(const cmat &m, std::complex<double> t);
  template const imat operator/(const imat &m, int t);
  template const smat operator/(const smat &m, short t);
  template const bmat operator/(const bmat &m, bin t);

  // ------------ Elementwise division -----------

  template const mat elem_div(const mat &m1, const mat &m2);
  template const cmat elem_div(const cmat &m1, const cmat &m2);
  template const imat elem_div(const imat &m1, const imat &m2);
  template const smat elem_div(const smat &m1, const smat &m2);
  template const bmat elem_div(const bmat &m1, const bmat &m2);

  // ------------- Concatenations -----------------

  template const mat concat_horizontal(const mat &m1, const mat &m2);
  template const cmat concat_horizontal(const cmat &m1, const cmat &m2);
  template const imat concat_horizontal(const imat &m1, const imat &m2);
  template const smat concat_horizontal(const smat &m1, const smat &m2);
  template const bmat concat_horizontal(const bmat &m1, const bmat &m2);

  template const mat concat_vertical(const mat &m1, const mat &m2);
  template const cmat concat_vertical(const cmat &m1, const cmat &m2);
  template const imat concat_vertical(const imat &m1, const imat &m2);
  template const smat concat_vertical(const smat &m1, const smat &m2);
  template const bmat concat_vertical(const bmat &m1, const bmat &m2);

  //----------- Output stream --------------

  template std::ostream &operator<<(std::ostream &os, const mat  &m);
  template std::ostream &operator<<(std::ostream &os, const cmat &m);
  template std::ostream &operator<<(std::ostream &os, const imat  &m);
  template std::ostream &operator<<(std::ostream &os, const smat  &m);
  template std::ostream &operator<<(std::ostream &os, const bmat  &m);

  //----------- Input stream -------------

  template std::istream &operator>>(std::istream &is, mat  &m);
  template std::istream &operator>>(std::istream &is, cmat &m);
  template std::istream &operator>>(std::istream &is, imat  &m);
  template std::istream &operator>>(std::istream &is, smat  &m);
  template std::istream &operator>>(std::istream &is, bmat  &m);

} // namespace itpp
