/*!
 * \file
 * \brief Sparse Vector Class implementation
 * \author Tony Ottosson and Tobias Ringstrom
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

#include <itpp/base/svec.h>


namespace itpp {

  /*--------------------------------------------------------------
   * Explicit initializations
   *--------------------------------------------------------------*/

  template class Sparse_Vec<int>;
  template class Sparse_Vec<double>;
  template class Sparse_Vec<std::complex<double> >;

  template sparse_ivec operator+(const sparse_ivec &, const sparse_ivec &);
  template sparse_vec operator+(const sparse_vec &, const sparse_vec &);
  template sparse_cvec operator+(const sparse_cvec &, const sparse_cvec &);

  template int operator*(const sparse_ivec &, const sparse_ivec &);
  template double operator*(const sparse_vec &, const sparse_vec &);
  template std::complex<double> operator*(const sparse_cvec &, const sparse_cvec &);

  template int operator*(const sparse_ivec &, const ivec &);
  template double operator*(const sparse_vec &, const vec &);
  template std::complex<double> operator*(const sparse_cvec &, const cvec &);

  template int operator*(const ivec &, const sparse_ivec &);
  template double operator*(const vec &, const sparse_vec &);
  template std::complex<double> operator*(const cvec &, const sparse_cvec &);

  template sparse_ivec elem_mult(const sparse_ivec &, const sparse_ivec &);
  template sparse_vec elem_mult(const sparse_vec &, const sparse_vec &);
  template sparse_cvec elem_mult(const sparse_cvec &, const sparse_cvec &);

  template ivec elem_mult(const sparse_ivec &, const ivec &);
  template vec elem_mult(const sparse_vec &, const vec &);
  template cvec elem_mult(const sparse_cvec &, const cvec &);

  template sparse_ivec elem_mult_s(const sparse_ivec &, const ivec &);
  template sparse_vec elem_mult_s(const sparse_vec &, const vec &);
  template sparse_cvec elem_mult_s(const sparse_cvec &, const cvec &);

  template ivec elem_mult(const ivec &, const sparse_ivec &);
  template vec elem_mult(const vec &, const sparse_vec &);
  template cvec elem_mult(const cvec &, const sparse_cvec &);

  template sparse_ivec elem_mult_s(const ivec &, const sparse_ivec &);
  template sparse_vec elem_mult_s(const vec &, const sparse_vec &);
  template sparse_cvec elem_mult_s(const cvec &, const sparse_cvec &);

} // namespace itpp
