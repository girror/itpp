/*!
 * \file
 * \brief Eigenvalue decomposition test program
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2008  (see AUTHORS file for a list of contributors)
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

#include <itpp/itstat.h>

using namespace itpp;
using namespace std;


int main(void)
{
  cout << "================================" << endl;
  cout << "  Test of eigenvalue routines   " << endl;
  cout << "================================" << endl;

  const double thresh = 1e-13;

  {
    cout << "Real symmetric matrix" << endl;
    mat A = randn(5, 5);
    A = transpose(A) * A; // make it symmetic
    mat V;
    vec d;
    eig_sym(A, d, V);

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A*V-V*diag(d)) = "
	 << round_to_zero(norm(A * V - V * diag(d)), thresh) << endl;
  }

  {
    cout << endl << "Real non-symmetric matrix" << endl;
    mat A = randn(5, 5);
    cmat V;
    cvec d;
    eig(A, d, V);

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A*V-V*diag(d)) = "
	 << round_to_zero(norm(A * V - V * diag(d)), thresh) << endl;
  }

  {
    cout << endl << "Complex hermitian matrix" << endl;
    cmat A = randn_c(5, 5);
    A = transpose(conj(A)) * A; // make it hermitian
    cmat V;
    vec d;
    eig_sym(A, d, V);

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A*V-V*diag(d)) = "
	 << round_to_zero(norm(A * V - V * diag(d)), thresh) << endl;
  }

  {
    cout << endl << "Complex non-hermitian matrix" << endl;
    cmat A = randn_c(5, 5);
    cmat V;
    cvec d;
    eig(A, d, V);

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A*V-V*diag(d)) = "
	 << round_to_zero(norm(A * V - V * diag(d)), thresh) << endl;
  }

  return 0;
}
