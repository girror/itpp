/*!
 * \file
 * \brief LU factorization routines test program
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
  cout << "=================================" << endl;
  cout << "Test of LU factorization routines" << endl;
  cout << "=================================" << endl;

  {
    cout << "Real matrix" << endl;
    mat X = randn(5, 5);
    mat L, U;
    ivec p;
    lu(X, L, U, p);

    mat P = to_mat(permutation_matrix(p));
    cout << "X = " << round_to_zero(X) << endl;
    cout << "norm(X - P^T*L*U) = "
         << round_to_zero(norm(X - transpose(P) * L * U)) << endl;
  }
  {
    cout << "Complex matrix" << endl;
    cmat X = randn_c(5, 5);
    cmat L, U;
    ivec p;
    lu(X, L, U, p);

    mat P = to_mat(permutation_matrix(p));
    cout << "X = " << round_to_zero(X) << endl;
    cout << "norm(X - P^H*L*U = "
         << round_to_zero(norm(X - hermitian_transpose(P) * L * U)) << endl;
  }

  return 0;
}
