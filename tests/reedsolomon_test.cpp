/*!
 * \file
 * \brief Reed-Solomon encoder/decoder class test program
 * \author Steve Peters and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
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

#include <itpp/itcomm.h>

using namespace itpp;
using std::cout;
using std::endl;


int main()
{
  cout << "==========================================" << endl;
  cout << "   Test of Reed-Solomon encoder/decoder   " << endl;
  cout << "==========================================" << endl;

  bmat u = randb(8, 15);
  bmat c(u.rows(), 21);
  bmat y(u.rows(), 21);
  bvec codeword, errorword;
  bmat decoded(u.rows(), u.cols());
  Reed_Solomon rs(3, 1);
  Reed_Solomon rs_sys(3, 1, true);

  bmat f = "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0";

  cout << "Non-systematic case" << endl;
  cout << "-------------------" << endl;
  for (int i = 0; i < u.rows(); i++) {
    cout << "Info word:       " << u.get_row(i) << endl;
    codeword = rs.encode(u.get_row(i));
    it_assert(c.cols() == length(codeword), "Error 1");
    c.set_row(i, rs.encode(u.get_row(i)));
    cout << "Encoded:         " << c.get_row(i) << endl;
    errorword = f.get_row(i) + c.get_row(i);
    it_assert(y.cols() == length(errorword), "Error 2");
    y.set_row(i, f.get_row(i) + c.get_row(i));
    cout << "One error added: " << y.get_row(i) << endl;
    decoded.set_row(i, rs.decode(y.get_row(i)));
    cout << "Decoded to:      " << decoded.get_row(i) << endl << endl;
  }

  cout << "Systematic case" << endl;
  cout << "---------------" << endl;
  for (int i = 0; i < u.rows(); i++) {
    c.set_row(i, rs_sys.encode(u.get_row(i)));
    cout << "Info word:       " << u.get_row(i) << endl;
    cout << "Encoded:         " << c.get_row(i) << endl;
    y.set_row(i, f.get_row(i) + c.get_row(i));
    cout << "One error added: " << y.get_row(i) << endl;
    decoded.set_row(i, rs_sys.decode(y.get_row(i)));
    cout << "Decoded to:      " << decoded.get_row(i) << endl << endl;
  }

  return 0;
}
