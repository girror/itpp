/*!
* \file 
* \brief Transforms test program
* \author Tony Ottosson, Thomas Eriksson, Simon Wood and Adam Piatyszek
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

#include <itpp/itbase.h>

using namespace itpp;
using namespace std;


#if defined(HAVE_FFT)

int main()
{
 cout << "========================" << endl;
 cout << "   Test of Transforms   " << endl;
 cout << "========================" << endl << endl;

 int N;
 {
   vec x, z;
   cvec y;

   N = 16;
   x = randn(N);
   cout << "Test 1: FFT/IFFT; Real input vector; N = " << N 
	<< ";  fft_real(x, y), ifft_real(y, z):" << endl << endl;
   
   cout << "x = " << round_to_zero(x) << endl;
   fft_real(x, y);
   cout << "y = " << round_to_zero(y) << endl;
   ifft_real(y, z);
   cout << "z = " << round_to_zero(z) << endl << endl;

   N = 11;
   x = randn(N);
   cout << "Test 2: FFT/IFFT; Real input vector; N = " << N 
	<< ";  fft_real(x, y), ifft_real(y, z):" << endl << endl;

   cout << "x = " << round_to_zero(x) << endl;
   fft_real(x, y);
   cout << "y = " << round_to_zero(y) << endl;
   ifft_real(y, z);
   cout << "z = " << round_to_zero(z) << endl << endl;
 }
 {
   cvec x, y, z;

   N = 32;
   x = randn_c(N);
   cout << "Test 3: FFT/IFFT; Complex input vector; N = " << N 
	<< ";  fft(x, y), ifft(y, z):" << endl << endl;
   
   cout << "x = " << round_to_zero(x) << endl;
   fft(x, y);
   cout << "y = " << round_to_zero(y) << endl;
   ifft(y, z);
   cout << "z = " << round_to_zero(z) << endl << endl;

   N = 7;
   x = randn_c(N);
   cout << "Test 4: FFT/IFFT; Complex input vector; N = " << N 
	<< ";  fft(x, y), ifft(y, z):" << endl << endl;
   
   cout << "x = " << round_to_zero(x) << endl;
   fft(x, y);
   cout << "y = " << round_to_zero(y) << endl;
   ifft(y, z);
   cout << "z = " << round_to_zero(z) << endl << endl;
 }
 {
   vec x, y, z;

   N = 8;
   x = randn(N);
   cout << "Test 5: DCT/IDCT; Real input vector; N = " << N 
	<< ";  dct(x, y), idct(y, z):" << endl << endl;
   
   cout << "x = " << round_to_zero(x) << endl;
   dct(x, y);
   cout << "y = " << round_to_zero(y) << endl;
   idct(y, z);
   cout << "z = " << round_to_zero(z) << endl << endl;

   N = 11;
   x = randn(N);
   cout << "Test 6: DCT/IDCT; Real input vector; N = " << N 
	<< ";  dct(x, y), idct(y, z):" << endl << endl;
   
   cout << "x = " << round_to_zero(x) << endl;
   dct(x, y);
   cout << "y = " << round_to_zero(y) << endl;
   idct(y, z);
   cout << "z = " << round_to_zero(z) << endl;
 }

 return 0;
}

#else

int main() { 
 cerr << "Error: FFTW library is needed to run this test program" << endl; 
 return 1;
}

#endif
