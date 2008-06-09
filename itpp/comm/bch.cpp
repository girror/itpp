/*!
 * \file
 * \brief Implementation of a BCH encoder/decoder class
 * \author Pal Frenger, Steve Peters, Adam Piatyszek and Stephan Ludwig
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2008  (see AUTHORS file for a list of contributors)
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

#include <itpp/comm/bch.h>
#include <itpp/base/binary.h>
#include <itpp/base/specmat.h>
#include <itpp/base/array.h>

namespace itpp
{

//---------------------- BCH -----------------------------------

BCH::BCH(int in_n, int in_k, int in_t, const ivec &genpolynom, bool sys):
    n(in_n), k(in_k), t(in_t), systematic(sys)
{
  //fix the generator polynomial g(x).
  ivec exponents(n - k + 1);
  bvec temp = oct2bin(genpolynom);
  for (int i = 0; i < temp.length(); i++) {
    exponents(i) = static_cast<int>(temp(temp.length() - i - 1)) - 1;
  }
  g.set(n + 1, exponents);
}

BCH::BCH(int in_n, int in_k, int in_t, bool sys):
    n(in_n), k(in_k), t(in_t), systematic(sys)
{
  // step 1: determine cyclotomic cosets
  // although we use elements in GF(n+1), we do not use GFX class, but ivec,
  // since we have to multiply by 2 and need the exponents in clear notation
  int m_tmp = int2bits(n);
  int two_pow_m = 1 << m_tmp;

  it_assert(two_pow_m == n + 1, "BCH::BCH(): (in_n + 1) is not a power of 2");
  it_assert(t > 0, "BCH::BCH(): in_t must be positive");

  Array<ivec> cyclo_sets(2*t + 1);
  // unfortunately it is not obvious how many cyclotomic cosets exist (?)
  // a bad guess is n/2, which can be a whole lot...
  // but we only need 2*t + 1 at maximum for step 2.
  // (since all elements are sorted ascending [cp. comment at 2.], the last
  // coset we need is the one with coset leader 2t. + coset {0})

  // start with {0} as first set
  int curr_coset_idx = 0;
  cyclo_sets(curr_coset_idx) = zeros_i(1);

  int cycl_element = 1;

  do {
    bool found = false;
    // find next element, which is not in a previous coset
    do {
      int i = 0;
      // we do not have to search the first coset, since this is always {0}
      found = false;
      while ((!found) && (i <= curr_coset_idx)) {
        int j = 0;
        while ((!found) && (j < cyclo_sets(i).length())) {
          if (cycl_element == cyclo_sets(i)(j)) {
            found = true;
          }
          j++;
        }
        i++;
      }
      cycl_element++;
    }
    while ((found) && (cycl_element <= 2*t));

    if (!found) {
      // found one
      cyclo_sets(++curr_coset_idx).set_size(m_tmp);
      // a first guess (we delete afterwards back to correct length):
      // there should be no more than m elements in one coset

      int element_index = 0;
      cyclo_sets(curr_coset_idx)(element_index) = cycl_element - 1;

      // multiply by two (mod 2^m - 1) as long as new elements are created
      while ((((cyclo_sets(curr_coset_idx)(element_index) * 2) % n)
              != cyclo_sets(curr_coset_idx)(0))
             && (element_index < m_tmp - 1)) {
        element_index++;
        cyclo_sets(curr_coset_idx)(element_index)
          = (cyclo_sets(curr_coset_idx)(element_index - 1) * 2) % n;
      }
      // delete unused digits
      if (element_index + 1 < m_tmp - 1) {
        cyclo_sets(curr_coset_idx).del(element_index + 1, m_tmp - 1);
      }
    }
  }
  while ((cycl_element <= 2*t) && (curr_coset_idx <= 2*t));

  // step 2: find all cosets that contain all the powers (1..2t) of alpha
  // this is pretty easy, since the cosets are in ascending order
  // (if regarding the first (=primitive) element for ordering) -
  // all due to the method, they have been constructed
  // Since we only calculated all cosets up to 2t, this is even trivial
  // => we take all curr_coset_idx Cosets

  // maximum index of cosets to be considered
  int max_coset_index = curr_coset_idx;

  // step 3: multiply the minimal polynomials corresponding to this sets
  // of powers
  g.set(two_pow_m, ivec("0"));   // = alpha^0 = 1
  ivec min_poly_exp(2);
  min_poly_exp(1) = 0;        // product of (x-alpha^cycl_element)

  for (int i = 1; i <= max_coset_index; i++) {
    for (int j = 0; j < cyclo_sets(i).length(); j++) {
      min_poly_exp(0) = cyclo_sets(i)(j);
      g *= GFX(two_pow_m, min_poly_exp);
    }
  }

  // finally check, whether k and t match, i.e. if degree(gen_poly) = n - k
  it_assert(g.get_true_degree() == n - k,
            "BCH::BCH(): parameters n, k, t do not match each other");
}


void BCH::encode(const bvec &uncoded_bits, bvec &coded_bits)
{
  int i, j, degree,
  itterations = floor_i(static_cast<double>(uncoded_bits.length()) / k);
  GFX m(n + 1, k);
  GFX c(n + 1, n);
  GFX r(n + 1, n - k);
  GFX uncoded_shifted(n + 1, n);
  coded_bits.set_size(itterations*n, false);
  bvec mbit(k), cbit(n);

  if (systematic)
    for (i = 0; i < n - k; i++)
      uncoded_shifted[i] = GF(n + 1, -1);

  for (i = 0; i < itterations; i++) {
    //Fix the message polynom m(x).
    mbit = uncoded_bits.mid(i * k, k);
    for (j = 0; j < k; j++) {
      degree = static_cast<int>(mbit(j)) - 1;
      m[j] = GF(n + 1, degree);
      if (systematic) {
        c[j] = m[j];
        uncoded_shifted[j+n-k] = m[j];
      }
    }
    //Fix the outputbits cbit.
    if (systematic) {
      r = modgfx(uncoded_shifted, g);
      for (j = k; j < n; j++) {
        c[j] = r[j-k];
      }
    }
    else {
      c = g * m;
    }
    for (j = 0; j < n; j++) {
      if (c[j] == GF(n + 1, 0)) {
        cbit(j) = 1;
      }
      else {
        cbit(j) = 0;
      }
    }
    coded_bits.replace_mid(i*n, cbit);
  }
}

bvec BCH::encode(const bvec &uncoded_bits)
{
  bvec coded_bits;
  encode(uncoded_bits, coded_bits);
  return coded_bits;
}

void BCH::decode(const bvec &coded_bits, bvec &decoded_bits)
{
  int j, i, degree, kk, foundzeros, cisvalid,
  itterations = floor_i(static_cast<double>(coded_bits.length()) / n);
  bvec rbin(n), mbin(k);
  decoded_bits.set_size(itterations*k, false);

  GFX r(n + 1, n - 1), c(n + 1, n - 1), m(n + 1, k - 1), S(n + 1, 2*t), Lambda(n + 1),
  OldLambda(n + 1), T(n + 1), Ohmega(n + 1), One(n + 1, (char*)"0");
  GF delta(n + 1), temp(n + 1);
  ivec errorpos;

  for (i = 0; i < itterations; i++) {
    //Fix the received polynomial r(x)
    rbin = coded_bits.mid(i * n, n);
    for (j = 0; j < n; j++) {
      degree = static_cast<int>(rbin(j)) - 1;
      r[j] = GF(n + 1, degree);
    }
    //Fix the syndrome polynomial S(x).
    S[0] = GF(n + 1, -1);
    for (j = 1; j <= 2*t; j++) {
      S[j] =  r(GF(n + 1, j));
    }
    if (S.get_true_degree() >= 1) { //Errors in the received word
      //Itterate to find Lambda(x).
      kk = 0;
      Lambda = GFX(n + 1, (char*)"0");
      T = GFX(n + 1, (char*)"0");
      while (kk < t) {
        Ohmega = Lambda * (S + One);
        delta = Ohmega[2*kk+1];
        OldLambda = Lambda;
        Lambda = OldLambda + delta * (GFX(n + 1, (char*)"-1 0") * T);
        if ((delta == GF(n + 1, -1)) || (OldLambda.get_degree() > kk)) {
          T = GFX(n + 1, (char*)"-1 -1 0") * T;
        }
        else {
          T = (GFX(n + 1, (char*)"-1 0") * OldLambda) / delta;
        }
        kk = kk + 1;
      }
      //Find the zeros to Lambda(x).
      errorpos.set_size(Lambda.get_true_degree(), true);
      foundzeros = 0;
      for (j = 0; j <= n - 1; j++) {
        temp = Lambda(GF(n + 1, j));
        if (Lambda(GF(n + 1, j)) == GF(n + 1, -1)) {
          errorpos(foundzeros) = (n - j) % n;
          foundzeros += 1;
          if (foundzeros >= Lambda.get_true_degree()) {
            break;
          }
        }
      }
      //Correct the codeword.
      for (j = 0; j < foundzeros; j++) {
        rbin(errorpos(j)) += 1;
      }
      //Reconstruct the corrected codeword.
      for (j = 0; j < n; j++) {
        degree = static_cast<int>(rbin(j)) - 1;
        c[j] = GF(n + 1, degree);
      }
      //Code word validation.
      S[0] = GF(n + 1, -1);
      for (j = 1; j <= 2*t; j++) {
        S[j] =  c(GF(n + 1, j));
      }
      if (S.get_true_degree() <= 0) { //c(x) is a valid codeword.
        cisvalid = true;
      }
      else {
        cisvalid = false;
      }
    }
    else {
      c = r;
      cisvalid = true;
    }
    //Construct the message bit vector.
    if (cisvalid) { //c(x) is a valid codeword.
      if (c.get_true_degree() > 1) {
        if (systematic) {
          for (j = 0; j < k; j++)
            m[j] = c[j];
        }
        else {
          m = divgfx(c, g);
        }
        mbin.clear();
        for (j = 0; j <= m.get_true_degree(); j++) {
          if (m[j] == GF(n + 1, 0)) {
            mbin(j) = 1;
          }
        }
      }
      else { //The zero word was transmitted
        mbin = zeros_b(k);
        m = GFX(n + 1, (char*)"-1");
      }
    }
    else { //Decoder failure.
      mbin = zeros_b(k);
      m = GFX(n + 1, (char*)"-1");
    }
    decoded_bits.replace_mid(i*k, mbin);
  }
}


bvec BCH::decode(const bvec &coded_bits)
{
  bvec decoded_bits;
  decode(coded_bits, decoded_bits);
  return decoded_bits;
}


// --- Soft-decision decoding is not implemented ---

void BCH::decode(const vec &, bvec &)
{
  it_error("BCH::decode(): Soft-decision decoding is not implemented");
}

bvec BCH::decode(const vec &)
{
  it_error("BCH::decode(): Soft-decision decoding is not implemented");
  return bvec();
}


} // namespace itpp
