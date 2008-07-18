/*!
 * \file
 * \brief Templated Vector Class Implementation
 * \author Tony Ottosson, Tobias Ringstrom, Adam Piatyszek and Conrad Sanderson
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

#include <itpp/base/vec.h>
#include <itpp/base/converters.h>
#include <limits>

//! \cond

namespace itpp
{


template<class Num_T>
std::string Vec<Num_T>::replace_commas(const std::string &str_in)
{
  // copy an input sting into a local variable str
  std::string str(str_in);
  // find first occurence of comma in string str
  std::string::size_type index = str.find(',', 0);
  while (index != std::string::npos) {
    // replace character at position index with space
    str.replace(index, 1, 1, ' ');
    // find next occurence of comma in string str
    index = str.find(',', index);
  }
  return str;
}


template<>
void Vec<double>::set(const std::string &str)
{
  std::istringstream buffer(replace_commas(str));
  double b = 0.0;
  double c = 0.0;
  double eps_margin;
  bool b_parsed = false;
  bool c_parsed = false;
  bool negative = false;
  bool nan_inf = false;
  int pos = 0, maxpos = 10;

  free();
  alloc(maxpos);

  while (buffer.peek() != EOF) {
    switch (buffer.peek()) {
      // skip spaces
    case ' ':
    case '\t':
      buffer.seekg(1, std::ios_base::cur);
      break;

      // skip '+' sign
    case '+':
      // check for not handled '-' sign
      it_assert(!negative, "Vec<double>::set(): Improper data string (-)");
      buffer.seekg(1, std::ios_base::cur);
      break;

      // check for '-' sign
    case '-':
      buffer.seekg(1, std::ios_base::cur);
      negative = true;
      break;

      // check for NaN
    case 'N':
    case 'n':
      buffer.seekg(1, std::ios_base::cur);
      it_assert((buffer.peek() == 'A') || (buffer.peek() == 'a'),
                "Vec<double>::set(): Improper data string (NaN)");
      buffer.seekg(1, std::ios_base::cur);
      it_assert((buffer.peek() == 'N') || (buffer.peek() == 'n'),
                "Vec<double>::set(): Improper data string (NaN)");
      buffer.seekg(1, std::ios_base::cur);
      it_assert(!negative, "Vec<double>::set(): Improper data string "
                "(-NaN not exist)");
      if (++pos > maxpos) {
        maxpos <<= 1;
        set_size(maxpos, true);
      }
      if (std::numeric_limits<double>::has_quiet_NaN) {
        data[pos-1] = std::numeric_limits<double>::quiet_NaN();
      }
      else if (std::numeric_limits<double>::has_signaling_NaN) {
        data[pos-1] = std::numeric_limits<double>::signaling_NaN();
      }
      else {
        it_error("Vec<double::set(): NaN not supported");
      }
      nan_inf = true;
      break; // case 'N'...

      // check for Inf
    case 'I':
    case 'i':
      buffer.seekg(1, std::ios_base::cur);
      it_assert((buffer.peek() == 'N') || (buffer.peek() == 'n'),
                "Vec<double>::set(): Improper data string (Inf)");
      buffer.seekg(1, std::ios_base::cur);
      it_assert((buffer.peek() == 'F') || (buffer.peek() == 'f'),
                "Vec<double>::set(): Improper data string (Inf)");
      buffer.seekg(1, std::ios_base::cur);
      it_assert(std::numeric_limits<double>::has_infinity,
                "Vec<double::set(): Inf not supported");
      if (++pos > maxpos) {
        maxpos <<= 1;
        set_size(maxpos, true);
      }
      if (negative) {
        data[pos-1] = -std::numeric_limits<double>::infinity();
        negative = false;
      }
      else {
        data[pos-1] = std::numeric_limits<double>::infinity();
      }
      nan_inf = true;
      break; // case 'I'...

    case ':': // reads format a:b:c or a:b
      it_assert(!negative, "Vec<double>::set(): Improper data string (-)");
      it_assert(!nan_inf, "Vec<double>::set(): Improper data string (Nan/Inf "
                " can not be used with a:b or a:b:c)");
      it_assert(pos == 1, "Vec<double>::set(): Improper data string (a:b)");
      buffer.seekg(1, std::ios_base::cur);
      // parse b
      while (buffer.peek() != EOF) {
        switch (buffer.peek()) {
        case ' ':
        case '\t':
          buffer.seekg(1, std::ios_base::cur);
          break;

        case ':':
          it_assert(b_parsed, "Vec<double>::set(): Improper data string "
                    "(a:b)");
          buffer.seekg(1, std::ios_base::cur);
          // parse c
          while (buffer.peek() != EOF) {
            switch (buffer.peek()) {
            case ' ':
            case '\t':
              buffer.seekg(1, std::ios_base::cur);
              break;

            default:
              it_assert(!c_parsed, "Vec<double>::set(): Improper data "
                        "string (a:b:c)");
              buffer.clear();
              buffer >> c;
              it_assert(!buffer.fail(), "Vec<double>::set(): Stream "
                        "operation failed (buffer >> c)");
              c_parsed = true;
            }
          }
          it_assert(c_parsed, "Vec<double>::set(): Improper data string "
                    "(a:b:c)");
          break;

        default:
          it_assert(!b_parsed, "Vec<double>::set(): Improper data string "
                    "(a:b)");
          buffer.clear();
          buffer >> b;
          it_assert(!buffer.fail(), "Vec<double>::set(): Stream operation "
                    "failed (buffer >> b)");
          b_parsed = true;
        }
      }
      it_assert(b_parsed, "Vec<double>::set(): Improper data string (a:b)");

      if (c_parsed) {
        // Adding this margin fixes precision problems in e.g. "0:0.2:3",
        // where the last value was 2.8 instead of 3.
        eps_margin = std::fabs((c - data[pos-1]) / b) * eps;
        if (b > 0 && c >= data[pos-1]) {
          while (data[pos-1] + b <= c + eps_margin) {
            if (++pos > maxpos) {
              maxpos <<= 1;
              set_size(maxpos, true);
            }
            data[pos-1] = data[pos-2] + b;
          }
        }
        else if (b < 0 && c <= data[pos-1]) {
          while (data[pos-1] + b >= c - eps_margin) {
            if (++pos > maxpos) {
              maxpos <<= 1;
              set_size(maxpos, true);
            }
            data[pos-1] = data[pos-2] + b;
          }
        }
        else if (b == 0 && c == data[pos-1]) {
          break;
        }
        else {
          it_error("Vec<double>::set(): Improper data string (a:b:c)");
        }
      } // if (c_parsed)
      else if (b_parsed) {
        eps_margin = std::fabs(b - data[pos-1]) * eps;
        if (b < data[pos-1]) {
          while (data[pos-1] - 1.0 >= b - eps_margin) {
            if (++pos > maxpos) {
              maxpos <<= 1;
              set_size(maxpos, true);
            }
            data[pos-1] = data[pos-2] - 1.0;
          }
        }
        else {
          while (data[pos-1] + 1.0 <= b + eps_margin) {
            if (++pos > maxpos) {
              maxpos <<= 1;
              set_size(maxpos, true);
            }
            data[pos-1] = data[pos-2] + 1.0;
          }
        }
      } // else if (b_parsed)
      else {
        it_error("Vec<double>::set(): Improper data string (a:b)");
      }
      break; // case ':'

    default:
      if (++pos > maxpos) {
        maxpos <<= 1;
        set_size(maxpos, true);
      }
      buffer >> data[pos-1];
      it_assert(!buffer.fail(), "Vec<double>::set(): Stream operation "
                "failed (buffer >> data)");
      if (negative) {
        data[pos-1] = -data[pos-1];
        negative = false;
      }
      break; // default
    }
  }
  set_size(pos, true);
}


template<>
void Vec<std::complex<double> >::set(const std::string &str)
{
  std::istringstream buffer(str);
  int pos = 0, maxpos = 10;

  free();
  alloc(maxpos);

  while (buffer.peek() != EOF) {
    switch (buffer.peek()) {
    case ':':
      it_error("Vec<complex>::set(): a:b:c and a:b expressions not valid "
               "for cvec");
      break;
    case ' ':
    case '\t':
    case ',':
      buffer.seekg(1, std::ios_base::cur);
      break;
    default:
      if (++pos > maxpos) {
        maxpos <<= 1;
        set_size(maxpos, true);
      }
      buffer >> data[pos-1];
      it_assert(!buffer.fail(), "Vec<complex>::set(): Stream operation "
                "failed (buffer >> data)");
    }
  }
  set_size(pos, true);
}


template<>
void Vec<bin>::set(const std::string &str)
{
  std::istringstream buffer(replace_commas(str));
  int pos = 0, maxpos = 10;

  free();
  alloc(maxpos);

  while (buffer.peek() != EOF) {
    switch (buffer.peek()) {
    case ':':
      it_error("Vec<bin>::set(): a:b:c and a:b expressions not valid "
               "for bvec");
      break;
    case ' ':
    case '\t':
      buffer.seekg(1, std::ios_base::cur);
      break;
    default:
      if (++pos > maxpos) {
        maxpos <<= 1;
        set_size(maxpos, true);
      }
      buffer >> data[pos-1];
      it_assert(!buffer.fail(), "Vec<bin>::set(): Stream operation failed "
                "(buffer >> data)");
    }
  }
  set_size(pos, true);
}


template<>
void Vec<int>::set(const std::string &str)
{
  std::istringstream buffer(replace_commas(str));
  int b = 0;
  int c = 0;
  bool b_parsed = false;
  bool c_parsed = false;
  bool negative = false;
  int pos = 0;
  int maxpos = 10;

  free();
  alloc(maxpos);

  while (buffer.peek() != EOF) {
    switch (buffer.peek()) {
      // skip spaces and tabs
    case ' ':
    case '\t':
      buffer.seekg(1, std::ios_base::cur);
      break;

      // skip '+' sign
    case '+':
      // check for not handled '-' sign
      it_assert(!negative, "Vec<double>::set(): Improper data string (-)");
      buffer.seekg(1, std::ios_base::cur);
      break;

      // check for '-' sign
    case '-':
      buffer.seekg(1, std::ios_base::cur);
      negative = true;
      break;

      // hexadecimal number or octal number or zero
    case '0':
      buffer.seekg(1, std::ios_base::cur);
      switch (buffer.peek()) {
        // hexadecimal number
      case 'x':
      case 'X':
        buffer.clear();
        buffer.seekg(-1, std::ios_base::cur);
        if (++pos > maxpos) {
          maxpos <<= 1;
          set_size(maxpos, true);
        }
        buffer >> std::hex >> data[pos-1];
        it_assert(!buffer.fail(), "Vec<int>::set(): Stream operation "
                  "failed (buffer >> hex >> data)");
        break; // case 'x'...

        // octal number
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
        buffer.clear();
        buffer.seekg(-1, std::ios_base::cur);
        if (++pos > maxpos) {
          maxpos <<= 1;
          set_size(maxpos, true);
        }
        buffer >> std::oct >> data[pos-1];
        it_assert(!buffer.fail(), "Vec<int>::set(): Stream operation "
                  "failed (buffer >> oct >> data)");
        break; // case '1'...

        // zero
      case EOF:
      case ' ':
      case '\t':
      case ':':
      case '0':
        buffer.clear();
        buffer.seekg(-1, std::ios_base::cur);
        if (++pos > maxpos) {
          maxpos <<= 1;
          set_size(maxpos, true);
        }
        buffer >> std::dec >> data[pos-1];
        it_assert(!buffer.fail(), "Vec<int>::set(): Stream operation "
                  "failed (buffer >> dec >> data)");
        break; // case EOF...

      default:
        it_error("Vec<int>::set(): Improper data string");
      }
      // check if just parsed data was negative
      if (negative) {
        data[pos-1] = -data[pos-1];
        negative = false;
      }
      break; // case '0'

      // decimal number
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      buffer.clear();
      if (++pos > maxpos) {
        maxpos <<= 1;
        set_size(maxpos, true);
      }
      buffer >> std::dec >> data[pos-1];
      it_assert(!buffer.fail(), "Vec<int>::set(): Stream operation "
                "failed (buffer >> dec >> data)");
      // check if just parsed data was negative
      if (negative) {
        data[pos-1] = -data[pos-1];
        negative = false;
      }
      break; // case '1'...

      // parse format a:b:c or a:b
    case ':':
      it_assert(pos == 1, "Vec<int>::set(): Improper data string (a:b)");
      buffer.seekg(1, std::ios_base::cur);
      // parse b
      while (buffer.peek() != EOF) {
        switch (buffer.peek()) {
        case ' ':
        case '\t':
          buffer.seekg(1, std::ios_base::cur);
          break;

          // skip '+' sign
        case '+':
          // check for not handled '-' sign
          it_assert(!negative, "Vec<double>::set(): Improper data string "
                    "(-)");
          buffer.seekg(1, std::ios_base::cur);
          break;

          // check for '-' sign
        case '-':
          buffer.seekg(1, std::ios_base::cur);
          negative = true;
          break;

          // hexadecimal number or octal number or zero
        case '0':
          it_assert(!b_parsed, "Vec<int>::set(): Improper data string "
                    "(a:b)");
          buffer.seekg(1, std::ios_base::cur);
          switch (buffer.peek()) {
            // hexadecimal number
          case 'x':
          case 'X':
            buffer.clear();
            buffer.seekg(-1, std::ios_base::cur);
            buffer >> std::hex >> b;
            it_assert(!buffer.fail(), "Vec<int>::set(): Stream operation "
                      "failed (buffer >> hex >> data)");
            break; // case 'x'...

            // octal number
          case '1':
          case '2':
          case '3':
          case '4':
          case '5':
          case '6':
          case '7':
            buffer.clear();
            buffer.seekg(-1, std::ios_base::cur);
            buffer >> std::oct >> b;
            it_assert(!buffer.fail(), "Vec<int>::set(): Stream operation "
                      "failed (buffer >> oct >> data)");
            break; // case '1'...

            // zero
          case EOF:
          case ' ':
          case '\t':
          case ':':
          case '0':
            buffer.clear();
            buffer.seekg(-1, std::ios_base::cur);
            buffer >> std::dec >> b;
            it_assert(!buffer.fail(), "Vec<int>::set(): Stream operation "
                      "failed (buffer >> dec >> data)");
            break; // case EOF...

          default:
            it_error("Vec<int>::set(): Improper data string (a:b)");
          } // switch (buffer.peek())
          // check if just parsed data was negative
          if (negative) {
            b = -b;
            negative = false;
          }
          b_parsed = true;
          break; // case '0'

          // decimal number
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          it_assert(!b_parsed, "Vec<int>::set(): Improper data string "
                    "(a:b)");
          buffer.clear();
          buffer >> std::dec >> b;
          it_assert(!buffer.fail(), "Vec<int>::set(): Stream operation "
                    "failed (buffer >> dec >> data)");
          // check if just parsed data was negative
          if (negative) {
            b = -b;
            negative = false;
          }
          b_parsed = true;
          break; // case '1'...

        case ':':
          it_assert(b_parsed, "Vec<int>::set(): Improper data string (a:b)");
          buffer.seekg(1, std::ios_base::cur);
          // parse c
          while (buffer.peek() != EOF) {
            switch (buffer.peek()) {
            case ' ':
            case '\t':
              buffer.seekg(1, std::ios_base::cur);
              break;

              // skip '+' sign
            case '+':
              // check for not handled '-' sign
              it_assert(!negative, "Vec<double>::set(): Improper data "
                        "string (-)");
              buffer.seekg(1, std::ios_base::cur);
              break;

              // check for '-' sign
            case '-':
              buffer.seekg(1, std::ios_base::cur);
              negative = true;
              break;

              // hexadecimal number or octal number or zero
            case '0':
              it_assert(!c_parsed, "Vec<int>::set(): Improper data string "
                        "(a:b:c)");
              buffer.seekg(1, std::ios_base::cur);
              switch (buffer.peek()) {
                // hexadecimal number
              case 'x':
              case 'X':
                buffer.clear();
                buffer.seekg(-1, std::ios_base::cur);
                buffer >> std::hex >> c;
                it_assert(!buffer.fail(), "Vec<int>::set(): Stream "
                          "operation failed (buffer >> hex >> data)");
                break; // case 'x'...

                // octal number
              case '1':
              case '2':
              case '3':
              case '4':
              case '5':
              case '6':
              case '7':
                buffer.clear();
                buffer.seekg(-1, std::ios_base::cur);
                buffer >> std::oct >> c;
                it_assert(!buffer.fail(), "Vec<int>::set(): Stream "
                          "operation failed (buffer >> oct >> data)");
                break; // case '1'...

                // zero
              case EOF:
              case ' ':
              case '\t':
              case '0':
                buffer.clear();
                buffer.seekg(-1, std::ios_base::cur);
                buffer >> std::dec >> c;
                it_assert(!buffer.fail(), "Vec<int>::set(): Stream "
                          "operation failed (buffer >> dec >> data)");
                break; // case EOF...

              default:
                it_error("Vec<int>::set(): Improper data string (a:b:c)");
              }
              c_parsed = true;
              break; // case '0'

              // decimal number
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
              it_assert(!c_parsed, "Vec<int>::set(): Improper data string "
                        "(a:b:c)");
              buffer.clear();
              buffer >> std::dec >> c;
              it_assert(!buffer.fail(), "Vec<int>::set(): Stream operation "
                        "failed (buffer >> dec >> data)");
              c_parsed = true;
              break;

            default:
              it_error("Vec<int>::set(): Improper data string (a:b:c)");
            } // switch (buffer.peek())
          } // while (buffer.peek() != EOF)
          // check if just parsed data was negative
          if (negative) {
            c = -c;
            negative = false;
          }
          it_assert(c_parsed, "Vec<int>::set(): Improper data string "
                    "(a:b:c)");
          break; // case ':'

        default:
          it_error("Vec<int>::set(): Improper data string (a:b)");
        } // switch (buffer.peek())
      } // while (buffer.peek() != EOF)

      if (c_parsed) {
        if (b > 0 && c >= data[pos-1]) {
          while (data[pos-1] + b <= c) {
            if (++pos > maxpos) {
              maxpos <<= 1;
              set_size(maxpos, true);
            }
            data[pos-1] = data[pos-2] + b;
          }
        }
        else if (b < 0 && c <= data[pos-1]) {
          while (data[pos-1] + b >= c) {
            if (++pos > maxpos) {
              maxpos <<= 1;
              set_size(maxpos, true);
            }
            data[pos-1] = data[pos-2] + b;
          }
        }
        else if (b == 0 && c == data[pos-1]) {
          break;
        }
        else {
          it_error("Vec<int>::set(): Improper data string (a:b:c)");
        }
      } // if (c_parsed)
      else if (b_parsed) {
        if (b < data[pos-1]) {
          while (data[pos-1] > b) {
            if (++pos > maxpos) {
              maxpos <<= 1;
              set_size(maxpos, true);
            }
            data[pos-1] = data[pos-2] - 1;
          }
        }
        else {
          while (data[pos-1] < b) {
            if (++pos > maxpos) {
              maxpos <<= 1;
              set_size(maxpos, true);
            }
            data[pos-1] = data[pos-2] + 1;
          }
        }
      } // else if (b_parsed)
      else {
        it_error("Vec<int>::set(): Improper data string (a:b)");
      }
      break; // case ':'

    default:
      it_error("Vec<int>::set(): Improper data string");
    }
  }
  // resize the parsed vector to its final length
  set_size(pos, true);
}

template<>
void Vec<short int>::set(const std::string &str)
{
  // parser for "short int" is the same as for "int", so reuse it here
  ivec iv(str);
  this->operator=(to_svec(iv));
}


template<>
bvec Vec<std::complex<double> >::operator==(std::complex<double>) const
{
  it_error("operator==: not implemented for complex");
  bvec temp;
  return temp;
}

template<>
bvec Vec<std::complex<double> >::operator!=(std::complex<double>) const
{
  it_error("operator!=: not implemented for complex");
  bvec temp;
  return temp;
}

template<>
bvec Vec<std::complex<double> >::operator<=(std::complex<double>) const
{
  it_error("operator<=: not implemented for complex");
  bvec temp;
  return temp;
}

template<>
bvec Vec<std::complex<double> >::operator>(std::complex<double>) const
{
  it_error("operator>: not implemented for complex");
  bvec temp;
  return temp;
}

template<>
bvec Vec<std::complex<double> >::operator<(std::complex<double>) const
{
  it_error("operator<: not implemented for complex");
  bvec temp;
  return temp;
}

template<>
bvec Vec<std::complex<double> >::operator>=(std::complex<double>) const
{
  it_error("operator>=: not implemented for complex");
  bvec temp;
  return temp;
}

template<>
Mat<std::complex<double> > Vec<std::complex<double> >::hermitian_transpose() const
{
  Mat<std::complex<double> > temp(1, datasize);
  for (int i = 0; i < datasize; i++)
    temp(i) = std::conj(data[i]);

  return temp;
}


//---------------------------------------------------------------------
// Instantiations
//---------------------------------------------------------------------

template class Vec<double>;
template class Vec<int>;
template class Vec<short int>;
template class Vec<std::complex<double> >;
template class Vec<bin>;

// addition operator

template vec operator+(const vec &v1, const vec &v2);
template cvec operator+(const cvec &v1, const cvec &v2);
template ivec operator+(const ivec &v1, const ivec &v2);
template svec operator+(const svec &v1, const svec &v2);
template bvec operator+(const bvec &v1, const bvec &v2);

template vec operator+(const vec &v1, double t);
template cvec operator+(const cvec &v1, std::complex<double> t);
template ivec operator+(const ivec &v1, int t);
template svec operator+(const svec &v1, short t);
template bvec operator+(const bvec &v1, bin t);

template vec operator+(double t, const vec &v1);
template cvec operator+(std::complex<double> t, const cvec &v1);
template ivec operator+(int t, const ivec &v1);
template svec operator+(short t, const svec &v1);
template bvec operator+(bin t, const bvec &v1);

// subraction operator

template vec operator-(const vec &v1, const vec &v2);
template cvec operator-(const cvec &v1, const cvec &v2);
template ivec operator-(const ivec &v1, const ivec &v2);
template svec operator-(const svec &v1, const svec &v2);
template bvec operator-(const bvec &v1, const bvec &v2);

template vec operator-(const vec &v, double t);
template cvec operator-(const cvec &v, std::complex<double> t);
template ivec operator-(const ivec &v, int t);
template svec operator-(const svec &v, short t);
template bvec operator-(const bvec &v, bin t);

template vec operator-(double t, const vec &v);
template cvec operator-(std::complex<double> t, const cvec &v);
template ivec operator-(int t, const ivec &v);
template svec operator-(short t, const svec &v);
template bvec operator-(bin t, const bvec &v);

// unary minus

template vec operator-(const vec &v);
template cvec operator-(const cvec &v);
template ivec operator-(const ivec &v);
template svec operator-(const svec &v);
template bvec operator-(const bvec &v);

// multiplication operator

#if !defined(HAVE_BLAS)
template double dot(const vec &v1, const vec &v2);
#if !(defined(HAVE_ZDOTUSUB) || defined(HAVE_ZDOTU_VOID))
template std::complex<double> dot(const cvec &v1, const cvec &v2);
#endif // !(HAVE_ZDOTUSUB || HAVE_ZDOTU_VOID)
#endif // HAVE_BLAS
template int dot(const ivec &v1, const ivec &v2);
template short dot(const svec &v1, const svec &v2);
template bin dot(const bvec &v1, const bvec &v2);

#if !defined(HAVE_BLAS)
template double operator*(const vec &v1, const vec &v2);
template std::complex<double> operator*(const cvec &v1, const cvec &v2);
#endif
template int operator*(const ivec &v1, const ivec &v2);
template short operator*(const svec &v1, const svec &v2);
template bin operator*(const bvec &v1, const bvec &v2);

#if !defined(HAVE_BLAS)
template mat outer_product(const vec &v1, const vec &v2, bool hermitian);
#endif
template imat outer_product(const ivec &v1, const ivec &v2, bool hermitian);
template smat outer_product(const svec &v1, const svec &v2, bool hermitian);
template bmat outer_product(const bvec &v1, const bvec &v2, bool hermitian);

template vec operator*(const vec &v, double t);
template cvec operator*(const cvec &v, std::complex<double> t);
template ivec operator*(const ivec &v, int t);
template svec operator*(const svec &v, short t);
template bvec operator*(const bvec &v, bin t);

template vec operator*(double t, const vec &v);
template cvec operator*(std::complex<double> t, const cvec &v);
template ivec operator*(int t, const ivec &v);
template svec operator*(short t, const svec &v);
template bvec operator*(bin t, const bvec &v);

// elementwise multiplication

template vec elem_mult(const vec &a, const vec &b);
template cvec elem_mult(const cvec &a, const cvec &b);
template ivec elem_mult(const ivec &a, const ivec &b);
template svec elem_mult(const svec &a, const svec &b);
template bvec elem_mult(const bvec &a, const bvec &b);

template void elem_mult_out(const vec &a, const vec &b, vec &out);
template void elem_mult_out(const cvec &a, const cvec &b, cvec &out);
template void elem_mult_out(const ivec &a, const ivec &b, ivec &out);
template void elem_mult_out(const svec &a, const svec &b, svec &out);
template void elem_mult_out(const bvec &a, const bvec &b, bvec &out);

template vec elem_mult(const vec &a, const vec &b, const vec &c);
template cvec elem_mult(const cvec &a, const cvec &b, const cvec &c);
template ivec elem_mult(const ivec &a, const ivec &b, const ivec &c);
template svec elem_mult(const svec &a, const svec &b, const svec &c);
template bvec elem_mult(const bvec &a, const bvec &b, const bvec &c);

template void elem_mult_out(const vec &a, const vec &b, const vec &c,
                            vec &out);
template void elem_mult_out(const cvec &a, const cvec &b, const cvec &c,
                            cvec &out);
template void elem_mult_out(const ivec &a, const ivec &b, const ivec &c,
                            ivec &out);
template void elem_mult_out(const svec &a, const svec &b, const svec &c,
                            svec &out);
template void elem_mult_out(const bvec &a, const bvec &b, const bvec &c,
                            bvec &out);

template vec elem_mult(const vec &a, const vec &b, const vec &c,
                       const vec &d);
template cvec elem_mult(const cvec &a, const cvec &b, const cvec &c,
                        const cvec &d);
template ivec elem_mult(const ivec &a, const ivec &b, const ivec &c,
                        const ivec &d);
template svec elem_mult(const svec &a, const svec &b, const svec &c,
                        const svec &d);
template bvec elem_mult(const bvec &a, const bvec &b, const bvec &c,
                        const bvec &d);

template void elem_mult_out(const vec &a, const vec &b, const vec &c,
                            const vec &d, vec &out);
template void elem_mult_out(const cvec &a, const cvec &b, const cvec &c,
                            const cvec &d, cvec &out);
template void elem_mult_out(const ivec &a, const ivec &b, const ivec &c,
                            const ivec &d, ivec &out);
template void elem_mult_out(const svec &a, const svec &b, const svec &c,
                            const svec &d, svec &out);
template void elem_mult_out(const bvec &a, const bvec &b, const bvec &c,
                            const bvec &d, bvec &out);

// in-place elementwise multiplication

template void elem_mult_inplace(const vec &a, vec &b);
template void elem_mult_inplace(const cvec &a, cvec &b);
template void elem_mult_inplace(const ivec &a, ivec &b);
template void elem_mult_inplace(const svec &a, svec &b);
template void elem_mult_inplace(const bvec &a, bvec &b);

// elementwise multiplication followed by summation

template double elem_mult_sum(const vec &a, const vec &b);
template std::complex<double> elem_mult_sum(const cvec &a, const cvec &b);
template int elem_mult_sum(const ivec &a, const ivec &b);
template short elem_mult_sum(const svec &a, const svec &b);
template bin elem_mult_sum(const bvec &a, const bvec &b);

// division operator

template vec operator/(const vec &v, double t);
template cvec operator/(const cvec &v, std::complex<double> t);
template ivec operator/(const ivec &v, int t);
template svec operator/(const svec &v, short t);
template bvec operator/(const bvec &v, bin t);

template vec operator/(double t, const vec &v);
template cvec operator/(std::complex<double> t, const cvec &v);
template ivec operator/(int t, const ivec &v);
template svec operator/(short t, const svec &v);
template bvec operator/(bin t, const bvec &v);

// elementwise division operator

template vec elem_div(const vec &a, const vec &b);
template cvec elem_div(const cvec &a, const cvec &b);
template ivec elem_div(const ivec &a, const ivec &b);
template svec elem_div(const svec &a, const svec &b);
template bvec elem_div(const bvec &a, const bvec &b);

template vec elem_div(double t, const vec &v);
template cvec elem_div(std::complex<double> t, const cvec &v);
template ivec elem_div(int t, const ivec &v);
template svec elem_div(short t, const svec &v);
template bvec elem_div(bin t, const bvec &v);

template void elem_div_out(const vec &a, const vec &b, vec &out);
template void elem_div_out(const cvec &a, const cvec &b, cvec &out);
template void elem_div_out(const ivec &a, const ivec &b, ivec &out);
template void elem_div_out(const svec &a, const svec &b, svec &out);
template void elem_div_out(const bvec &a, const bvec &b, bvec &out);

// elementwise division followed by summation

template double elem_div_sum(const vec &a, const vec &b);
template std::complex<double> elem_div_sum(const cvec &a, const cvec &b);
template int elem_div_sum(const ivec &a, const ivec &b);
template short elem_div_sum(const svec &a, const svec &b);
template bin elem_div_sum(const bvec &a, const bvec &b);

// concat operator

template vec concat(const vec &v, double a);
template cvec concat(const cvec &v, std::complex<double> a);
template ivec concat(const ivec &v, int a);
template svec concat(const svec &v, short a);
template bvec concat(const bvec &v, bin a);

template vec concat(double a, const vec &v);
template cvec concat(std::complex<double> a, const cvec &v);
template ivec concat(int a, const ivec &v);
template svec concat(short a, const svec &v);
template bvec concat(bin a, const bvec &v);

template vec concat(const vec &v1, const vec &v2);
template cvec concat(const cvec &v1, const cvec &v2);
template ivec concat(const ivec &v1, const ivec &v2);
template svec concat(const svec &v1, const svec &v2);
template bvec concat(const bvec &v1, const bvec &v2);

template vec concat(const vec &v1, const vec &v2, const vec &v3);
template cvec concat(const cvec &v1, const cvec &v2, const cvec &v3);
template ivec concat(const ivec &v1, const ivec &v2, const ivec &v3);
template svec concat(const svec &v1, const svec &v2, const svec &v3);
template bvec concat(const bvec &v1, const bvec &v2, const bvec &v3);

template vec concat(const vec &v1, const vec &v2,
                    const vec &v3, const vec &v4);
template cvec concat(const cvec &v1, const cvec &v2,
                     const cvec &v3, const cvec &v4);
template ivec concat(const ivec &v1, const ivec &v2,
                     const ivec &v3, const ivec &v4);
template svec concat(const svec &v1, const svec &v2,
                     const svec &v3, const svec &v4);
template bvec concat(const bvec &v1, const bvec &v2,
                     const bvec &v3, const bvec &v4);

template vec concat(const vec &v1, const vec &v2, const vec &v3,
                    const vec &v4, const vec &v5);
template cvec concat(const cvec &v1, const cvec &v2, const cvec &v3,
                     const cvec &v4, const cvec &v5);
template ivec concat(const ivec &v1, const ivec &v2, const ivec &v3,
                     const ivec &v4, const ivec &v5);
template svec concat(const svec &v1, const svec &v2, const svec &v3,
                     const svec &v4, const svec &v5);
template bvec concat(const bvec &v1, const bvec &v2, const bvec &v3,
                     const bvec &v4, const bvec &v5);

// I/O streams

template std::ostream &operator<<(std::ostream& os, const vec &vect);
template std::ostream &operator<<(std::ostream& os, const cvec &vect);
template std::ostream &operator<<(std::ostream& os, const svec &vect);
template std::ostream &operator<<(std::ostream& os, const ivec &vect);
template std::ostream &operator<<(std::ostream& os, const bvec &vect);
template std::istream &operator>>(std::istream& is, vec &vect);
template std::istream &operator>>(std::istream& is, cvec &vect);
template std::istream &operator>>(std::istream& is, svec &vect);
template std::istream &operator>>(std::istream& is, ivec &vect);
template std::istream &operator>>(std::istream& is, bvec &vect);

} // namespace itpp

//! \endcond
