/*!
 * \file
 * \brief Binary file formats implementations
 * \author Tony Ottosson, Thomas Eriksson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
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

#include <itpp/base/binfile.h>
#include <itpp/base/math/misc.h>
#include <cstring>


using std::ofstream;
using std::ifstream;
using std::fstream;
using std::ios;


namespace itpp
{

//! Read binary data and optionally switch endianness
template<typename T1, typename T2> inline
void read_endian(T1& st, T2& data, bool switch_endian = false)
{
  int bytes = sizeof(T2);
  char *c = reinterpret_cast<char *>(&data);
  if (!switch_endian)
    st.read(c, bytes);
  else
    for (int i = bytes - 1; i >= 0; i--)
      st.get(c[i]);
}

//! Write binary data and optionally switch endianness
template<typename T1, typename T2> inline
void write_endian(T1& st, T2 data, bool switch_endian = false)
{
  int bytes = sizeof(T2);
  char *c = reinterpret_cast<char *>(&data);
  if (!switch_endian)
    st.write(c, bytes);
  else
    for (int i = bytes - 1; i >= 0; i--)
      st.put(c[i]);
}

// ----------------------------------------------------------------------

bool exist(const std::string& name)
{
  bool file_exists = false;
  ifstream file(name.c_str(), ios::in);
  if (file.is_open()) {
    file_exists = true;
  }
  file.close();
  return file_exists;
}

// ----------------------------------------------------------------------
// bfstream_base
// ----------------------------------------------------------------------

bfstream_base::bfstream_base(endian e):
    switch_endianity(false),
    native_endianity(is_bigendian() ? b_endian : l_endian)
{
  if (native_endianity != e)
    switch_endianity = true;
}

// ----------------------------------------------------------------------
// bofstream
// ----------------------------------------------------------------------

bofstream::bofstream(const std::string& name, endian e) :
    bfstream_base(e), ofstream(name.c_str(), ios::out | ios::binary) {}

bofstream::bofstream() : bfstream_base(), ofstream() {}

void bofstream::open(const std::string& name, endian e)
{
  if (native_endianity != e)
    switch_endianity = true;
  else
    switch_endianity = false;
  ofstream::open(name.c_str(), ios::out | ios::binary);
}

bofstream& bofstream::operator<<(char a)
{
  put(a);
  return *this;
}

bofstream& bofstream::operator<<(unsigned char a)
{
  put(static_cast<char>(a));
  return *this;
}

bofstream& bofstream::operator<<(int16_t a)
{
  write_endian<bofstream, int16_t>(*this, a, switch_endianity);
  return *this;
}

bofstream& bofstream::operator<<(uint16_t a)
{
  write_endian<bofstream, uint16_t>(*this, a, switch_endianity);
  return *this;
}

bofstream& bofstream::operator<<(int32_t a)
{
  write_endian<bofstream, int32_t>(*this, a, switch_endianity);
  return *this;
}

bofstream& bofstream::operator<<(uint32_t a)
{
  write_endian<bofstream, uint32_t>(*this, a, switch_endianity);
  return *this;
}

bofstream& bofstream::operator<<(int64_t a)
{
  write_endian<bofstream, int64_t>(*this, a, switch_endianity);
  return *this;
}

bofstream& bofstream::operator<<(uint64_t a)
{
  write_endian<bofstream, uint64_t>(*this, a, switch_endianity);
  return *this;
}

bofstream& bofstream::operator<<(float a)
{
  write_endian<bofstream, float>(*this, a, switch_endianity);
  return *this;
}

bofstream& bofstream::operator<<(double a)
{
  write_endian<bofstream, double>(*this, a, switch_endianity);
  return *this;
}

bofstream& bofstream::operator<<(const char *a)
{
  write(a, strlen(a) + 1);
  return *this;
}

bofstream& bofstream::operator<<(const std::string& a)
{
  write(a.c_str(), a.size() + 1);
  return *this;
}

// ----------------------------------------------------------------------
// bifstream
// ----------------------------------------------------------------------

bifstream::bifstream(const std::string& name, endian e) :
    bfstream_base(e), ifstream(name.c_str(), ios::in | ios::binary) {}

bifstream::bifstream() : bfstream_base(), ifstream() {}

void bifstream::open(const std::string& name, endian e)
{
  if (native_endianity != e)
    switch_endianity = true;
  else
    switch_endianity = false;
  ifstream::open(name.c_str(), ios::in | ios::binary);
}

int bifstream::length() // in bytes
{
  std::streampos pos1, len;
  pos1 = tellg();
  seekg(0, ios::end);
  len = tellg();
  seekg(pos1);
  return int(len);
}

bifstream& bifstream::operator>>(char& a)
{
  get(a);
  return *this;
}

bifstream& bifstream::operator>>(unsigned char& a)
{
  char tmp;
  get(tmp);
  a = tmp;
  return *this;
}

bifstream& bifstream::operator>>(int16_t& a)
{
  read_endian<bifstream, int16_t>(*this, a, switch_endianity);
  return *this;
}

bifstream& bifstream::operator>>(uint16_t& a)
{
  read_endian<bifstream, uint16_t>(*this, a, switch_endianity);
  return *this;
}

bifstream& bifstream::operator>>(int32_t& a)
{
  read_endian<bifstream, int32_t>(*this, a, switch_endianity);
  return *this;
}

bifstream& bifstream::operator>>(uint32_t& a)
{
  read_endian<bifstream, uint32_t>(*this, a, switch_endianity);
  return *this;
}

bifstream& bifstream::operator>>(int64_t& a)
{
  read_endian<bifstream, int64_t>(*this, a, switch_endianity);
  return *this;
}

bifstream& bifstream::operator>>(uint64_t& a)
{
  read_endian<bifstream, uint64_t>(*this, a, switch_endianity);
  return *this;
}

bifstream& bifstream::operator>>(float& a)
{
  read_endian<bifstream, float>(*this, a, switch_endianity);
  return *this;
}

bifstream& bifstream::operator>>(double& a)
{
  read_endian<bifstream, double>(*this, a, switch_endianity);
  return *this;
}

bifstream& bifstream::operator>>(char *a)
{
  getline(a, '\0');
  return *this;
}

bifstream& bifstream::operator>>(std::string& a)
{
  std::getline(*this, a, '\0');
  return *this;
}

// ----------------------------------------------------------------------
// bfstream
// ----------------------------------------------------------------------

bfstream::bfstream(const std::string& name, endian e) :
    bfstream_base(e), fstream(name.c_str(), ios::in | ios::out | ios::binary)
{}

bfstream::bfstream() : bfstream_base(), fstream() {}

void bfstream::open(const std::string& name, bool trnc, endian e)
{
  if (native_endianity != e)
    switch_endianity = true;
  else
    switch_endianity = false;

  if (trnc)
    fstream::open(name.c_str(), ios::in | ios::out | ios::binary
                  | ios::trunc);
  else
    fstream::open(name.c_str(), ios::in | ios::out | ios::binary);
}

void bfstream::open_readonly(const std::string& name, endian e)
{
  if (native_endianity != e)
    switch_endianity = true;
  else
    switch_endianity = false;
  fstream::open(name.c_str(), ios::in | ios::binary);
}

int bfstream::length() // in bytes
{
  std::streampos pos1, len;
  pos1 = tellg();
  seekg(0, ios::end);
  len = tellg();
  seekg(pos1);
  return int(len);
}

bfstream& bfstream::operator<<(char a)
{
  put(a);
  return *this;
}

bfstream& bfstream::operator<<(unsigned char a)
{
  put(static_cast<char>(a));
  return *this;
}

bfstream& bfstream::operator<<(int16_t a)
{
  write_endian<bfstream, int16_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator<<(uint16_t a)
{
  write_endian<bfstream, uint16_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator<<(int32_t a)
{
  write_endian<bfstream, int32_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator<<(uint32_t a)
{
  write_endian<bfstream, uint32_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator<<(int64_t a)
{
  write_endian<bfstream, int64_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator<<(uint64_t a)
{
  write_endian<bfstream, uint64_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator<<(float a)
{
  write_endian<bfstream, float>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator<<(double a)
{
  write_endian<bfstream, double>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator<<(const char *a)
{
  write(a, strlen(a) + 1);
  return *this;
}

bfstream& bfstream::operator<<(const std::string& a)
{
  write(a.c_str(), a.size() + 1);
  return *this;
}


bfstream& bfstream::operator>>(char& a)
{
  get(a);
  return *this;
}

bfstream& bfstream::operator>>(unsigned char& a)
{
  char tmp;
  get(tmp);
  a = tmp;
  return *this;
}

bfstream& bfstream::operator>>(int16_t& a)
{
  read_endian<bfstream, int16_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator>>(uint16_t& a)
{
  read_endian<bfstream, uint16_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator>>(int32_t& a)
{
  read_endian<bfstream, int32_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator>>(uint32_t& a)
{
  read_endian<bfstream, uint32_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator>>(int64_t& a)
{
  read_endian<bfstream, int64_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator>>(uint64_t& a)
{
  read_endian<bfstream, uint64_t>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator>>(float& a)
{
  read_endian<bfstream, float>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator>>(double& a)
{
  read_endian<bfstream, double>(*this, a, switch_endianity);
  return *this;
}

bfstream& bfstream::operator>>(char *a)
{
  getline(a, '\0');
  return *this;
}

bfstream& bfstream::operator>>(std::string& a)
{
  std::getline(*this, a, '\0');
  return *this;
}

} // namespace itpp
