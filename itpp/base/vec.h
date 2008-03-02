/*!
 * \file
 * \brief Templated Vector Class Definitions
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

#ifndef VEC_H
#define VEC_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <itpp/base/itassert.h>
#include <itpp/base/math/misc.h>
#include <itpp/base/copy_vector.h>
#include <itpp/base/factory.h>


namespace itpp {

  // Declaration of Vec
  template<class Num_T> class Vec;
  // Declaration of Mat
  template<class Num_T> class Mat;
  // Declaration of bin
  class bin;

  //-----------------------------------------------------------------------------------
  // Declaration of Vec Friends
  //-----------------------------------------------------------------------------------

  //! Addition of two vectors
  template<class Num_T> const Vec<Num_T> operator+(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
  //! Addition of a vector and a scalar
  template<class Num_T> const Vec<Num_T> operator+(const Vec<Num_T> &v, const Num_T t);
  //! Addition of a scalar and a vector
  template<class Num_T> const Vec<Num_T> operator+(const Num_T t, const Vec<Num_T> &v);

  //! Subtraction of a vector from a vector
  template<class Num_T> const Vec<Num_T> operator-(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
  //! Subtraction of a scalar from a vector
  template<class Num_T> const Vec<Num_T> operator-(const Vec<Num_T> &v, const Num_T t);
  //! Subtraction of vector from scalar. Results in a vector
  template<class Num_T> const Vec<Num_T> operator-(const Num_T t, const Vec<Num_T> &v);
  //! Negation of vector
  template<class Num_T> const Vec<Num_T> operator-(const Vec<Num_T> &v);

  //! Inner (dot) product of two vectors v1 and v2
  template<class Num_T> Num_T dot(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
  //! Inner (dot) product of two vectors v1 and v2
  template<class Num_T> Num_T operator*(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
  /*!
   * \brief Outer product of two vectors v1 and v2
   *
   * When \a v1 and \a v2 are complex vectors (cvec), the third boolean
   * argument \a hermitian can be set to \a true to conjugate \a v2
   * (Matlab's v1 * v2' operation). This parameter is ignored for types
   * other then cvec.
   */
  template<class Num_T>
  const Mat<Num_T> outer_product(const Vec<Num_T> &v1, const Vec<Num_T> &v2,
				 bool hermitian = false);
  //! Multiplication of a vector and a scalar
  template<class Num_T> const Vec<Num_T> operator*(const Vec<Num_T> &v, const Num_T t);
  //! Multiplication of a scalar and a vector. Results in a vector
  template<class Num_T> const Vec<Num_T> operator*(const Num_T t, const Vec<Num_T> &v);

  //! Element-wise multiplication of the two vectors.  Same functionality as Matlab/Octave expression a .* b
  template<class Num_T> const Vec<Num_T> elem_mult(const Vec<Num_T> &a, const Vec<Num_T> &b);
  //! Element-wise multiplication of the three vectors
  template<class Num_T> const Vec<Num_T> elem_mult(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c);
  //! Element-wise multiplication of the four vectors
  template<class Num_T> const Vec<Num_T> elem_mult(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c, const Vec<Num_T> &d);

  //! Element-wise multiplication of the two vectors, storing the result in vector \c out (which is re-sized if necessary)
  template<class Num_T> void elem_mult_out(const Vec<Num_T> &a, const Vec<Num_T> &b, Vec<Num_T> &out);
  //! Element-wise multiplication of the three vectors, storing the result in vector \c out (which is re-sized if necessary)
  template<class Num_T> void elem_mult_out(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c, Vec<Num_T> &out);
  //! Element-wise multiplication of the four vectors, storing the result in vector \c out (which is re-sized if necessary)
  template<class Num_T> void elem_mult_out(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c, const Vec<Num_T> &d, Vec<Num_T> &out);

  //! In-place element-wise multiplication of two vectors. Fast version of b = elem_mult(a,b)
  template<class Num_T> void elem_mult_inplace(const Vec<Num_T> &a, Vec<Num_T> &b);
  //! Element-wise multiplication of two vectors, followed by summation of the resultant elements. Fast version of sum(elem_mult(a,b))
  template<class Num_T> Num_T elem_mult_sum(const Vec<Num_T> &a, const Vec<Num_T> &b);

  //! Division of all elements in \c v with \c t
  template<class Num_T> const Vec<Num_T> operator/(const Vec<Num_T> &v, const Num_T t);
  //! Division of \c t with all elements in \c v
  template<class Num_T> const Vec<Num_T> operator/(const Num_T t, const Vec<Num_T> &v);

  //! Elementwise division of two vectors. Same functionality as Matlab/Octave expression a ./ b
  template<class Num_T> const Vec<Num_T> elem_div(const Vec<Num_T> &a, const Vec<Num_T> &b);
  //! Elementwise division of scalar \c t and vector \c v
  template<class Num_T> const Vec<Num_T> elem_div(const Num_T t, const Vec<Num_T> &v);
  //! Elementwise division of two vectors, storing the result in vector \c out (which is re-sized if necessary)
  template<class Num_T> void elem_div_out(const Vec<Num_T> &a, const Vec<Num_T> &b, Vec<Num_T> &out);
  //! Elementwise division of two vectors, followed by summation of the resultant elements. Fast version of sum(elem_div(a,b))
  template<class Num_T> Num_T elem_div_sum(const Vec<Num_T> &a, const Vec<Num_T> &b);

  //! Append element \c a to the end of the vector \c v
  template<class Num_T> const Vec<Num_T> concat(const Vec<Num_T> &v, const Num_T a);
  //! Concat element \c a to the beginning of the vector \c v
  template<class Num_T> const Vec<Num_T> concat(const Num_T a, const Vec<Num_T> &v);
  //! Concat vectors \c v1 and \c v2
  template<class Num_T> const Vec<Num_T> concat(const Vec<Num_T> &v1,const Vec<Num_T> &v2);
  //! Concat vectors \c v1, \c v2 and \c v3
  template<class Num_T> const Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3);
  //! Concat vectors \c v1, \c v2, \c v3 and \c v4
  template<class Num_T> const Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3, const Vec<Num_T> &v4);
  //! Concat vectors \c v1, \c v2 \c v3, \c v4 and \c v5
  template<class Num_T> const Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3, const Vec<Num_T> &v4, const Vec<Num_T> &v5);

  //-----------------------------------------------------------------------------------
  // Declaration of Vec
  //-----------------------------------------------------------------------------------

  /*!
    \ingroup arr_vec_mat
    \brief Vector Class (Templated)
    \author Tony Ottosson, Tobias Ringstrom, Adam Piatyszek and Conrad Sanderson

    Vectors can be of arbitrarily types, but conversions and functions are
    prepared for \c bin, \c short, \c int, \c double, and \c complex<double>
    vectors and these are predefined as: \c bvec, \c svec, \c ivec, \c vec,
    and \c cvec. \c double and \c complex<double> are \c double and
    \c complex<double> respectively.

    Examples:

    Vector Constructors:
    When constructing a vector without a length (memory) use
    \code vec temp; \endcode
    For construction of a vector of a given length use
    \code vec temp(length); \endcode
    It is also possible to assign the constructed vector the value and size
    of another vector by
    \code vec temp(invector); \endcode
    If you have explicit values you would like to assign to the vector it is
    possible to do this using strings as:
    \code
    vec a("0 0.7 5 9.3"); // that is a = [0, 0.7, 5, 9.3]
    vec a="0 0.7 5 9.3";  // the constructor are called implicitly
    ivec b="0:5";  // that is b = [0, 1, 2, 3, 4, 5]
    vec c="3:2.5:13";  // that is c = [3, 5.5, 8, 10.5, 13]
    \endcode
    It is also possible to change length by
    \code temp.set_size(new_length, false); \endcode
    where \c false is used to indicate that the old values in \c temp
    is not copied. If you like to preserve the values use \c true.

    There are a number of methods to access parts of a vector. Examples are
    \code
    a(5);     // Element number 5
    a(5,9);  // Elements 5, 6, 7, 8, and 9
    a.left(10);  // The 10 most left elements (the first)
    a.right(10); // The 10 most right elements (the last)
    a.mid(5, 7); // 7 elements starting from element 5
    \endcode

    It is also possible to modify parts of a vector as e.g. in
    \code
    a.del(5);    // deletes element number 5
    a.ins(3.4, 9); // inserts the element 3.4 at position 9
    a.replace_mid(12, b); // replaces elements from 12 with the vector b
    \endcode

    It is of course also possible to perform the common linear algebra
    methods such as addition, subtraction, and scalar product (*). Observe
    though, that vectors are assumed to be column-vectors in operations with
    matrices.

    Most elementary functions such as sin(), cosh(), log(), abs(), ..., are
    available as operations on the individual elements of the vectors. Please
    see the individual functions for more details.

    By default, the Vec elements are created using the default constructor for
    the element type. This can be changed by specifying a suitable Factory in
    the Vec constructor call; see Detailed Description for Factory.
  */
  template<class Num_T>
  class Vec {
  public:
    //! The type of the vector values
    typedef Num_T value_type;

    //! Default constructor. An element factory \c f can be specified.
    explicit Vec(const Factory &f = DEFAULT_FACTORY);
    //! Constructor with size parameter. An element factory \c f can be specified.
    explicit Vec(int size, const Factory &f = DEFAULT_FACTORY);
    //! Copy constructor
    Vec(const Vec<Num_T> &v);
    //! Copy constructor, which takes an element factory \c f as an additional argument.
    Vec(const Vec<Num_T> &v, const Factory &f);
    //! Constructor taking a char string as input. An element factory \c f can be specified.
    Vec(const char *values, const Factory &f = DEFAULT_FACTORY);
    //! Constructor taking a string as input. An element factory \c f can be specified.
    Vec(const std::string &values, const Factory &f = DEFAULT_FACTORY);
    //! Constructor taking a C-array as input. Copies all data. An element factory \c f can be specified.
    Vec(const Num_T *c_array, int size, const Factory &f = DEFAULT_FACTORY);

    //! Destructor
    ~Vec();

    //! The size of the vector
    int length() const { return datasize; }
    //! The size of the vector
    int size() const { return datasize; }

    //! Set length of vector. if copy = true then keeping the old values
    void set_size(int size, const bool copy=false);
    //! Set length of vector. if copy = true then keeping the old values
    void set_length(int size, const bool copy=false) { set_size(size, copy); }
    //! Set the vector to the all zero vector
    void zeros();
    //! Set the vector to the all zero vector
    void clear() { zeros(); }
    //! Set the vector to the all one vector
    void ones();
    //! Set the vector equal to the values in the \c str string
    void set(const char *str);
    //! Set the vector equal to the values in the \c str string
    void set(const std::string &str);

    //! C-style index operator. First element is 0
    const Num_T &operator[](int i) const;
    //! Index operator. First element is 0
    const Num_T &operator()(int i) const;
    //! C-style index operator. First element is 0
    Num_T &operator[](int i);
    //! Index operator. First element is 0
    Num_T &operator()(int i);
    //! Sub-vector with elements from \c i1 to \c i2. Index -1 indicates the last element.
    const Vec<Num_T> operator()(int i1, int i2) const;
    //! Sub-vector where the elements are given by the list \c indexlist
    const Vec<Num_T> operator()(const Vec<int> &indexlist) const;

    //! Accessor-style method. First element is 0
    const Num_T &get(int i) const;
    //! Sub-vector with elements from \c i1 to \c i2. Index -1 indicates the last element.
    const Vec<Num_T> get(int i1, int i2) const;
    //! Modifier-style method. First element is 0
    void set(int i, const Num_T &v);

    //! Matrix transpose. Converts to a matrix with the vector in the first row
    Mat<Num_T> transpose() const;
    //! Matrix transpose. Converts to a matrix with the vector in the first row
    Mat<Num_T> T() const { return this->transpose(); }
    //! Hermitian matrix transpose. Converts to a matrix with the conjugate of the vector in the first row
    Mat<Num_T> hermitian_transpose() const;
    //! Hermitian matrix transpose. Converts to a matrix with the conjugate of the vector in the first row
    Mat<Num_T> H() const { return this->hermitian_transpose(); }

    //! Addition of vector
    Vec<Num_T>& operator+=(const Vec<Num_T> &v);
    //! Addition of scalar
    Vec<Num_T>& operator+=(const Num_T t);
    //! Addition of two vectors
    friend const Vec<Num_T> operator+<>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Addition of a vector and a scalar
    friend const Vec<Num_T> operator+<>(const Vec<Num_T> &v, const Num_T t);
    //! Addition of a scalar and a vector
    friend const Vec<Num_T> operator+<>(const Num_T t, const Vec<Num_T> &v);

    //! Subtraction of vector
    Vec<Num_T>& operator-=(const Vec<Num_T> &v);
    //! Subtraction of scalar
    Vec<Num_T>& operator-=(const Num_T t);
    //! Subtraction of \c v2 from \c v1
    friend const Vec<Num_T> operator-<>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Subtraction of scalar from vector
    friend const Vec<Num_T> operator-<>(const Vec<Num_T> &v, const Num_T t);
    //! Sutraction of vector from scalar
    friend const Vec<Num_T> operator-<>(const Num_T t, const Vec<Num_T> &v);
    //! Negation of vector
    friend const Vec<Num_T> operator-<>(const Vec<Num_T> &v);

    //! Multiply with a scalar
    Vec<Num_T>& operator*=(const Num_T t);
    //! Inner (dot) product
    friend Num_T operator*<>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Inner (dot) product
    friend Num_T dot <>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Outer product of two vectors v1 and v2
    friend const Mat<Num_T> outer_product <>(const Vec<Num_T> &v1,
					     const Vec<Num_T> &v2,
					     bool hermitian);
    //! Elementwise multiplication of vector and scalar
    friend const Vec<Num_T> operator*<>(const Vec<Num_T> &v, const Num_T t);
    //! Elementwise multiplication of vector and scalar
    friend const Vec<Num_T> operator*<>(const Num_T t, const Vec<Num_T> &v);

    //! Elementwise multiplication
    friend const Vec<Num_T> elem_mult <>(const Vec<Num_T> &a, const Vec<Num_T> &b);
    //! Elementwise multiplication of three vectors
    friend const Vec<Num_T> elem_mult <>(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c);
    //! Elementwise multiplication of four vectors
    friend const Vec<Num_T> elem_mult <>(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c, const Vec<Num_T> &d);

    //! Elementwise multiplication, storing the result in vector \c out (which is re-sized if necessary)
    friend void elem_mult_out <>(const Vec<Num_T> &a, const Vec<Num_T> &b, Vec<Num_T> &out);
    //! Elementwise multiplication of three vectors, storing the result in vector \c out (which is re-sized if necessary)
    friend void elem_mult_out <>(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c, Vec<Num_T> &out);
    //! Elementwise multiplication of four vectors, storing the result in vector \c out (which is re-sized if necessary)
    friend void elem_mult_out <>(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c, const Vec<Num_T> &d, Vec<Num_T> &out);

    //! In-place element-wise multiplication of two vectors. Fast version of b = elem_mult(a,b)
    friend void elem_mult_inplace <>(const Vec<Num_T> &a, Vec<Num_T> &b);
    //! Element-wise multiplication of two vectors, followed by summation of the resultant elements. Fast version of sum(elem_mult(a,b))
    friend Num_T elem_mult_sum <>(const Vec<Num_T> &a, const Vec<Num_T> &b);

    //! Elementwise division
    Vec<Num_T>& operator/=(const Num_T t);
    //! Elementwise division
    Vec<Num_T>& operator/=(const Vec<Num_T> &v);

    //! Elementwise division
    friend const Vec<Num_T> operator/<>(const Vec<Num_T> &v, const Num_T t);
    //! Elementwise division
    friend const Vec<Num_T> operator/<>(const Num_T t, const Vec<Num_T> &v);

    //! Elementwise division
    friend const Vec<Num_T> elem_div <>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Elementwise division
    friend const Vec<Num_T> elem_div <>(const Num_T t, const Vec<Num_T> &v);
    //! Elementwise division
    friend void elem_div_out <>(const Vec<Num_T> &v1, const Vec<Num_T> &v2, Vec<Num_T> &out);
    //! Elementwise division, followed by summation of the resultant elements. Fast version of sum(elem_mult(a,b))
    friend Num_T elem_div_sum <>(const Vec<Num_T> &a, const Vec<Num_T> &b);

    //! Get the elements in the vector where \c binlist is \c 1
    Vec<Num_T> get(const Vec<bin> &binlist) const;
    //! Get the right \c nr elements from the vector
    Vec<Num_T> right(int nr) const;
    //! Get the left \c nr elements from the vector
    Vec<Num_T> left(int nr) const;
    //! Get the middle part of vector from \c start including \c nr elements
    Vec<Num_T> mid(int start, int nr) const;
    //! Split the vector into two parts at element \c pos. Return the first part and keep the second.
    Vec<Num_T> split(int pos);
    //! Shift in element \c In at position 0 \c n times
    void shift_right(const Num_T In, int n=1);
    //! Shift in vector \c In at position 0
    void shift_right(const Vec<Num_T> &In);
    //! Shift out the \c n left elements and a the same time shift in the element \c at last position \c n times
    void shift_left(const Num_T In, int n=1);
    //! Shift in vector \c In at last position
    void shift_left(const Vec<Num_T> &In);

    //! Append element \c a to the end of the vector \c v
    friend const Vec<Num_T> concat<>(const Vec<Num_T> &v, const Num_T a);
    //! Concat element \c a to the beginning of the vector \c v
    friend const Vec<Num_T> concat<>(const Num_T a, const Vec<Num_T> &v);
    //! Concat vectors \c v1 and \c v2
    friend const Vec<Num_T> concat<>(const Vec<Num_T> &v1,const Vec<Num_T> &v2);
    //! Concat vectors \c v1, \c v2 and \c v3
    friend const Vec<Num_T> concat<>(const Vec<Num_T> &v1, const Vec<Num_T> &v2,
				     const Vec<Num_T> &v3);
    //! Concat vectors \c v1, \c v2, \c v3 and \c v4
    friend const Vec<Num_T> concat<>(const Vec<Num_T> &v1, const Vec<Num_T> &v2,
				     const Vec<Num_T> &v3, const Vec<Num_T> &v4);
    //! Concat vectors \c v1, \c v2, \c v3, \c v4 and \c v5
    friend const Vec<Num_T> concat<>(const Vec<Num_T> &v1, const Vec<Num_T> &v2,
				     const Vec<Num_T> &v3, const Vec<Num_T> &v4,
				     const Vec<Num_T> &v5);

    //! Set subvector defined by indicies \c i1 to \c i2 to vector \c v
    void set_subvector(int i1, int i2, const Vec<Num_T> &v);
    //! Set subvector defined by first index \c i and size of vector \c v to \c v
    void set_subvector(int i, const Vec<Num_T> &v);
    //! Set subvector defined by indicies i1 to i2 to constant t
    void set_subvector(int i1, int i2, const Num_T t);
    //! Replace the elements from \c pos by the vector \c v
    void replace_mid(int pos, const Vec<Num_T> &v);
    //! Delete element number \c index
    void del(int index);
    //! Delete elements from \c i1 to \c i2
    void del(int i1, int i2);
    //! Insert element \c in at \c index
    void ins(int index, const Num_T in);
    //! Insert vector \c in at \c index
    void ins(int index, const Vec<Num_T> &in);

    //! Assign all elements in vector to \c t
    Vec<Num_T>& operator=(const Num_T t);
    //! Assign vector the value and length of \c v
    Vec<Num_T>& operator=(const Vec<Num_T> &v);
    //! Assign vector equal to the 1-dimensional matrix \c m
    Vec<Num_T>& operator=(const Mat<Num_T> &m);
    //! Assign vector the values in the string \c values
    Vec<Num_T>& operator=(const char *values);

    //! Elementwise equal to the scalar
    Vec<bin> operator==(const Num_T value) const;
    //! Elementwise not-equal to the scalar
    Vec<bin> operator!=(const Num_T value) const;
    //! Elementwise less than the scalar
    Vec<bin> operator<(const Num_T value) const;
    //! Elementwise less than and equal to the scalar
    Vec<bin> operator<=(const Num_T value) const;
    //! Elementwise greater than the scalar
    Vec<bin> operator>(const Num_T value) const;
    //! Elementwise greater than and equal to the scalar
    Vec<bin> operator>=(const Num_T value) const;

    //! Compare two vectors. False if wrong sizes or different values
    bool operator==(const Vec<Num_T> &v) const;
    //! Compare two vectors. True if different
    bool operator!=(const Vec<Num_T> &v) const;

    //! Index operator without boundary check. Not recommended to use.
    Num_T &_elem(int i) { return data[i]; }
    //! Index operator without boundary check. Not recommended to use.
    const Num_T &_elem(int i) const { return data[i]; }

    //! Get the pointer to the internal structure. Not recommended to use.
    Num_T *_data() { return data; }
    //! Get the pointer to the internal structure. Not recommended to use.
    const Num_T *_data() const { return data; }

  protected:
    //! Allocate storage for a vector of length \c size.
    void alloc(int size);
    //! Free the storage space allocated by the vector
    void free();

    //! The current number of elements in the vector
    int datasize;
    //! A pointer to the data area
    Num_T *data;
    //! Element factory (set to DEFAULT_FACTORY to use Num_T default constructors only)
    const Factory &factory;
  private:
    // This function is used in set() methods to replace commas with spaces
    std::string replace_commas(const std::string &str);
    //! Check whether index \c i is in the allowed range
    bool in_range(int i) const { return ((i < datasize) && (i >= 0)); }
  };

  //-----------------------------------------------------------------------------------
  // Type definitions of vec, cvec, ivec, svec, and bvec
  //-----------------------------------------------------------------------------------

  /*!
    \relates Vec
    \brief Definition of double vector type
  */
  typedef Vec<double> vec;

  /*!
    \relates Vec
    \brief Definition of complex<double> vector type
  */
  typedef Vec<std::complex<double> > cvec;

  /*!
    \relates Vec
    \brief Definition of integer vector type
  */
  typedef Vec<int> ivec;

  /*!
    \relates Vec
    \brief Definition of short vector type
  */
  typedef Vec<short int> svec;

  /*!
    \relates Vec
    \brief Definition of binary vector type
  */
  typedef Vec<bin> bvec;

} //namespace itpp


#include <itpp/base/mat.h>

namespace itpp {

  //-----------------------------------------------------------------------------------
  // Declaration of input and output streams for Vec
  //-----------------------------------------------------------------------------------

  /*!
    \relates Vec
    \brief Stream output of vector
  */
  template<class Num_T>
  std::ostream &operator<<(std::ostream &os, const Vec<Num_T> &v);

  /*!
    \relates Vec
    \brief Stream input of vector

    The input can be on the form "1 2 3" or "[1 2 3]", i.e. with or without
    brackets. The first form is compatible with the set method, while the
    second form is compatible with the ostream operator. The elements can be
    separated by blank space or commas. "[]" means an empty vector. "1:4"
    means "1 2 3 4". "1:3:10" means every third integer from 1 to 10, i.e.
    "1 4 7 10".
  */
  template<class Num_T>
  std::istream &operator>>(std::istream &is, Vec<Num_T> &v);

  //-----------------------------------------------------------------------------------
  // Implementation of templated Vec members and friends
  //-----------------------------------------------------------------------------------

  template<class Num_T> inline
  void Vec<Num_T>::alloc(int size)
  {
    if (size > 0) {
      create_elements(data, size, factory);
      datasize = size;
    }
    else {
      data = 0;
      datasize = 0;
    }
  }

  template<class Num_T> inline
  void Vec<Num_T>::free()
  {
    destroy_elements(data, datasize);
    datasize = 0;
  }


  template<class Num_T> inline
  Vec<Num_T>::Vec(const Factory &f) : datasize(0), data(0), factory(f) {}

  template<class Num_T> inline
  Vec<Num_T>::Vec(int size, const Factory &f) : datasize(0), data(0), factory(f)
  {
    it_assert_debug(size>=0, "Negative size in Vec::Vec(int)");
    alloc(size);
  }

  template<class Num_T> inline
  Vec<Num_T>::Vec(const Vec<Num_T> &v) : datasize(0), data(0), factory(v.factory)
  {
    alloc(v.datasize);
    copy_vector(datasize, v.data, data);
  }

  template<class Num_T> inline
  Vec<Num_T>::Vec(const Vec<Num_T> &v, const Factory &f) : datasize(0), data(0), factory(f)
  {
    alloc(v.datasize);
    copy_vector(datasize, v.data, data);
  }

  template<class Num_T> inline
  Vec<Num_T>::Vec(const char *values, const Factory &f) : datasize(0), data(0), factory(f)
  {
    set(values);
  }

  template<class Num_T> inline
  Vec<Num_T>::Vec(const std::string &values, const Factory &f) : datasize(0), data(0), factory(f)
  {
    set(values);
  }

  template<class Num_T> inline
  Vec<Num_T>::Vec(const Num_T *c_array, int size, const Factory &f) : datasize(0), data(0), factory(f)
  {
    alloc(size);
    copy_vector(size, c_array, data);
  }

  template<class Num_T> inline
  Vec<Num_T>::~Vec()
  {
    free();
  }

  template<class Num_T>
  void Vec<Num_T>::set_size(int size, bool copy)
  {
    it_assert_debug(size >= 0, "Vec::set_size(): New size must not be negative");
    if (datasize == size)
      return;
    if (copy) {
      // create a temporary pointer to the allocated data
      Num_T* tmp = data;
      // store the current number of elements
      int old_datasize = datasize;
      // check how many elements we need to copy
      int min = datasize < size ? datasize : size;
      // allocate new memory
      alloc(size);
      // copy old elements into a new memory region
      copy_vector(min, tmp, data);
      // initialize the rest of resized vector
      for (int i = min; i < size; ++i)
	data[i] = Num_T(0);
      // delete old elements
      destroy_elements(tmp, old_datasize);
    }
    else {
      free();
      alloc(size);
    }
  }

  template<class Num_T> inline
  const Num_T& Vec<Num_T>::operator[](int i) const
  {
    it_assert_debug(in_range(i), "Vec<>::operator[]: Index out of range");
    return data[i];
  }

  template<class Num_T> inline
  const Num_T& Vec<Num_T>::operator()(int i) const
  {
    it_assert_debug(in_range(i), "Vec<>::operator(): Index out of range");
    return data[i];
  }

  template<class Num_T> inline
  Num_T& Vec<Num_T>::operator[](int i)
  {
    it_assert_debug(in_range(i), "Vec<>::operator[]: Index out of range");
    return data[i];
  }

  template<class Num_T> inline
  Num_T& Vec<Num_T>::operator()(int i)
  {
    it_assert_debug(in_range(i), "Vec<>::operator(): Index out of range");
    return data[i];
  }

  template<class Num_T> inline
  const Vec<Num_T> Vec<Num_T>::operator()(int i1, int i2) const
  {
    if (i1 == -1) i1 = datasize-1;
    if (i2 == -1) i2 = datasize-1;

    it_assert_debug((i1 >= 0) && (i1 <= i2) && (i2 < datasize),
                    "Vec<>::operator()(i1, i2): Indexing out of range");

    Vec<Num_T> s(i2-i1+1);
    copy_vector(s.datasize, data+i1, s.data);

    return s;
  }

  template<class Num_T>
  const Vec<Num_T> Vec<Num_T>::operator()(const Vec<int> &indexlist) const
  {
    Vec<Num_T> temp(indexlist.length());
    for (int i=0;i<indexlist.length();i++) {
      it_assert_debug(in_range(indexlist(i)), "Vec<>::operator()(ivec &): "
                      "Index i=" << i << " out of range");
      temp(i)=data[indexlist(i)];
    }
    return temp;
  }


  template<class Num_T> inline
  const Num_T& Vec<Num_T>::get(int i) const
  {
    it_assert_debug(in_range(i), "Vec<>::get(): Index out of range");
    return data[i];
  }

  template<class Num_T> inline
  const Vec<Num_T> Vec<Num_T>::get(int i1, int i2) const
  {
    return (*this)(i1, i2);
  }


  template<class Num_T> inline
  void Vec<Num_T>::zeros()
  {
    for (int i = 0; i < datasize; i++)
      data[i] = Num_T(0);
  }

  template<class Num_T> inline
  void Vec<Num_T>::ones()
  {
    for (int i = 0; i < datasize; i++)
      data[i] = Num_T(1);
  }

  template<class Num_T> inline
  void Vec<Num_T>::set(int i, const Num_T &v)
  {
    it_assert_debug(in_range(i), "Vec<>::set(i, t): Index out of range");
    data[i] = v;
  }

  //! \cond
  template<>
  void Vec<double>::set(const std::string &str);
  template<>
  void Vec<std::complex<double> >::set(const std::string &str);
  template<>
  void Vec<int>::set(const std::string &str);
  template<>
  void Vec<short int>::set(const std::string &str);
  template<>
  void Vec<bin>::set(const std::string &str);
  //! \endcond

  template<class Num_T>
  void Vec<Num_T>::set(const std::string &str)
  {
    it_error("Vec::set(): Only `double', `complex<double>', `int', "
	     "`short int' and `bin' types supported");
  }

  template<class Num_T>
  void Vec<Num_T>::set(const char *str)
  {
    set(std::string(str));
  }


  template<class Num_T>
  Mat<Num_T> Vec<Num_T>::transpose() const
  {
    Mat<Num_T> temp(1, datasize);
    copy_vector(datasize, data, temp._data());
    return temp;
  }

  //! \cond
  template<>
  Mat<std::complex<double> > Vec<std::complex<double> >::hermitian_transpose() const;
  //! \endcond

  template<class Num_T>
  Mat<Num_T> Vec<Num_T>::hermitian_transpose() const
  {
    Mat<Num_T> temp(1, datasize);
    copy_vector(datasize, data, temp._data());
    return temp;
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator+=(const Vec<Num_T> &v)
  {
    if (datasize == 0) { // if not assigned a size.
      if (this != &v) { // check for self addition
	alloc(v.datasize);
	copy_vector(datasize, v.data, data);
      }
    } else {
      it_assert_debug(datasize == v.datasize, "Vec::operator+=: Wrong sizes");
      for (int i = 0; i < datasize; i++)
	data[i] += v.data[i];
    }
    return *this;
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator+=(const Num_T t)
  {
    for (int i=0;i<datasize;i++)
      data[i]+=t;
    return *this;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator+(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int i;
    Vec<Num_T> r(v1.datasize);

    it_assert_debug(v1.datasize==v2.datasize, "Vec::operator+: wrong sizes");
    for (i=0; i<v1.datasize; i++)
      r.data[i] = v1.data[i] + v2.data[i];

    return r;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator+(const Vec<Num_T> &v, const Num_T t)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = v.data[i] + t;

    return r;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator+(const Num_T t, const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = t + v.data[i];

    return r;
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator-=(const Vec<Num_T> &v)
  {
    if (datasize == 0) { // if not assigned a size.
      if (this != &v) { // check for self decrementation
	alloc(v.datasize);
	for (int i = 0; i < v.datasize; i++)
	  data[i] = -v.data[i];
      }
    } else {
      it_assert_debug(datasize == v.datasize, "Vec::operator-=: Wrong sizes");
      for (int i = 0; i < datasize; i++)
	data[i] -= v.data[i];
    }
    return *this;
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator-=(const Num_T t)
  {
    for (int i=0;i<datasize;i++)
      data[i]-=t;
    return *this;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator-(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int i;
    Vec<Num_T> r(v1.datasize);

    it_assert_debug(v1.datasize==v2.datasize, "Vec::operator-: wrong sizes");
    for (i=0; i<v1.datasize; i++)
      r.data[i] = v1.data[i] - v2.data[i];

    return r;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator-(const Vec<Num_T> &v, const Num_T t)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = v.data[i] - t;

    return r;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator-(const Num_T t, const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = t - v.data[i];

    return r;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator-(const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = -v.data[i];

    return r;
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator*=(const Num_T t)
  {
    scal_vector(datasize, t, data);
    return *this;
  }

#if defined(HAVE_BLAS)
  template<> inline
  double dot(const vec &v1, const vec &v2)
  {
    it_assert_debug(v1.datasize == v2.datasize, "vec::dot: wrong sizes");
    int incr = 1;
    return blas::ddot_(&v1.datasize, v1.data, &incr, v2.data, &incr);
  }

#if defined(HAVE_ZDOTUSUB) || defined(HAVE_ZDOTU_VOID)
  template<> inline
  std::complex<double> dot(const cvec &v1, const cvec &v2)
  {
    it_assert_debug(v1.datasize == v2.datasize, "cvec::dot: wrong sizes");
    int incr = 1;
    std::complex<double> output;
    blas::zdotusub_(&output, &v1.datasize, v1.data, &incr, v2.data, &incr);
    return output;
  }
#endif // HAVE_ZDOTUSUB || HAVE_ZDOTU_VOID
#endif // HAVE_BLAS

  template<class Num_T> inline
  Num_T dot(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int i;
    Num_T r=Num_T(0);

    it_assert_debug(v1.datasize==v2.datasize, "Vec::dot: wrong sizes");
    for (i=0; i<v1.datasize; i++)
      r += v1.data[i] * v2.data[i];

    return r;
  }

  template<class Num_T> inline
  Num_T operator*(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    return dot(v1, v2);
  }


#if defined(HAVE_BLAS)
  template<> inline
  const mat outer_product(const vec &v1, const vec &v2, bool hermitian)
  {
    it_assert_debug((v1.datasize > 0) && (v2.datasize > 0),
		    "Vec::outer_product():: Input vector of zero size");

    mat out(v1.datasize, v2.datasize);
    out.zeros();
    double alpha = 1.0;
    int incr = 1;
    blas::dger_(&v1.datasize, &v2.datasize, &alpha, v1.data, &incr,
		v2.data, &incr, out._data(), &v1.datasize);
    return out;
  }

  template<> inline
  const cmat outer_product(const cvec &v1, const cvec &v2, bool hermitian)
  {
    it_assert_debug((v1.datasize > 0) && (v2.datasize > 0),
		    "Vec::outer_product():: Input vector of zero size");

    cmat out(v1.datasize, v2.datasize);
    out.zeros();
    std::complex<double> alpha(1.0);
    int incr = 1;
    if (hermitian) {
      blas::zgerc_(&v1.datasize, &v2.datasize, &alpha, v1.data, &incr,
		   v2.data, &incr, out._data(), &v1.datasize);
    }
    else {
      blas::zgeru_(&v1.datasize, &v2.datasize, &alpha, v1.data, &incr,
		   v2.data, &incr, out._data(), &v1.datasize);
    }
    return out;
  }
#else
  //! Outer product of two vectors v1 and v2
  template<> inline
  const cmat outer_product(const cvec &v1, const cvec &v2, bool hermitian)
  {
    it_assert_debug((v1.datasize > 0) && (v2.datasize > 0),
		    "Vec::outer_product():: Input vector of zero size");

    cmat out(v1.datasize, v2.datasize);
    if (hermitian) {
      for (int i = 0; i < v1.datasize; ++i) {
	for (int j = 0; j < v2.datasize; ++j) {
	  out(i, j) = v1.data[i] * conj(v2.data[j]);
	}
      }
    }
    else {
      for (int i = 0; i < v1.datasize; ++i) {
	for (int j = 0; j < v2.datasize; ++j) {
	  out(i, j) = v1.data[i] * v2.data[j];
	}
      }
    }
    return out;
  }
#endif // HAVE_BLAS

  template<class Num_T> inline
  const Mat<Num_T> outer_product(const Vec<Num_T> &v1, const Vec<Num_T> &v2,
				 bool hermitian)
  {
    int i, j;

    it_assert_debug((v1.datasize > 0) && (v2.datasize > 0),
		    "Vec::outer_product:: Input vector of zero size");

    Mat<Num_T> r(v1.datasize, v2.datasize);
    for (i=0; i<v1.datasize; i++) {
      for (j=0; j<v2.datasize; j++) {
	r(i,j) = v1.data[i] * v2.data[j];
      }
    }

    return r;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator*(const Vec<Num_T> &v, const Num_T t)
  {
    int i;
    Vec<Num_T> r(v.datasize);
    for (i=0; i<v.datasize; i++)
      r.data[i] = v.data[i] * t;

    return r;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator*(const Num_T t, const Vec<Num_T> &v)
  {
    return operator*(v, t);
  }

  template<class Num_T> inline
  const Vec<Num_T> elem_mult(const Vec<Num_T> &a, const Vec<Num_T> &b)
  {
    Vec<Num_T> out;
    elem_mult_out(a,b,out);
    return out;
  }

  template<class Num_T> inline
  const Vec<Num_T> elem_mult(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c)
  {
    Vec<Num_T> out;
    elem_mult_out(a,b,c,out);
    return out;
  }

  template<class Num_T> inline
  const Vec<Num_T> elem_mult(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c, const Vec<Num_T> &d)
  {
    Vec<Num_T> out;
    elem_mult_out(a,b,c,d,out);
    return out;
  }

  template<class Num_T> inline
  void elem_mult_out(const Vec<Num_T> &a, const Vec<Num_T> &b, Vec<Num_T> &out)
  {
    it_assert_debug(a.datasize==b.datasize, "Vec::elem_mult_out: wrong sizes");

    if(out.datasize != a.datasize)
      out.set_size(a.size());

    for(int i=0; i<a.datasize; i++)
      out.data[i] = a.data[i] * b.data[i];
  }

  template<class Num_T> inline
  void elem_mult_out(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c, Vec<Num_T> &out)
  {
    it_assert_debug(a.datasize==b.datasize==c.datasize, "Vec::elem_mult_out: wrong sizes");

    if(out.datasize != a.datasize)
      out.set_size(a.size());

    for(int i=0; i<a.datasize; i++)
      out.data[i] = a.data[i] * b.data[i] * c.data[i];
  }

  template<class Num_T> inline
  void elem_mult_out(const Vec<Num_T> &a, const Vec<Num_T> &b, const Vec<Num_T> &c, const Vec<Num_T> &d, Vec<Num_T> &out)
  {
    it_assert_debug(a.datasize==b.datasize==c.datasize==d.datasize, "Vec::elem_mult_out: wrong sizes");

    if(out.datasize != a.datasize)
      out.set_size(a.size());

    for(int i=0; i<a.datasize; i++)
      out.data[i] = a.data[i] * b.data[i] * c.data[i] * d.data[i];
  }

  template<class Num_T> inline
  void elem_mult_inplace(const Vec<Num_T> &a, Vec<Num_T> &b)
  {
    it_assert_debug(a.datasize==b.datasize, "Vec::elem_mult_inplace: wrong sizes");

    for(int i=0; i<a.datasize; i++)
      b.data[i] *= a.data[i];
  }

  template<class Num_T> inline
  Num_T elem_mult_sum(const Vec<Num_T> &a, const Vec<Num_T> &b)
  {
    it_assert_debug(a.datasize==b.datasize, "Vec::elem_mult_sum: wrong sizes");

    Num_T acc = 0;

    for(int i=0; i<a.datasize; i++)
      acc += a.data[i] * b.data[i];

    return acc;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator/(const Vec<Num_T> &v, const Num_T t)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = v.data[i] / t;

    return r;
  }

  template<class Num_T> inline
  const Vec<Num_T> operator/(const Num_T t, const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = t / v.data[i];

    return r;
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator/=(const Num_T t)
  {
    for (int i = 0; i < datasize; ++i) {
      data[i] /= t;
    }
    return *this;
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator/=(const Vec<Num_T> &v)
  {
    it_assert_debug(datasize == v.datasize, "Vec::operator/=(): wrong sizes");
    for (int i = 0; i < datasize; ++i) {
      data[i] /= v.data[i];
    }
    return *this;
  }

  template<class Num_T> inline
  const Vec<Num_T> elem_div(const Vec<Num_T> &a, const Vec<Num_T> &b)
  {
    Vec<Num_T> out;
    elem_div_out(a,b,out);
    return out;
  }

  template<class Num_T> inline
  const Vec<Num_T> elem_div(const Num_T t, const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = t / v.data[i];

    return r;
  }

  template<class Num_T> inline
  void elem_div_out(const Vec<Num_T> &a, const Vec<Num_T> &b, Vec<Num_T> &out)
  {
    it_assert_debug(a.datasize==b.datasize, "Vecelem_div_out: wrong sizes");

    if(out.datasize != a.datasize)
      out.set_size(a.size());

    for(int i=0; i<a.datasize; i++)
      out.data[i] = a.data[i] / b.data[i];
  }

  template<class Num_T> inline
  Num_T elem_div_sum(const Vec<Num_T> &a, const Vec<Num_T> &b)
  {
    it_assert_debug(a.datasize==b.datasize, "Vec::elem_div_sum: wrong sizes");

    Num_T acc = 0;

    for(int i=0; i<a.datasize; i++)
      acc += a.data[i] / b.data[i];

    return acc;
  }

  template<class Num_T>
  Vec<Num_T> Vec<Num_T>::get(const Vec<bin> &binlist) const
  {
    int size = binlist.size();
    it_assert_debug(datasize == size, "Vec::get(bvec &): wrong sizes");
    Vec<Num_T> temp(size);
    int j = 0;
    for (int i = 0; i < size; ++i) {
      if (binlist(i) == bin(1)) {
	temp(j) = data[i];
	j++;
      }
    }
    temp.set_size(j, true);
    return temp;
  }

  template<class Num_T> inline
  Vec<Num_T> Vec<Num_T>::right(int nr) const
  {
    it_assert_debug(nr <= datasize, "Vec::right(): index out of range");
    Vec<Num_T> temp(nr);
    if (nr > 0) {
      copy_vector(nr, &data[datasize-nr], temp.data);
    }
    return temp;
  }

  template<class Num_T> inline
  Vec<Num_T> Vec<Num_T>::left(int nr) const
  {
    it_assert_debug(nr <= datasize, "Vec::left(): index out of range");
    Vec<Num_T> temp(nr);
    if (nr > 0) {
      copy_vector(nr, data, temp.data);
    }
    return temp;
  }

  template<class Num_T> inline
  Vec<Num_T> Vec<Num_T>::mid(int start, int nr) const
  {
    it_assert_debug((start >= 0) && ((start+nr) <= datasize),
		    "Vec::mid(): indexing out of range");
    Vec<Num_T> temp(nr);
    if (nr > 0) {
      copy_vector(nr, &data[start], temp.data);
    }
    return temp;
  }

  template<class Num_T>
  Vec<Num_T> Vec<Num_T>::split(int pos)
  {
    it_assert_debug(in_range(pos), "Vec<>::split(): Index out of range");
    Vec<Num_T> temp1(pos);
    Vec<Num_T> temp2(datasize-pos);
    copy_vector(pos, data, temp1.data);
    copy_vector(datasize-pos, &data[pos], temp2.data);
    (*this) = temp2;
    return temp1;
  }

  template<class Num_T>
  void Vec<Num_T>::shift_right(const Num_T In, int n)
  {
    int i=datasize;

    it_assert_debug(n>=0, "Vec::shift_right: index out of range");
    while (--i >= n)
      data[i] = data[i-n];
    while (i >= 0)
      data[i--] = In;
  }

  template<class Num_T>
  void Vec<Num_T>::shift_right(const Vec<Num_T> &In)
  {
    int	i;

    for (i=datasize-1; i>=In.datasize; i--)
      data[i]=data[i-In.datasize];
    for (i=0; i<In.datasize; i++)
      data[i]=In[i];
  }

  template<class Num_T>
  void Vec<Num_T>::shift_left(const Num_T In, int n)
  {
    int i;

    it_assert_debug(n>=0, "Vec::shift_left: index out of range");
    for (i=0; i<datasize-n; i++)
      data[i] = data[i+n];
    while (i < datasize)
      data[i++] = In;
  }

  template<class Num_T>
  void Vec<Num_T>::shift_left(const Vec<Num_T> &In)
  {
    int	i;

    for (i=0; i<datasize-In.datasize; i++)
      data[i]=data[i+In.datasize];
    for (i=datasize-In.datasize; i<datasize; i++)
      data[i]=In[i-datasize+In.datasize];
  }

  template<class Num_T>
  const Vec<Num_T> concat(const Vec<Num_T> &v, const Num_T a)
  {
    int size = v.size();
    Vec<Num_T> temp(size+1);
    copy_vector(size, v.data, temp.data);
    temp(size) = a;
    return temp;
  }

  template<class Num_T>
  const Vec<Num_T> concat(const Num_T a, const Vec<Num_T> &v)
  {
    int size = v.size();
    Vec<Num_T> temp(size+1);
    temp(0) = a;
    copy_vector(size, v.data, &temp.data[1]);
    return temp;
  }

  template<class Num_T>
  const Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int size1 = v1.size();
    int size2 = v2.size();
    Vec<Num_T> temp(size1+size2);
    copy_vector(size1, v1.data, temp.data);
    copy_vector(size2, v2.data, &temp.data[size1]);
    return temp;
  }

  template<class Num_T>
  const Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2,
			  const Vec<Num_T> &v3)
  {
    int size1 = v1.size();
    int size2 = v2.size();
    int size3 = v3.size();
    Vec<Num_T> temp(size1+size2+size3);
    copy_vector(size1, v1.data, temp.data);
    copy_vector(size2, v2.data, &temp.data[size1]);
    copy_vector(size3, v3.data, &temp.data[size1+size2]);
    return temp;
  }

  template<class Num_T>
  const Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2,
			  const Vec<Num_T> &v3, const Vec<Num_T> &v4)
  {
    int size1 = v1.size();
    int size2 = v2.size();
    int size3 = v3.size();
    int size4 = v4.size();
    Vec<Num_T> temp(size1+size2+size3+size4);
    copy_vector(size1, v1.data, temp.data);
    copy_vector(size2, v2.data, &temp.data[size1]);
    copy_vector(size3, v3.data, &temp.data[size1+size2]);
    copy_vector(size4, v4.data, &temp.data[size1+size2+size3]);
    return temp;
  }

  template<class Num_T>
  const Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2,
			  const Vec<Num_T> &v3, const Vec<Num_T> &v4,
			  const Vec<Num_T> &v5)
  {
    int size1 = v1.size();
    int size2 = v2.size();
    int size3 = v3.size();
    int size4 = v4.size();
    int size5 = v5.size();
    Vec<Num_T> temp(size1+size2+size3+size4+size5);
    copy_vector(size1, v1.data, temp.data);
    copy_vector(size2, v2.data, &temp.data[size1]);
    copy_vector(size3, v3.data, &temp.data[size1+size2]);
    copy_vector(size4, v4.data, &temp.data[size1+size2+size3]);
    copy_vector(size5, v5.data, &temp.data[size1+size2+size3+size4]);
    return temp;
  }

  template<class Num_T> inline
  void Vec<Num_T>::set_subvector(int i1, int i2, const Vec<Num_T> &v)
  {
    if (i1 == -1) i1 = datasize-1;
    if (i2 == -1) i2 = datasize-1;

    it_assert_debug(i1>=0 && i2>=0 && i1<datasize && i2<datasize, "Vec::set_subvector(): indicies out of range");
    it_assert_debug(i2>=i1, "Vec::set_subvector(): i2 >= i1 necessary");
    it_assert_debug(i2-i1+1 == v.datasize, "Vec::set_subvector(): wrong sizes");

    copy_vector(v.datasize, v.data, data+i1);
  }

  template<class Num_T> inline
  void Vec<Num_T>:: set_subvector(int i, const Vec<Num_T> &v)
  {
    it_assert_debug((i >= 0) && (i + v.datasize <= datasize),
                    "Vec<>::set_subvector(int, const Vec<> &): "
                    "Index out of range or too long input vector");
    copy_vector(v.datasize, v.data, data+i);
  }

  template<class Num_T>
  void Vec<Num_T>::set_subvector(int i1, int i2, const Num_T t)
  {
    if (i1 == -1) i1 = datasize-1;
    if (i2 == -1) i2 = datasize-1;
    it_assert_debug((i1 >= 0) && (i1 <= i2) && (i2 < datasize),
                    "Vec<>::set_subvector(int, int, Num_T): Indexing out "
                    "of range");
    for (int i=i1;i<=i2;i++)
      data[i] = t;
  }

  template<class Num_T>
  void Vec<Num_T>::replace_mid(int pos, const Vec<Num_T> &v)
  {
    it_assert_debug((pos>=0) && ((pos+v.length())<=datasize), "Vec::replace_mid: indexing out of range");
    copy_vector(v.datasize, v.data, &data[pos]);
  }

  template<class Num_T>
  void Vec<Num_T>::del(int index)
  {
    it_assert_debug(in_range(index), "Vec<>::del(int): Index out of range");
    Vec<Num_T> temp(*this);
    set_size(datasize-1, false);
    copy_vector(index, temp.data, data);
    copy_vector(datasize-index, &temp.data[index+1], &data[index]);
  }

  template<class Num_T>
  void  Vec<Num_T>::del(int i1, int i2)
  {
    if (i1 == -1) i1 = datasize-1;
    if (i2 == -1) i2 = datasize-1;
    it_assert_debug((i1 >= 0) && (i1 <= i2) && (i2 < datasize),
                    "Vec<>::del(int, int): Indexing out of range");
    Vec<Num_T> temp(*this);
    int new_size = datasize-(i2-i1+1);
    set_size(new_size, false);
    copy_vector(i1, temp.data, data);
    copy_vector(datasize-i1, &temp.data[i2+1], &data[i1]);
  }

  template<class Num_T>
  void Vec<Num_T>::ins(int index, const Num_T in)
  {
    it_assert_debug((index >= 0) && (index <= datasize),
                    "Vec<>::ins(): Index out of range");
    Vec<Num_T> Temp(*this);

    set_size(datasize+1, false);
    copy_vector(index, Temp.data, data);
    data[index]=in;
    copy_vector(Temp.datasize-index, Temp.data+index, data+index+1);
  }

  template<class Num_T>
  void Vec<Num_T>::ins(int index, const Vec<Num_T> &in)
  {
    it_assert_debug((index >= 0) && (index <= datasize),
                    "Vec<>::ins(): Index out of range");
    Vec<Num_T> Temp(*this);

    set_size(datasize+in.length(), false);
    copy_vector(index, Temp.data, data);
    copy_vector(in.size(), in.data, &data[index]);
    copy_vector(Temp.datasize-index, Temp.data+index, data+index+in.size());
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator=(const Num_T t)
  {
    for (int i=0;i<datasize;i++)
      data[i] = t;
    return *this;
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator=(const Vec<Num_T> &v)
  {
    if (this != &v) {
      set_size(v.datasize, false);
      copy_vector(datasize, v.data, data);
    }
    return *this;
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator=(const Mat<Num_T> &m)
  {
    if (m.cols() == 1) {
      set_size(m.rows(), false);
      copy_vector(m.rows(), m._data(), data);
    }
    else if (m.rows() == 1) {
      set_size(m.cols(), false);
      copy_vector(m.cols(), m._data(), m.rows(), data, 1);
    }
    else
      it_error("Vec<>::operator=(Mat<Num_T> &): Wrong size of input matrix");
    return *this;
  }

  template<class Num_T> inline
  Vec<Num_T>& Vec<Num_T>::operator=(const char *values)
  {
    set(values);
    return *this;
  }

  //! \cond
  template<>
  bvec Vec<std::complex<double> >::operator==(std::complex<double>) const;
  //! \endcond

  template<class Num_T>
  bvec Vec<Num_T>::operator==(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec::operator==: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)==value);

    return temp;
  }

  //! \cond
  template<>
  bvec Vec<std::complex<double> >::operator!=(std::complex<double>) const;
  //! \endcond

  template<class Num_T>
  bvec Vec<Num_T>::operator!=(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec::operator!=: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)!=value);

    return temp;
  }

  //! \cond
  template<>
  bvec Vec<std::complex<double> >::operator<(std::complex<double>) const;
  //! \endcond

  template<class Num_T>
  bvec Vec<Num_T>::operator<(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec::operator<: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)<value);

    return temp;
  }

  //! \cond
  template<>
  bvec Vec<std::complex<double> >::operator<=(std::complex<double>) const;
  //! \endcond

  template<class Num_T>
  bvec Vec<Num_T>::operator<=(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec::operator<=: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)<=value);

    return temp;
  }

  //! \cond
  template<>
  bvec Vec<std::complex<double> >::operator>(std::complex<double>) const;
  //! \endcond

  template<class Num_T>
  bvec Vec<Num_T>::operator>(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec::operator>: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)>value);

    return temp;
  }

  //! \cond
  template<>
  bvec Vec<std::complex<double> >::operator>=(std::complex<double>) const;
  //! \endcond

  template<class Num_T>
  bvec Vec<Num_T>::operator>=(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec::operator>=: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)>=value);

    return temp;
  }

  template<class Num_T>
  bool Vec<Num_T>::operator==(const Vec<Num_T> &invector) const
  {
    // OBS ! if wrong size, return false
    if (datasize!=invector.datasize) return false;
    for (int i=0;i<datasize;i++) {
      if (data[i]!=invector.data[i]) return false;
    }
    return true;
  }

  template<class Num_T>
  bool Vec<Num_T>::operator!=(const Vec<Num_T> &invector) const
  {
    if (datasize!=invector.datasize) return true;
    for (int i=0;i<datasize;i++) {
      if (data[i]!=invector.data[i]) return true;
    }
    return false;
  }

  //! Output stream operator of a vector \c v
  template<class Num_T>
  std::ostream &operator<<(std::ostream &os, const Vec<Num_T> &v)
  {
    int i, sz=v.length();

    os << "[" ;
    for (i=0; i<sz; i++) {
      os << v(i) ;
      if (i < sz-1)
	os << " ";
    }
    os << "]" ;

    return os;
  }

  //! Input stream operator to read a vector
  template<class Num_T>
  std::istream &operator>>(std::istream &is, Vec<Num_T> &v)
  {
    std::ostringstream buffer;
    bool started = false;
    bool finished = false;
    bool brackets = false;
    char c;

    while (!finished) {
      if (is.eof()) {
        finished = true;
      } else {
        c = is.get();

        if (is.eof() || (c == '\n')) {
          if (brackets) {
            // Right bracket missing
            is.setstate(std::ios_base::failbit);
            finished = true;
          } else if (!((c == '\n') && !started)) {
            finished = true;
          }
        } else if ((c == ' ') || (c == '\t')) {
          if (started) {
            buffer << ' ';
          }
        } else if (c == '[') {
          if (started) {
            // Unexpected left bracket
            is.setstate(std::ios_base::failbit);
            finished = true;
          } else {
            started = true;
            brackets = true;
          }
        } else if (c == ']') {
          if (!started || !brackets) {
            // Unexpected right bracket
            is.setstate(std::ios_base::failbit);
            finished = true;
          } else {
            finished = true;
          }
          while (!is.eof() && (((c = is.peek()) == ' ') || (c == '\t'))) {
            is.get();
          }
          if (!is.eof() && (c == '\n')) {
            is.get();
          }
        } else {
          started = true;
          buffer << c;
        }
      }
    }

    if (!started) {
      v.set_size(0, false);
    } else {
      v.set(buffer.str());
    }

    return is;
  }

  //! \cond

  // ----------------------------------------------------------------------
  // Instantiations
  // ----------------------------------------------------------------------

#ifdef HAVE_EXTERN_TEMPLATE

  extern template class Vec<double>;
  extern template class Vec<int>;
  extern template class Vec<short int>;
  extern template class Vec<std::complex<double> >;
  extern template class Vec<bin>;

  // addition operator

  extern template const vec operator+(const vec &v1, const vec &v2);
  extern template const cvec operator+(const cvec &v1, const cvec &v2);
  extern template const ivec operator+(const ivec &v1, const ivec &v2);
  extern template const svec operator+(const svec &v1, const svec &v2);
  extern template const bvec operator+(const bvec &v1, const bvec &v2);

  extern template const vec operator+(const vec &v1, double t);
  extern template const cvec operator+(const cvec &v1, std::complex<double> t);
  extern template const ivec operator+(const ivec &v1, int t);
  extern template const svec operator+(const svec &v1, short t);
  extern template const bvec operator+(const bvec &v1, bin t);

  extern template const vec operator+(double t, const vec &v1);
  extern template const cvec operator+(std::complex<double> t, const cvec &v1);
  extern template const ivec operator+(int t, const ivec &v1);
  extern template const svec operator+(short t, const svec &v1);
  extern template const bvec operator+(bin t, const bvec &v1);

  // subraction operator

  extern template const vec operator-(const vec &v1, const vec &v2);
  extern template const cvec operator-(const cvec &v1, const cvec &v2);
  extern template const ivec operator-(const ivec &v1, const ivec &v2);
  extern template const svec operator-(const svec &v1, const svec &v2);
  extern template const bvec operator-(const bvec &v1, const bvec &v2);

  extern template const vec operator-(const vec &v, double t);
  extern template const cvec operator-(const cvec &v, std::complex<double> t);
  extern template const ivec operator-(const ivec &v, int t);
  extern template const svec operator-(const svec &v, short t);
  extern template const bvec operator-(const bvec &v, bin t);

  extern template const vec operator-(double t, const vec &v);
  extern template const cvec operator-(std::complex<double> t, const cvec &v);
  extern template const ivec operator-(int t, const ivec &v);
  extern template const svec operator-(short t, const svec &v);
  extern template const bvec operator-(bin t, const bvec &v);

  // unary minus

  extern template const vec operator-(const vec &v);
  extern template const cvec operator-(const cvec &v);
  extern template const ivec operator-(const ivec &v);
  extern template const svec operator-(const svec &v);
  extern template const bvec operator-(const bvec &v);

  // multiplication operator

#if !defined(HAVE_BLAS)
  extern template double dot(const vec &v1, const vec &v2);
#if !(defined(HAVE_ZDOTUSUB) || defined(HAVE_ZDOTU_VOID))
  extern template std::complex<double> dot(const cvec &v1, const cvec &v2);
#endif // !(HAVE_ZDOTUSUB || HAVE_ZDOTU_VOID)
#endif // HAVE_BLAS
  extern template int dot(const ivec &v1, const ivec &v2);
  extern template short dot(const svec &v1, const svec &v2);
  extern template bin dot(const bvec &v1, const bvec &v2);

#if !defined(HAVE_BLAS)
  extern template double operator*(const vec &v1, const vec &v2);
  extern template std::complex<double> operator*(const cvec &v1,
                                                 const cvec &v2);
#endif
  extern template int operator*(const ivec &v1, const ivec &v2);
  extern template short operator*(const svec &v1, const svec &v2);
  extern template bin operator*(const bvec &v1, const bvec &v2);

#if !defined(HAVE_BLAS)
  extern template const mat outer_product(const vec &v1, const vec &v2,
					  bool hermitian);
#endif
  extern template const imat outer_product(const ivec &v1, const ivec &v2,
					   bool hermitian);
  extern template const smat outer_product(const svec &v1, const svec &v2,
					   bool hermitian);
  extern template const bmat outer_product(const bvec &v1, const bvec &v2,
					   bool hermitian);

  extern template const vec operator*(const vec &v, double t);
  extern template const cvec operator*(const cvec &v, std::complex<double> t);
  extern template const ivec operator*(const ivec &v, int t);
  extern template const svec operator*(const svec &v, short t);
  extern template const bvec operator*(const bvec &v, bin t);

  extern template const vec operator*(double t, const vec &v);
  extern template const cvec operator*(std::complex<double> t, const cvec &v);
  extern template const ivec operator*(int t, const ivec &v);
  extern template const svec operator*(short t, const svec &v);
  extern template const bvec operator*(bin t, const bvec &v);

  // elementwise multiplication

  extern template const vec elem_mult(const vec &a, const vec &b);
  extern template const cvec elem_mult(const cvec &a, const cvec &b);
  extern template const ivec elem_mult(const ivec &a, const ivec &b);
  extern template const svec elem_mult(const svec &a, const svec &b);
  extern template const bvec elem_mult(const bvec &a, const bvec &b);

  extern template void elem_mult_out(const vec &a, const vec &b, vec &out);
  extern template void elem_mult_out(const cvec &a, const cvec &b, cvec &out);
  extern template void elem_mult_out(const ivec &a, const ivec &b, ivec &out);
  extern template void elem_mult_out(const svec &a, const svec &b, svec &out);
  extern template void elem_mult_out(const bvec &a, const bvec &b, bvec &out);

  extern template const vec elem_mult(const vec &a, const vec &b,
                                      const vec &c);
  extern template const cvec elem_mult(const cvec &a, const cvec &b,
                                       const cvec &c);
  extern template const ivec elem_mult(const ivec &a, const ivec &b,
                                       const ivec &c);
  extern template const svec elem_mult(const svec &a, const svec &b,
                                       const svec &c);
  extern template const bvec elem_mult(const bvec &a, const bvec &b,
                                       const bvec &c);

  extern template void elem_mult_out(const vec &a, const vec &b, const vec &c,
                                     vec &out);
  extern template void elem_mult_out(const cvec &a, const cvec &b,
                                     const cvec &c, cvec &out);
  extern template void elem_mult_out(const ivec &a, const ivec &b,
                                     const ivec &c, ivec &out);
  extern template void elem_mult_out(const svec &a, const svec &b,
                                     const svec &c, svec &out);
  extern template void elem_mult_out(const bvec &a, const bvec &b,
                                     const bvec &c, bvec &out);

  extern template const vec elem_mult(const vec &a, const vec &b, const vec &c,
                                      const vec &d);
  extern template const cvec elem_mult(const cvec &a, const cvec &b,
                                       const cvec &c, const cvec &d);
  extern template const ivec elem_mult(const ivec &a, const ivec &b,
                                       const ivec &c, const ivec &d);
  extern template const svec elem_mult(const svec &a, const svec &b,
                                       const svec &c, const svec &d);
  extern template const bvec elem_mult(const bvec &a, const bvec &b,
                                       const bvec &c, const bvec &d);

  extern template void elem_mult_out(const vec &a, const vec &b, const vec &c,
                                     const vec &d, vec &out);
  extern template void elem_mult_out(const cvec &a, const cvec &b,
                                     const cvec &c, const cvec &d, cvec &out);
  extern template void elem_mult_out(const ivec &a, const ivec &b,
                                     const ivec &c, const ivec &d, ivec &out);
  extern template void elem_mult_out(const svec &a, const svec &b,
                                     const svec &c, const svec &d, svec &out);
  extern template void elem_mult_out(const bvec &a, const bvec &b,
                                     const bvec &c, const bvec &d, bvec &out);

  // in-place elementwise multiplication

  extern template void elem_mult_inplace(const vec &a, vec &b);
  extern template void elem_mult_inplace(const cvec &a, cvec &b);
  extern template void elem_mult_inplace(const ivec &a, ivec &b);
  extern template void elem_mult_inplace(const svec &a, svec &b);
  extern template void elem_mult_inplace(const bvec &a, bvec &b);

  // elementwise multiplication followed by summation

  extern template double elem_mult_sum(const vec &a, const vec &b);
  extern template std::complex<double> elem_mult_sum(const cvec &a,
                                                     const cvec &b);
  extern template int elem_mult_sum(const ivec &a, const ivec &b);
  extern template short elem_mult_sum(const svec &a, const svec &b);
  extern template bin elem_mult_sum(const bvec &a, const bvec &b);

  // division operator

  extern template const vec operator/(const vec &v, double t);
  extern template const cvec operator/(const cvec &v, std::complex<double> t);
  extern template const ivec operator/(const ivec &v, int t);
  extern template const svec operator/(const svec &v, short t);
  extern template const bvec operator/(const bvec &v, bin t);

  extern template const vec operator/(double t, const vec &v);
  extern template const cvec operator/(std::complex<double> t, const cvec &v);
  extern template const ivec operator/(int t, const ivec &v);
  extern template const svec operator/(short t, const svec &v);
  extern template const bvec operator/(bin t, const bvec &v);

  // elementwise division operator

  extern template const vec elem_div(const vec &a, const vec &b);
  extern template const cvec elem_div(const cvec &a, const cvec &b);
  extern template const ivec elem_div(const ivec &a, const ivec &b);
  extern template const svec elem_div(const svec &a, const svec &b);
  extern template const bvec elem_div(const bvec &a, const bvec &b);

  extern template const vec elem_div(double t, const vec &v);
  extern template const cvec elem_div(std::complex<double> t, const cvec &v);
  extern template const ivec elem_div(int t, const ivec &v);
  extern template const svec elem_div(short t, const svec &v);
  extern template const bvec elem_div(bin t, const bvec &v);

  extern template void elem_div_out(const vec &a, const vec &b, vec &out);
  extern template void elem_div_out(const cvec &a, const cvec &b, cvec &out);
  extern template void elem_div_out(const ivec &a, const ivec &b, ivec &out);
  extern template void elem_div_out(const svec &a, const svec &b, svec &out);
  extern template void elem_div_out(const bvec &a, const bvec &b, bvec &out);

  // elementwise division followed by summation

  extern template double elem_div_sum(const vec &a, const vec &b);
  extern template std::complex<double> elem_div_sum(const cvec &a,
                                                    const cvec &b);
  extern template int elem_div_sum(const ivec &a, const ivec &b);
  extern template short elem_div_sum(const svec &a, const svec &b);
  extern template bin elem_div_sum(const bvec &a, const bvec &b);

  // concat operator

  extern template const vec concat(const vec &v, double a);
  extern template const cvec concat(const cvec &v, std::complex<double> a);
  extern template const ivec concat(const ivec &v, int a);
  extern template const svec concat(const svec &v, short a);
  extern template const bvec concat(const bvec &v, bin a);

  extern template const vec concat(double a, const vec &v);
  extern template const cvec concat(std::complex<double> a, const cvec &v);
  extern template const ivec concat(int a, const ivec &v);
  extern template const svec concat(short a, const svec &v);
  extern template const bvec concat(bin a, const bvec &v);

  extern template const vec concat(const vec &v1, const vec &v2);
  extern template const cvec concat(const cvec &v1, const cvec &v2);
  extern template const ivec concat(const ivec &v1, const ivec &v2);
  extern template const svec concat(const svec &v1, const svec &v2);
  extern template const bvec concat(const bvec &v1, const bvec &v2);

  extern template const vec concat(const vec &v1, const vec &v2,
                                   const vec &v3);
  extern template const cvec concat(const cvec &v1, const cvec &v2,
                                    const cvec &v3);
  extern template const ivec concat(const ivec &v1, const ivec &v2,
                                    const ivec &v3);
  extern template const svec concat(const svec &v1, const svec &v2,
                                    const svec &v3);
  extern template const bvec concat(const bvec &v1, const bvec &v2,
                                    const bvec &v3);

  extern template const vec concat(const vec &v1, const vec &v2, const vec &v3,
                                   const vec &v4);
  extern template const cvec concat(const cvec &v1, const cvec &v2,
                                    const cvec &v3, const cvec &v4);
  extern template const ivec concat(const ivec &v1, const ivec &v2,
                                    const ivec &v3, const ivec &v4);
  extern template const svec concat(const svec &v1, const svec &v2,
                                    const svec &v3, const svec &v4);
  extern template const bvec concat(const bvec &v1, const bvec &v2,
                                    const bvec &v3, const bvec &v4);

  extern template const vec concat(const vec &v1, const vec &v2, const vec &v3,
                                   const vec &v4, const vec &v5);
  extern template const cvec concat(const cvec &v1, const cvec &v2,
                                    const cvec &v3, const cvec &v4,
                                    const cvec &v5);
  extern template const ivec concat(const ivec &v1, const ivec &v2,
                                    const ivec &v3, const ivec &v4,
                                    const ivec &v5);
  extern template const svec concat(const svec &v1, const svec &v2,
                                    const svec &v3, const svec &v4,
                                    const svec &v5);
  extern template const bvec concat(const bvec &v1, const bvec &v2,
                                    const bvec &v3, const bvec &v4,
                                    const bvec &v5);

  // I/O streams

  extern template std::ostream &operator<<(std::ostream& os, const vec &vect);
  extern template std::ostream &operator<<(std::ostream& os, const cvec &vect);
  extern template std::ostream &operator<<(std::ostream& os, const svec &vect);
  extern template std::ostream &operator<<(std::ostream& os, const ivec &vect);
  extern template std::ostream &operator<<(std::ostream& os, const bvec &vect);
  extern template std::istream &operator>>(std::istream& is, vec &vect);
  extern template std::istream &operator>>(std::istream& is, cvec &vect);
  extern template std::istream &operator>>(std::istream& is, svec &vect);
  extern template std::istream &operator>>(std::istream& is, ivec &vect);
  extern template std::istream &operator>>(std::istream& is, bvec &vect);

#endif // HAVE_EXTERN_TEMPLATE

  //! \endcond

} // namespace itpp

#endif // #ifndef VEC_H
