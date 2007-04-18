/*!
 * \file 
 * \brief Definitions of Bit Error Rate Counter (BERC) and 
 *        BLock Error Rate Counter (BLERC) classes
 * \author Pal Frenger and Adam Piatyszek
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
 
#ifndef ERROR_COUNTERS_H
#define ERROR_COUNTERS_H

#include <itpp/base/vec.h>


namespace itpp {

  /*! 
    \brief Bit Error Rate Counter (BERC) Class

    Example:
    \code
    #include <itpp/itcomm.h>

    int main() {
      //Initiate the Bit Error Counter
      BERC berc;

      //Initiate a Binary Symetric Channel with cross-over probability 0.1
      BSC binary_symetric_channel(0.1);
  
      bvec transmitted_bits = randb(100);
      bvec received_bits = binary_symetric_channel(transmitted_bits);

      //Count the number of bit errors
      berc.count(transmitted_bits, received_bits);

      cout << "Estimated bit error probability is " << berc.get_errorrate() 
           << endl;

      return 0;
    } 
    \endcode
  */
  class BERC {
  public:
    /*! 
      \brief Constructor for the berc class.

      <ul>
      <li> \a delay is positive if \a in2 is a delayed replica of 
      \a in1 and negative otherwise. </li>
      <li> \a ignorefirst and \a ignorelast may be used if errors in
      the begining and/or the end is to be ignored.</li> 
      </ul>
    */
    BERC(int indelay = 0, int inignorefirst = 0, int inignorelast = 0);
    //! Cumulative error counter
    void count(const bvec &in1, const bvec &in2);
    //! Run this member function if the delay between \a in1 and
    //! \a in2 is unknown. 
    void estimate_delay(const bvec &in1, const bvec &in2, int mindelay = -100,
			int maxdelay = 100);
    //! Clears the bit error counter
    void clear() { errors = 0; corrects = 0; }
    //! Writes an error report
    void report();
    //! Return the \a delay, assumed or estimated, between \a in1 and \a in2.
    int get_delay() { return delay; }
    //! Returns the counted number of bit errors
    double get_errors() { return errors; }
    //! Returns the counted number of corectly received bits
    double get_corrects() { return corrects; }
    //! Returns the total number of bits processed
    double get_total_bits() { return (errors+corrects); }
    //! Returns the estimated bit error rate.
    double get_errorrate() { return (errors / (corrects + errors)); }
    /*!
      \brief static function to allow simple and fast count of bit-errors

      Returns the number of errors between in1 and in2. Typical usage:
      \code
      bvec in1 = randb(100);
      bvec in2 = randb(100);
      double errors = BERC::count_errors(in1, in2);
      \endcode
    */
    static double count_errors(const bvec &in1, const bvec &in2, 
			       int indelay = 0, int inignorefirst = 0, 
			       int inignorelast = 0);

  private:
    int delay;
    int ignorefirst;
    int ignorelast;
    double errors;
    double corrects;
  };

  /*!
    \brief Class for counting block error rates.
  
    Use this class to count block errors in binary vectors.
  */
  class BLERC {
  public:
    //! Class constructor
    BLERC(void);
    //! Specialised constructor
    BLERC(int blocksize);
    //! Set the block size
    void set_blocksize(int inblocksize, bool clear = true);
    //! Calculate the number of block errors between \a in1 and \a in2
    void count(const bvec &in1, const bvec &in2);
    //! Clear the block error counter
    void clear() { errors = 0; corrects = 0; }
    //! Returns the number of block errors
    double get_errors() { return errors; }
    //! Returns the number of correct blocks
    double get_corrects() { return corrects; }
    //! Returns the total number of block processed
    double get_total_blocks() { return (errors+corrects); }
    //! Returns the block error rate
    double get_errorrate() { return (errors / (corrects + errors)); }

    //protected:
  private:
    bool setup_done;
    int blocksize;
    double errors;
    double corrects;
    bool CORR;
  };

} // namespace itpp

#endif // #ifndef ERROR_COUNTERS_H
