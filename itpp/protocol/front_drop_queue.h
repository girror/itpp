/*!
 * \file
 * \brief Definitions of a Front Drop Queue class
 * \author Anders Persson and Tony Ottosson
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2009  (see AUTHORS file for a list of contributors)
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

#ifndef FRONT_DROP_QUEUE_H
#define FRONT_DROP_QUEUE_H

#include <itpp/protocol/packet.h>
#include <itpp/protocol/events.h>


namespace itpp
{

//! \addtogroup protocol
//@{

//! ADD DOCUMENTATION HERE
#define DEFAULT_MAX_BYTES_IN_QUEUE 24000

//! ADD DOCUMENTATION HERE
class Front_Drop_Queue : public virtual std::queue<Packet*>
{
public:
  //! ADD DOCUMENTATION HERE
  Front_Drop_Queue(const int max_bytes = DEFAULT_MAX_BYTES_IN_QUEUE)  {
    max_bytes_in_queue = max_bytes;
    bytes_in_queue = 0;
    debug = false;
  }

  // TODO destructor
  //  ~FrontDropQueue() { }

  //! ADD DOCUMENTATION HERE
  void set_debug(const bool enable_debug = true) {
    debug = enable_debug;
  }

  //! ADD DOCUMENTATION HERE
  void push(Packet *packet);
  //! ADD DOCUMENTATION HERE
  void pop();

  //! ADD DOCUMENTATION HERE
  void set_max_byte_size(int max_bytes) { max_bytes_in_queue = max_bytes; }
  //! ADD DOCUMENTATION HERE
  int max_byte_size() { return max_bytes_in_queue; }
  //! ADD DOCUMENTATION HERE
  int byte_size() { return bytes_in_queue; }

private:
  int max_bytes_in_queue;
  int bytes_in_queue;
  int debug;
};

//@}

} // namespace itpp

#endif // #ifndef FRONT_DROP_QUEUE_H

