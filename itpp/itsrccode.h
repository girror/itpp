/*!
 * \file
 * \brief Include file for the IT++ source coding module
 * \author Tony Ottosson
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

#ifndef ITSRCCODE_H
#define ITSRCCODE_H

/*!
 * \defgroup srccode Source Coding Module
 * @{
 */

//! \defgroup audio Audio
//! \defgroup image Image Functions and Classes
//! \defgroup lpc LPC-related Functions
//! \defgroup sourcecoding Source Coding Routines

/*!
 * @}
 */

#include <itpp/itsignal.h>
#include <itpp/srccode/audiofile.h>
#include <itpp/srccode/gmm.h>
#include <itpp/srccode/lpcfunc.h>
#include <itpp/srccode/vq.h>
#include <itpp/srccode/vqtrain.h>

#endif // #ifndef ITSRCCODE_H
