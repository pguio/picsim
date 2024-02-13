/**************************************************************************
 *
 * $Id: picsim-defs.h,v 1.58 2011/03/26 15:36:08 patrick Exp $
 *
 * Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
 * All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2.  of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef PICSIM_DEFS_H
#define PICSIM_DEFS_H

#if defined(HAVE_CONFIG_H)
#include <picsim-config.h>
#endif

#include <list>
#include <vector>

#include <blitz/array.h>
#include <blitz/numinquire.h>

#include <range-spec.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#define DEBUG_LEVEL
#include <picsim-debug.h>

namespace picsim {

#define PICSIM_COPYRIGHT \
"Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>\n\n"\
"This is free software; see the source for copying conditions.  There is NO\n"\
"warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"

#if defined(FLOAT_COORDINATE)
  typedef float  real;
#elif defined(DOUBLE_COORDINATE)
  typedef double real;
#else
#error macros FLOAT_COORDINATE or DOUBLE_COORDINATE must be defined
#endif


#if !defined(DIMR)
#error macros DIMR must be defined
#endif
#if !defined(DIMV)
#error macros DIMV must be defined
#endif

#define DIMRV DIMR+DIMV

  typedef blitz::TinyVector<bool,DIMR>       RVectorb;
  typedef blitz::TinyVector<bool,DIMV>       VVectorb;
  typedef blitz::TinyVector<bool,DIMRV>      RVVectorb;

  typedef blitz::TinyVector<int,DIMR>        RVectori;
  typedef blitz::TinyVector<int,DIMR-1>      BVectori;
  typedef blitz::TinyVector<int,DIMV>        VVectori;
  typedef blitz::TinyVector<int,DIMRV>       RVVectori;

  typedef blitz::TinyVector<real,DIMR>      RVectorr;
  typedef blitz::TinyVector<real,DIMV>      VVectorr;
  typedef blitz::TinyVector<real,DIMRV>     RVVectorr;

  typedef blitz::TinyVector<double,DIMR>     RVectord;

  typedef blitz::Array<real,DIMRV>          RVArrayr;

  typedef std::vector<real> VecReal;

  typedef std::vector<int>         VecInt;
  typedef VecInt::iterator         VecIntIter;
  typedef VecInt::const_iterator   VecIntConstIter;

  // Define string DIMRSTR and DIMVSTR
#if (DIMR == 2)
#define DIMRSTR "2"
#define DIMRVSTR "2+3"
#elif (DIMR == 3)
#define DIMRSTR "3"
#define DIMRVSTR "3+3"
#endif

#define DIMVSTR "3"

#if (DIMR == 2)
  const RVectorr DEFAULT_E(0.0, 0.0);
  const VVectorr DEFAULT_B(0.0, 1.0, 0.0);
  const RVectorr DEFAULT_G(0.0, -5.0e-2);
#elif (DIMR == 3)
  const RVectorr DEFAULT_E(0.0, 0.0, 0.0);
  const RVectorr DEFAULT_B(0.0, 1.0, 0.0);
  const RVectorr DEFAULT_G(0.0, -5.0e-2, 0.0);
#endif

  // Define the operator<< template for TinyVector<T,DIM*>
  template<class T>
  std::ostream& operator<<(std::ostream &os, const blitz::TinyVector<T,DIMR> &v)
  {
    for (int d=0; d<DIMR-1; ++d)
      os << v(d) << ", ";
    return os << v(DIMR-1);
  }

#if (DIMR != DIMV)
  template<class T>
  std::ostream& operator<<(std::ostream &os, const blitz::TinyVector<T,DIMV> &v)
  {
    for (int d=0; d<DIMV-1; ++d)
      os << v(d) << ", ";
    return os << v(DIMV-1);
  }
#endif

  template<class T>
  std::ostream& operator<<(std::ostream &os, const blitz::TinyVector<T,DIMRV> &v)
  {
    for (int d=0; d<DIMRV-1; ++d)
      os << v(d) << ", ";
    return os << v(DIMRV-1);
  }

  typedef std::vector<std::string>    NameList;
  typedef NameList::iterator          NameListIter;
  typedef NameList::const_iterator    NameListConstIter;

  template<class T>
  std::ostream& operator<<(std::ostream &os, const std::vector<T> &v)
  {
    if (!v.empty()) {
      typename std::vector<T>::const_iterator i = v.begin();
      typename std::vector<T>::const_iterator end = v.end();
      os << *i;
      for ( ++i; i != end; ++i )
        os << ", " << *i;
    }
    return os;
  }


  typedef blitz::TinyVector<int,3> IterationSpec; // start, stride, end
  const IterationSpec DEFAULT_SaveIter(0, 10, 50);

#if defined(HAVE_MPI)
#if defined(FLOAT_COORDINATE)
#define MPI_COORDINATE MPI_FLOAT
#elif defined(DOUBLE_COORDINATE)
#define MPI_COORDINATE MPI_DOUBLE
#endif

#endif // defined(HAVE_MPI)

}


template<class T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &v)
{
  if (!v.empty()) {
    typename std::vector<T>::const_iterator i = v.begin();
    typename std::vector<T>::const_iterator end = v.end();
    os << *i;
    for ( ++i; i != end; ++i )
      os << ", " << *i;
  }
  return os;
}


#include <parser.h>

#endif // PICSIM_DEFS_H
