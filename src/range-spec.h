/**************************************************************************
 *
 * $Id: range-spec.h,v 1.5 2011/03/26 15:36:08 patrick Exp $
 *
 * Copyright (c) 2003-2011 Patrick Guio <patrick.guio@gmail.com>
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

#ifndef RANGE_SPEC_H
#define RANGE_SPEC_H

#include <blitz/range.h>
#include <iostream>

namespace picsim {

  template <class P_numtype>
  class RangeSpec {
  public:

    typedef blitz::TinyVector<P_numtype,3> Triplet;

    friend std::ostream& operator<<(std::ostream &os, const RangeSpec & r) {
      if (r.size() == 1)
        return os << r.start();
      else
        return os << r.start() << ":" << r.stride() << ":" << r.end();
    }

    RangeSpec(P_numtype val=0) : data(val)
    {}
    RangeSpec(Triplet val) : data(val)
    {}
    RangeSpec(P_numtype start, P_numtype stride, P_numtype end) :
      data(start, stride, end)
    {}

    P_numtype start()  const {
      return data(0);
    }
    P_numtype stride() const {
      return data(1);
    }
    P_numtype end()    const {
      return data(2);
    }

    void start(P_numtype _x) {
      data(0) = _x;
    }
    void stride(P_numtype _x) {
      data(1) = _x;
    }
    void end(P_numtype _x) {
      data(2) = _x;
    }

    int size()          const {
      return static_cast<int>((data(2)-data(0))/data(1)+1);
    }

    Triplet data;

  };

  inline
  blitz::Range toRange(const RangeSpec<int> &spec)
  {
    return blitz::Range(spec.start(), spec.end(), spec.stride());
  }

}

#endif
