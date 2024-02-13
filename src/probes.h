/**************************************************************************
 *
 * $Id: probes.h,v 1.31 2011/03/26 15:36:08 patrick Exp $
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

#ifndef PROBES_H
#define PROBES_H

#include <scheduler.h>

namespace picsim {

  class Probes : public parser::Parser {
  public:

    typedef std::string string;
    typedef std::ostream ostream;

    typedef blitz::Range Range;
    typedef blitz::TinyVector<RangeSpec<int>  , DIMR> RRangeSpeci;
    typedef blitz::TinyVector<RangeSpec<real>, DIMR> RRangeSpecr;

    friend ostream &operator<<(ostream &os, const Probes &p);

    Probes(int nargs, char *args[], const string name="");
    ~Probes();

    void initialise(Scheduler &scheduler);

    const string &name()  const {
      return probesName;
    }

    int   size(int d)     const {
      return irange(d).size();
    }

    real start(int d )   const {
      return rrange(d).start();
    }
    real stride(int d)   const {
      return rrange(d).stride();
    }
    real end(int d)      const {
      return rrange(d).end();
    }

    Range range(int d)    const {
      return toRange(irange(d));
    }

  protected:

    typedef parser::Parser Parser;

  private:

#if (DIMR==2)

    enum parser_enum { _ix=1, _jy };
#elif (DIMR==3)

    enum parser_enum { _ix=1, _jy, _kz };
#endif

    const string probesName;

    RRangeSpeci irange;
    RRangeSpecr rrange;

    void initParsing (int nargs, char *args[]);
    void paramParsing();

  };

}

#endif // PROBES_H
