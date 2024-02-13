/**************************************************************************
 *
 * $Id: probes.cpp,v 1.34 2011/03/26 15:36:08 patrick Exp $
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

#include <probes.h>

namespace picsim {

  using std::ostream;

#define ID "$Id: probes.cpp,v 1.34 2011/03/26 15:36:08 patrick Exp $"

  ostream &operator<<(ostream &os, const Probes &p)
  {
    return os << "Probes name = " << p.probesName << '\n'
           << "Indices     = (" << p.irange << ")\n"
           << "Coordinates = (" << p.rrange << ")\n";
  }


  Probes::Probes(int nargs, char *args[], const string name)
    : Parser(nargs, args), probesName(name)
  {
    initParsing(nargs, args);
    paramParsing();
  }

  Probes::~Probes()
  {}

  void Probes::initialise(Scheduler &scheduler)
  {
    RVectorr a(scheduler.getDomainSize()/(scheduler.getGridSize()-1));
    RVectorr b(scheduler.getDomainMinBoundary());

    for (int d=0; d<DIMR; ++d) {
      rrange(d).start( a(d)*irange(d).start() +b(d));
      rrange(d).stride(a(d)*irange(d).stride()     );
      rrange(d).end(   a(d)*irange(d).end()   +b(d));
    }
  }

  void Probes::initParsing(int nargs, char *args[])
  {
    registerClass(probesName);
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    string prefix(probesName);
    prefix += "::";
    setPrefix(prefix.c_str());

    parseLevelDebugOption("dl");

#if 0
    const char ivect[] = "int[3]";
#endif

    using parser::types::intVect;

    insertOption(_ix, "x", intVect, "Range spec for x coordinates", Any(irange(0).data));
    insertOption(_jy, "y", intVect, "Range spec for y coordinates", Any(irange(1).data));
#if (DIMR==3)

    insertOption(_kz, "z", intVect, "Range spec for z coordinates", Any(irange(2).data));
#endif
  }

  void Probes::paramParsing()
  {
    parseOption(_ix, irange(0).data);
    parseOption(_jy, irange(1).data);
#if (DIMR==3)

    parseOption(_kz, irange(2).data);
#endif
  }


#undef ID

}
