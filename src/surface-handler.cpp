/**************************************************************************
 *
 * $Id: surface-handler.cpp,v 1.24 2011/03/26 15:36:08 patrick Exp $
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

#include <surface-handler.h>

namespace picsim {

  using parser::header;

  using std::ostream;

#define ID "$Id: surface-handler.cpp,v 1.24 2011/03/26 15:36:08 patrick Exp $"


  ostream& operator<<(ostream& os, const SurfaceHandler &s)
  {
    if ( ! s.empty() ) {
#if defined(HAVE_MPI)
      if (s.rankProc == masterProc) {
#endif
        os << header("Surface Handler setup");
        SurfaceHandler::SurfaceListConstIter i=s.begin(), e=s.end();
        for ( ; i!=e; ++i) {
          os << "\n## Surface " << (**i).name() << '\n' << **i << '\n';
        }
#if defined(HAVE_MPI)

      }
#endif

    }
    return os;
  }

  SurfaceHandler::SurfaceHandler(int nargs, char *args[])
    : Parser(nargs, args)
  {
    initParsing(nargs, args);
    paramParsing();

    for (NameListIter i=surfaceNames.begin(), e=surfaceNames.end(); i!=e; ++i) {
      insert(nargs, args, *i);
    }
  }

  SurfaceHandler::~SurfaceHandler()
  {
    for (SurfaceListIter i=surfaces.begin(), e=surfaces.end(); i!=e; ++i) {

      BEGIN_DEBUG_OUTPUT(30)
      HEADER_DEBUG_OUTPUT1("PhaseSpaceHandler::~PhaseSpaceHandler")
      VAR_DEBUG_OUTPUT1((**i).name())
      END_DEBUG_OUTPUT

      delete *i;
    }
  }

  void SurfaceHandler::insert(int nargs, char *args[], const string &name)
  {
    surfaces.insert(surfaces.end(), new Surface(nargs, args, name.c_str()) );
  }

  void SurfaceHandler::initParsing(int nargs, char *args[])
  {
    registerClass("SurfaceHandler");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    parseLevelDebugOption("SurfaceHandler::dl");

    using parser::types::stringVect;

    insertOption(_surfaces, "Surfaces", stringVect, "Vector of Surfaces names", Any(surfaceNames));
  }

  void SurfaceHandler::paramParsing()
  {
    parseOption(_surfaces, surfaceNames);
  }

#if defined(HAVE_MPI)

#define PARSE(Fun)                                                  \
bool SurfaceHandler::Fun() const                                    \
{                                                                   \
  if (rankProc == masterProc && Parser::Fun()) {                    \
    SurfaceListConstIter i=surfaces.begin(), e=surfaces.end();      \
    for ( ; i!=e; ++i) {                                            \
      (**i).Fun();                                                  \
    }                                                               \
    return true;                                                    \
  }                                                                 \
  return false;                                                     \
}

#else

#define PARSE(Fun)                                                  \
bool SurfaceHandler::Fun() const                                    \
{                                                                   \
  if ( Parser::Fun() ) {                                            \
    SurfaceListConstIter i=surfaces.begin(), e=surfaces.end();      \
    for ( ; i!=e; ++i) {                                            \
      (**i).Fun();                                                  \
    }                                                               \
    return true;                                                    \
  }                                                                 \
  return false;                                                     \
}

#endif


  PARSE(parseHelp)
  PARSE(parseVersion)
  PARSE(parseTemplate)

#undef PARSE

  void SurfaceHandler::initialise(const Scheduler &scheduler)
  {
    for (SurfaceListIter i=surfaces.begin(), e=surfaces.end() ; i!=e; ++i) {
      (**i).initialise(scheduler);
    }
  }


  void SurfaceHandler::processSpecies(SpeciesHandler &species)
  {
    for (SurfaceListIter i=surfaces.begin(), e=surfaces.end() ; i!=e; ++i) {
      (**i).processSpecies(species);
    }
  }


  void SurfaceHandler::processDensity(Field &rho)
  {
    for (SurfaceListIter i=surfaces.begin(), e=surfaces.end() ; i!=e; ++i) {
      (**i).processDensity(rho);
    }
  }


#undef ID

}
