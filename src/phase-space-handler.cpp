/**************************************************************************
 *
 * $Id: phase-space-handler.cpp,v 1.22 2011/03/26 15:36:08 patrick Exp $
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

#include <phase-space-handler.h>

namespace picsim {

  using parser::header;

  using std::ostream;

#define ID "$Id: phase-space-handler.cpp,v 1.22 2011/03/26 15:36:08 patrick Exp $"

  ostream& operator<<(ostream& os, const PhaseSpaceHandler &h)
  {
    if ( ! h.empty() ) {
#if defined(HAVE_MPI)
      if (h.rankProc == masterProc) {
#endif
        os << header("Phase Space Handler setup");
        PhaseSpaceHandler::PhaseSpaceListConstIter i=h.begin(), e=h.end();
        for ( ; i!=e; ++i) {
          os << "\n## Phase Space " << (**i).name() << '\n' << **i << '\n';
        }
#if defined(HAVE_MPI)

      }
#endif

    }
    return os;
  }

  PhaseSpaceHandler::PhaseSpaceHandler(int nargs, char *args[])
    : Parser(nargs, args)
  {
    initParsing(nargs, args);
    paramParsing();

    for (NameListIter i=spaceNames.begin(), e=spaceNames.end(); i!=e; ++i) {
      insert(nargs, args, *i);
    }
  }

  PhaseSpaceHandler::~PhaseSpaceHandler()
  {
    for (PhaseSpaceListIter i=spaces.begin(), e=spaces.end(); i!=e; ++i) {

      BEGIN_DEBUG_OUTPUT(30)
      HEADER_DEBUG_OUTPUT1("PhaseSpaceHandler::~PhaseSpaceHandler")
      VAR_DEBUG_OUTPUT1((**i).name())
      END_DEBUG_OUTPUT

      delete *i;
    }
  }


#if defined(HAVE_MPI)

#define PARSE(Fun)                                                  \
bool PhaseSpaceHandler::Fun() const                                 \
{                                                                   \
  if (rankProc == masterProc && Parser::Fun()) {                    \
    PhaseSpaceListConstIter i=spaces.begin(), e=spaces.end();       \
    for ( ; i!=e; ++i) {                                            \
      (**i).Fun();                                                  \
    }                                                               \
    return true;                                                    \
  }                                                                 \
  return false;                                                     \
}

#else

#define PARSE(Fun)                                                  \
bool PhaseSpaceHandler::Fun() const                                 \
{                                                                   \
  if ( Parser::Fun() ) {                                            \
    PhaseSpaceListConstIter i=spaces.begin(), e=spaces.end();       \
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


  void PhaseSpaceHandler::initialise()
  {
    PhaseSpaceListIter i=spaces.begin(), e=spaces.end();
    for ( ; i!=e; ++i) {
      (**i).initialise();
    }
  }

  void PhaseSpaceHandler::insert(int nargs, char *args[], const string &name)
  {
    spaces.insert(spaces.end(), new PhaseSpace(nargs, args, name.c_str()) );
  }

  void PhaseSpaceHandler::initParsing(int nargs, char *args[])
  {
    registerClass("PhaseSpaceHandler");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    parseLevelDebugOption("PhaseSpaceHandler::dl");

    using parser::types::stringVect;

    insertOption(_spaces, "PhaseSpace", stringVect, "Vector of Phase Space names", Any(spaceNames));
  }

  void PhaseSpaceHandler::paramParsing()
  {
    parseOption(_spaces, spaceNames);
  }

#undef ID

}
