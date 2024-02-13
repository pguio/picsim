/**************************************************************************
 *
 * $Id: probes-handler.cpp,v 1.14 2011/03/26 15:36:08 patrick Exp $
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

#include <probes-handler.h>

namespace picsim {

  using parser::header;

  using std::ostream;


#define ID "$Id: probes-handler.cpp,v 1.14 2011/03/26 15:36:08 patrick Exp $"


  ostream& operator<<(ostream& os, const ProbesHandler &p)
  {
    if ( ! p.empty() ) {
      os << header("Probes Handler setup");
      ProbesHandler::ProbesListConstIter i=p.begin(), e=p.end();
      for ( ; i!=e; ++i) {
        os << "\n## Probes " << (**i).name() << '\n' << **i << '\n';
      }
    }
    return os;
  }

  ProbesHandler::ProbesHandler(int nargs, char *args[])
    : Parser(nargs, args)
  {
    initParsing(nargs, args);
    paramParsing();

    for (NameListIter i=probesNames.begin(), e=probesNames.end(); i!=e; ++i) {
      insert(nargs, args, *i);
    }
  }

  ProbesHandler::~ProbesHandler()
  {
    for (ProbesListIter i=probes.begin(), e=probes.end(); i!=e; ++i) {

      BEGIN_DEBUG_OUTPUT(30)
      HEADER_DEBUG_OUTPUT1("ProbesHandler::~ProbesHandler")
      VAR_DEBUG_OUTPUT1((**i).name())
      END_DEBUG_OUTPUT

      delete *i;
    }
  }

  void ProbesHandler::insert(int nargs, char *args[], const string &name)
  {
    probes.insert(probes.end(), new Probes(nargs, args, name.c_str()) );
  }

  void ProbesHandler::initParsing(int nargs, char *args[])
  {
    registerClass("ProbesHandler");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    parseLevelDebugOption("ProbesHandler::dl");

    using parser::types::stringVect;

    insertOption(_probes, "Probes", stringVect, "Vector of Probes names", Any(probesNames));
  }

  void ProbesHandler::paramParsing()
  {
    parseOption(_probes, probesNames);
  }

#if defined(HAVE_MPI)

#define PARSE(Fun)                                                  \
bool ProbesHandler::Fun() const                                     \
{                                                                   \
  if (rankProc == masterProc  && Parser::Fun()) {                   \
    ProbesListConstIter i=probes.begin(), e=probes.end();           \
    for ( ; i!=e; ++i) {                                            \
      (**i).Fun();                                                  \
    }                                                               \
    return true;                                                    \
  }                                                                 \
  return false;                                                     \
}

#else

#define PARSE(Fun)                                                  \
bool ProbesHandler::Fun() const                                     \
{                                                                   \
  if ( Parser::Fun() ) {                                            \
    ProbesListConstIter i=probes.begin(), e=probes.end();           \
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


#undef ID

}
