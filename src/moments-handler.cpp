/**************************************************************************
 *
 * $Id: moments-handler.cpp,v 1.15 2011/03/26 15:36:08 patrick Exp $
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

#include <moments-handler.h>

namespace picsim {

  using parser::header;

  using std::ostream;

#define ID "$Id: moments-handler.cpp,v 1.15 2011/03/26 15:36:08 patrick Exp $"

  ostream& operator<<(ostream& os, const MomentsHandler &h)
  {
    if ( ! h.empty() ) {
#if defined(HAVE_MPI)
      if (h.rankProc == masterProc) {
#endif
        os << header("Moments Handler setup");
        MomentsHandler::MomentsListConstIter i=h.begin(), e=h.end();
        for ( ; i!=e; ++i) {
          os << "\n## Moment " << (**i).name() << '\n' << **i << '\n';
        }
#if defined(HAVE_MPI)

      }
#endif

    }
    return os;
  }

  MomentsHandler::MomentsHandler(int nargs, char *args[])
    : Parser(nargs, args)
  {
    initParsing(nargs, args);
    paramParsing();

    for (NameListIter i=momentsNames.begin(), e=momentsNames.end(); i!=e; ++i) {
      insert(nargs, args, *i);
    }
  }

  MomentsHandler::~MomentsHandler()
  {
    for (MomentsListIter i=moments.begin(), e=moments.end(); i!=e; ++i) {

      BEGIN_DEBUG_OUTPUT(30)
      HEADER_DEBUG_OUTPUT1("MomentsHandler::~MomentsHandler")
      VAR_DEBUG_OUTPUT1((**i).name())
      END_DEBUG_OUTPUT

      delete *i;
    }
  }


#if defined(HAVE_MPI)

#define PARSE(Fun)                                                  \
bool MomentsHandler::Fun() const                                    \
{                                                                   \
  if (rankProc == masterProc && Parser::Fun()) {                    \
    MomentsListConstIter i=moments.begin(), e=moments.end();        \
    for ( ; i!=e; ++i) {                                            \
      (**i).Fun();                                                  \
    }                                                               \
    return true;                                                    \
  }                                                                 \
  return false;                                                     \
}

#else

#define PARSE(Fun)                                                  \
bool MomentsHandler::Fun() const                                    \
{                                                                   \
  if ( Parser::Fun() ) {                                            \
    MomentsListConstIter i=moments.begin(), e=moments.end();        \
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


  void MomentsHandler::initialise(const Scheduler & scheduler)
  {
    MomentsListIter i=moments.begin(), e=moments.end();
    for ( ; i!=e; ++i) {
      (**i).initialise(scheduler);
    }
  }

  void MomentsHandler::insert(int nargs, char *args[], const string &name)
  {
    moments.insert(moments.end(), new Moments(nargs, args, name.c_str()) );
  }

  void MomentsHandler::initParsing(int nargs, char *args[])
  {
    registerClass("MomentsHandler");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    parseLevelDebugOption("MomentsHandler::dl");

    using parser::types::stringVect;

    insertOption(_moments, "Moments", stringVect, "Vector of Moments names", Any(momentsNames));
  }

  void MomentsHandler::paramParsing()
  {
    parseOption(_moments, momentsNames);
  }

#undef ID

}
