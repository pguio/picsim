/**************************************************************************
 *
 * $Id: spectra-handler.cpp,v 1.11 2011/03/26 15:36:08 patrick Exp $
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

#include <spectra-handler.h>

namespace picsim {

  using parser::header;

  using std::ostream;

#define ID "$Id: spectra-handler.cpp,v 1.11 2011/03/26 15:36:08 patrick Exp $"


  ostream& operator<<(ostream& os, const SpectraHandler &s)
  {
    if ( ! s.empty() ) {
      os << header("Spectra Handler setup");
      SpectraHandler::SpectraListConstIter i=s.begin(), e=s.end();
      for ( ; i!=e; ++i) {
        os << "\n## Spectra " << (**i).name() << '\n' << **i << '\n';
      }
    }
    return os;
  }

  SpectraHandler::SpectraHandler(int nargs, char *args[])
    : Parser(nargs, args)
  {
    initParsing(nargs, args);
    paramParsing();

    for (NameListIter i=spectraNames.begin(), e=spectraNames.end(); i!=e; ++i) {
      insert(nargs, args, *i);
    }
  }

  SpectraHandler::~SpectraHandler()
  {
    for (SpectraListIter i=spectra.begin(), e=spectra.end(); i!=e; ++i) {

      BEGIN_DEBUG_OUTPUT(30)
      HEADER_DEBUG_OUTPUT1("SpectraHandler::~SpectraHandler")
      VAR_DEBUG_OUTPUT1((**i).name())
      END_DEBUG_OUTPUT

      delete *i;
    }
  }

  void SpectraHandler::insert(int nargs, char *args[], const string &name)
  {
    spectra.insert(spectra.end(), new Spectra(nargs, args, name.c_str()) );
  }

  void SpectraHandler::initParsing(int nargs, char *args[])
  {
    registerClass("SpectraHandler");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    parseLevelDebugOption("SpectraHandler::dl");

    using parser::types::stringVect;

    insertOption(_spectra,"Spectra",stringVect, "Vector of Spectra names",Any(spectraNames));
  }

  void SpectraHandler::paramParsing()
  {
    parseOption(_spectra, spectraNames);
  }

#if defined(HAVE_MPI)

#define PARSE(Fun)                                                  \
bool SpectraHandler::Fun() const                                    \
{                                                                   \
  bool parsed = false;                                              \
  if (rankProc == masterProc) {                                     \
    parsed = Parser::Fun();                                         \
  }                                                                 \
  MPI_Bcast((void *)&parsed, 1, MPI_BYTE, mpiRoot, MPI_COMM_WORLD); \
  if (parsed) {                                                     \
    for (int ip=0; ip<nbProc; ++ip) {                               \
      if (ip == rankProc) {                                         \
        SpectraListConstIter  i = spectra.begin();                  \
        (**i).Fun();                                                \
      }                                                             \
      MPI_Barrier(MPI_COMM_WORLD);                                  \
    }                                                               \
    return true;                                                    \
  }                                                                 \
  return false;                                                     \
}

#else

#define PARSE(Fun)                                                  \
bool SpectraHandler::Fun() const                                    \
{                                                                   \
  if ( Parser::Fun() ) {                                            \
    SpectraListConstIter i=spectra.begin(), e=spectra.end();        \
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
