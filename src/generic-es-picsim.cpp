/**************************************************************************
 *
 * $Id: generic-es-picsim.cpp,v 1.14 2017/12/23 18:52:30 patrick Exp $
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


#ifndef GENERIC_ES_PICSIM_CPP
#define GENERIC_ES_PICSIM_CPP

namespace picsim {

  using std::endl;

  template <class Solver>
  GenericEsPicSim<Solver>::GenericEsPicSim(int nargs, char *args[])
    : GenericPicSim(nargs, args), solver(nargs, args), diagnostics(nargs, args)
  {}

  template <class Solver>
  GenericEsPicSim<Solver>::~GenericEsPicSim()
  {}

#if defined(HAVE_MPI)

#define PARSE(Fun)                                      \
template <class Solver>                                 \
bool GenericEsPicSim<Solver>::Fun() const               \
{                                                       \
	bool parsed = GenericPicSim::Fun();                   \
	if (rankProc == masterProc && parsed) {               \
		solver.Fun();                                       \
		diagnostics.Fun();                                  \
	}                                                     \
	return parsed;                                        \
}

#else

#define PARSE(Fun)                                      \
template <class Solver>                                 \
bool GenericEsPicSim<Solver>::Fun() const               \
{                                                       \
	if ( GenericPicSim::Fun() ) {                         \
		solver.Fun();                                       \
		diagnostics.Fun();                                  \
		return true;                                        \
	}                                                     \
	return false;                                         \
}

#endif

  PARSE(parseHelp)
  PARSE(parseVersion)
  PARSE(parseTemplate)

#undef PARSE

  template <class Solver>
  void GenericEsPicSim<Solver>::initialise()
  {
    GenericPicSim::initialise();

    mudfas::FieldVeci gridSize(solver.getGridSize());
    Rho.resize(gridSize);
    Phi.resize(gridSize);

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT1("GenericEsPicSim<Solver>::initialise")
    VAR_DEBUG_OUTPUT1(gridSize)
    VAR_DEBUG_OUTPUT2(Rho.shape(),Phi.shape())
    END_DEBUG_OUTPUT

#if defined(HAVE_MPI)
    if (isEforce()) {
#if defined(SOL1)
      Boundary<real,DIMR> alfa(species.getAlfa());
      if (rankProc == masterProc) {
        solver.initialise();
        //solver.setGbdr(alfa);
      }
#else // SOL2
      Boundary<real,DIMR> alfa(species.getAlfa());
      solver.initialise();
      //solver.setGbdr(alfa);
#endif
    }
#else
    if (isEforce()) {
      solver.initialise();
      Boundary<mudfas::real,DIMR> alfa(species.getAlfa());
      //solver.setGbdr(alfa);
    }
#endif

    initDiagnostics();
  }

  template <class Solver>
  void GenericEsPicSim<Solver>::computeSources()
  {
    species.depositCharge(scheduler, Rho);
    surfaces.processDensity(Rho);
  }

  template<class Solver>
  void GenericEsPicSim<Solver>::printOn(ostream& os) const
  {
    GenericPicSim::printOn(os);

#if defined(HAVE_MPI)

    if (rankProc == masterProc) {
      if ( isEforce() ) {
        os << solver << endl;
      }
      os << diagnostics << endl;
    }
#else
    if ( isEforce() ) {
      os << solver << endl;
    }
    os << diagnostics << endl;
#endif
  }

}

#endif // GENERIC_ES_PICSIM_CPP

