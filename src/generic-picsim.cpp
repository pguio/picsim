/**************************************************************************
 *
 * $Id: generic-picsim.cpp,v 1.60 2011/11/10 10:08:07 patrick Exp $
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

#include <generic-picsim.h>

namespace picsim {

  using parser::header;

  GenericPicSim::GenericPicSim(int nargs, char *args[]) :
    Parser(nargs, args),
    draw(nargs, args), scheduler(nargs, args),
    species(nargs, args), surfaces(nargs, args),
    numParticles(0), kineticEnergy(0), kineticTemp(0)
  {
    parseLevelDebugOption("GenericPicSim:dl");
  }

  GenericPicSim::~GenericPicSim ()
  {}

#if defined(HAVE_MPI)

#define PARSE(FUNC)                                                 \
bool GenericPicSim::FUNC() const                                    \
{                                                                   \
  bool parsed = false;                                              \
  if (rankProc == masterProc) {                                     \
    if ( (parsed=Parser::FUNC()) ) {                                \
			draw.FUNC();                                                  \
      scheduler.FUNC();                                             \
		}                                                               \
	}                                                                 \
	MPI_Bcast((void *)&parsed, 1, MPI_BYTE, mpiRoot, MPI_COMM_WORLD); \
	if (parsed) {                                                     \
		species.FUNC();                                                 \
    surfaces.FUNC();                                                \
	}                                                                 \
	return parsed;                                                    \
}

#else

#define PARSE(FUNC)                                                 \
bool GenericPicSim::FUNC() const                                    \
{                                                                   \
	if ( Parser::FUNC() ) {                                           \
		draw.FUNC();                                                    \
		scheduler.FUNC();                                               \
		species.FUNC();                                                 \
		surfaces.FUNC();                                                \
		return true;                                                    \
	}                                                                 \
	return false;                                                     \
}

#endif

  PARSE(parseHelp)
  PARSE(parseVersion)
  PARSE(parseTemplate)


  void GenericPicSim::initialise()
  {
    draw.initialise();
    scheduler.initialise();
    species.initialise(draw,scheduler);
    surfaces.initialise(scheduler);
    surfaces.processSpecies(species);
#if defined(HAVE_MPI)

    if (rankProc == masterProc) {
      int numSpecies = species.getSpeciesNumber();
      numParticles.resize(numSpecies);
      kineticEnergy.resize(numSpecies);
      kineticTemp.resize(numSpecies);
    }
#else
    int numSpecies = species.getSpeciesNumber();
    numParticles.resize(numSpecies);
    kineticEnergy.resize(numSpecies);
    kineticTemp.resize(numSpecies);
#endif

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT1("GenericPicSim::initialise")
    VAR_DEBUG_OUTPUT2(numParticles.size(),kineticEnergy.size())
    END_DEBUG_OUTPUT
  }

  bool GenericPicSim::integrateOneStepForward()
  {
    computeSources();

    solveFields();

    species.moveParticles(draw, scheduler);
    surfaces.processSpecies(species);
    species.getState(numParticles, kineticEnergy,kineticTemp);

    processDiagnostics();

    scheduler.update();

    return scheduler.isEndOfSimulation();
  }

  void GenericPicSim::printOn(ostream& os) const
  {
    const string hstr(DIMRSTR "R" DIMVSTR "V PIC simulation Version " VERSION);

#if defined(HAVE_MPI)
    if (rankProc == masterProc) {

      Parser::printCmd(os);

      os << header(hstr, 1, 1)
         << draw << std::endl
         << scheduler << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    os << species << std::endl
       << surfaces << std::endl;
#else

    Parser::printCmd(os);

    os << header(hstr, 1, 1)
       << draw << '\n'
       << scheduler << '\n'
       << species << '\n'
       << surfaces << '\n';
#endif
  }

}
