/**************************************************************************
 *
 * $Id: phase-space.h,v 1.31 2011/03/26 15:36:08 patrick Exp $
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

#ifndef PHASE_SPACE_H
#define PHASE_SPACE_H

#include <scheduler.h>
#include <species-handler.h>

namespace picsim {

  class PhaseSpace : public parser::Parser {
  public:

    typedef std::ostream ostream;

    typedef Species::ParticleListConstIter  ParticleListConstIter;
    typedef SpeciesHandler::SpeciesListIter SpeciesListIter;
    typedef SpeciesHandler::SpeciesListConstIter SpeciesListConstIter;

    friend ostream &operator<<(ostream &os, const PhaseSpace &V);

    PhaseSpace(int nargs, char *args[], const string name="");
    ~PhaseSpace();

    void initialise();

    void computePhaseSpace(const Scheduler & scheduler,
                           const SpeciesHandler &species,
                           RVArrayr &space);

    const string & name() const;

    unsigned  numDim()   const;
    RVVectori dim()      const;
    int       dim(int d) const;

    bool isIntegrated(int d) const;

    mudfas::real getStart(int d)  const;
    mudfas::real getEnd(int d)    const;
    mudfas::real getStride(int d) const;

  protected:

    typedef parser::Parser Parser;

  private:

    enum parser_enum { _speciesNames=1,
                       _start, _end, _stride, _integrate, _speed
                     };

    const string spaceName;

    NameList speciesNames;

    RVVectorr start, end, stride;
    RVVectorb integrate;
    RVectorr speed;

    unsigned spaceNumDim;
    RVVectori spaceDim;

    void integrateParticles(real t, const Species &species, RVArrayr &space);
    bool isSpeciesinSpace(const string & name) const;

    void initParsing(int nargs, char *args[]);
    void paramParsing();
  };

}

#endif // PHASE_SPACE_H
