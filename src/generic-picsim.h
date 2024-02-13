/**************************************************************************
 *
 * $Id: generic-picsim.h,v 1.48 2011/03/26 15:36:08 patrick Exp $
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

/*
 * Abstract Base Class GenericPicSim
 *
 * The Derived Classes should provide at least
 * the following member functions
 *
 * 		computeSources()
 * 		solveFields()
 * 		processDiagnostics()
 */

#ifndef GENERIC_PICSIM_H
#define GENERIC_PICSIM_H

#include <species-handler.h>
#include <surface-handler.h>

namespace picsim {

  class GenericPicSim : public parser::Parser {
  public:

    typedef std::ostream ostream;

    typedef pdf::NumberGenerator<real> NumberGenerator;

    GenericPicSim(int nargs, char *args[]);
    virtual ~GenericPicSim();

    virtual bool parseHelp() const;
    virtual bool parseVersion() const;
    virtual bool parseTemplate() const;

    void virtual initialise();

    bool integrateOneStepForward();

  protected:

    typedef parser::Parser Parser;
    typedef Parser::LUT LUT;
    typedef Parser::LUTPair Pair;

    NumberGenerator draw;
    Scheduler       scheduler;
    SpeciesHandler  species;
    SurfaceHandler  surfaces;

    // Particles macro state
    VecReal numParticles, kineticEnergy, kineticTemp;

    virtual void printOn(ostream& os) const;

    virtual void computeSources() = 0;
    virtual void solveFields() = 0;
    virtual void processDiagnostics() = 0;

  private:

  };

}

#endif // GENERIC_PICSIM_H

