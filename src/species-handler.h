/**************************************************************************
 *
 * $Id: species-handler.h,v 1.33 2011/03/26 15:36:08 patrick Exp $
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

#ifndef SPECIES_HANDLER_H
#define SPECIES_HANDLER_H

#include <species.h>
#include <species-factory.h>
#include <map>

namespace picsim {

  class SpeciesHandler : public parser::Parser {
  public:

    typedef std::string string;
    typedef std::ostream ostream;
    typedef std::vector<SpeciesFactory<Species>::BasePtr> SpeciesList;

    typedef mudfas::Field Field;

    typedef SpeciesList::iterator       SpeciesListIter;
    typedef SpeciesList::const_iterator SpeciesListConstIter;

    typedef pdf::NumberGenerator<real> NumberGenerator;

    friend ostream& operator<<(ostream& os, const SpeciesHandler &h);

    SpeciesHandler(int nargs, char *args[]);
    ~SpeciesHandler();

    virtual bool parseHelp() const;
    virtual bool parseVersion() const;
    virtual bool parseTemplate() const;

    void initialise(NumberGenerator &draw, Scheduler &scheduler);

    void depositCharge(Scheduler &scheduler, Field &Rho);
    void moveParticles(NumberGenerator &draw, Scheduler &scheduler);
    void getState(VecReal &numParticles, VecReal &kineticEnergy,
                  VecReal &kineticTemp);

    Boundary<real,DIMR> getAlfa()  const;

    int getSpeciesNumber()         const;

    SpeciesListIter      begin();
    SpeciesListIter      end();
    SpeciesListConstIter begin()   const;
    SpeciesListConstIter end()     const;

    int size()                     const;

  protected:

    typedef parser::Parser Parser;

  private:

    const int speciesStartKey;

    typedef SpeciesFactory<Species> Factory;
    typedef Factory::IDKeyType IDKeyType;
    typedef std::map<IDKeyType, NameList> SpeciesTypes;

    SpeciesList species;
    int numSpecies;
    SpeciesTypes types;

#if defined(HAVE_MPI)

    typedef std::map<IDKeyType, VecInt>  SpeciesNumCpu;

    SpeciesNumCpu numCpu;

#endif

    real rhoFactor;

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;

  };

}

#endif

