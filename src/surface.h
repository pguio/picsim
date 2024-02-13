/**************************************************************************
 *
 * $Id: surface.h,v 1.26 2017/09/26 17:12:32 patrick Exp $
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

#ifndef SURFACE_H
#define SURFACE_H

#include <species-handler.h>
#include <mudfas-defs.h>
#if !defined(BLITZ2)
#include <blitz/tinyvec-et.h>
#endif

namespace picsim {

  class Surface : public parser::Parser {
  public:

    typedef std::ostream ostream;
    typedef std::string string;

    typedef mudfas::Field Field;

    typedef SpeciesHandler::SpeciesListIter SpeciesListIter;
    typedef Species::ParticleListIter       ParticleListIter;

    friend ostream& operator<<(ostream& os, const Surface &s);

    Surface(int nargs, char *args[], const string name="");
    ~Surface();

    void initialise(const Scheduler &scheduler);

    const string &name() const {
      return surfaceName;
    }

    void processSpecies(SpeciesHandler &species);
    void processDensity(Field &rho);

  protected:

    typedef parser::Parser Parser;
    typedef Parser::LUT LUT;
    typedef Parser::LUTPair Pair;

    enum type_enum { cartesian, spherical, cylindrical };
    enum dir_enum { xDir=0, yDir, zDir };
    enum property_enum { absorbing, emitting };

  private:

    enum parser_enum { _type=1, _cylDir, _property,
                       _surfRho, _rmin, _rmax
                     };

    const string surfaceName;

    int type;
    int cylDir;
    int property;
    real surfRho;
    Boundary<real,DIMR> limits;

    LUT typeMap;
    LUT dirMap;
    LUT propertyMap;

    RVectorr center;
    RVectorr scale;

    Field rhos;

    void processParticles(Species &species);

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };

}

#endif
