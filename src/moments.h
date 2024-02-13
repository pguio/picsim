/**************************************************************************
 *
 * $Id: moments.h,v 1.14 2011/03/26 15:36:08 patrick Exp $
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

#ifndef MOMENTS_H
#define MOMENTS_H

#include <species-handler.h>

namespace picsim {

  class Moments : public parser::Parser {
  public:

    typedef std::ostream ostream;

    typedef Species::Coordinate              Coordinate;
    typedef std::list<const Coordinate *>    ParticlePList;
    typedef ParticlePList::const_iterator    ParticlePListConstIter;

    typedef Species::ParticleListConstIter   ParticleListConstIter;
    typedef SpeciesHandler::SpeciesListIter  SpeciesListIter;

    typedef blitz::Array<ParticlePList,DIMR> ParticlePListArray;

    typedef mudfas::Field Field;

    typedef blitz::TinyVector<Field,DIMV>    VectorField;


    friend ostream &operator<<(ostream &os, const Moments &m);

    Moments(int nargs, char *args[], const string name="");
    ~Moments();

    void initialise(const Scheduler &scheduler);

    void computeMoments(SpeciesHandler &species, Field &dens,
                        VectorField &vel, VectorField &Temp);

    const string &name() const;

    bool isVel(int d)  const;
    bool isTemp(int d) const;

    mudfas::real getStart(int d)  const;
    mudfas::real getEnd(int d)    const;
    mudfas::real getStride(int d) const;

    RVectori dim()      const;
    int      dim(int d) const;

  protected:

    typedef parser::Parser Parser;

  private:

    enum parser_enum { _speciesNames=1,
                       _start, _end, _stride,
                       _vel, _Temp
                     };

    const string momentName;

    NameList speciesNames;

    RVectorr start, end, stride;

    VVectorb uFlag, TFlag;

    Field n;
    real rhoFactor;

    VectorField u, T;

    ParticlePListArray grid;
    RVectori gridSize;

    void computeMoments(Species &species);
    bool isSpeciesinMoments(const string &name) const;

    void initParsing(int nargs, char *args[]);
    void paramParsing();
  };

}

#endif // MOMENTS_H
