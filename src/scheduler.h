/**************************************************************************
 *
 * $Id: scheduler.h,v 1.93 2011/11/07 18:40:34 patrick Exp $
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


#ifndef _SCHEDULER_H
#define _SCHEDULER_H


#if defined(HAVE_CONFIG_H)
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#endif
#include <picsim-defs.h>
#include <mudfas.h>
#include <time-utils.h>
#if !defined(BLITZ2)
#include <blitz/tinyvec-et.h>
#endif

namespace picsim {

  class Scheduler : public parser::Parser {
  public:

    typedef std::ostream ostream;
    typedef blitz::Range Range;
    typedef mudfas::Field Field;
    typedef mudfas::VectorField VectorField;

    enum gradient_type { momentum_conserving, energy_conserving };

    friend ostream& operator<<(ostream& os, const Scheduler &s);

    Scheduler(int nargs, char *args[]);
    ~Scheduler();

    void initialise();
    void update();

    bool isEforce() const {
      return Eforce;
    }
    bool isBforce() const {
      return Bforce;
    }
    bool isGforce() const {
      return Gforce;
    }
    int  getEMode() const {
      return EMode;
    }

    real getRefMass()   const {
      return RefMass;
    }
    real getRefLambda() const {
      return RefLambda;
    }
    real getRefTemp()   const {
      return RefTemperature;
    }

    Boundary<real,DIMR> getDomainBoundary() const {
      return Domain;
    }
    RVectori getGridSize()          const {
      return GridSize;
    }
    RVectorr getDomainMinBoundary() const {
      return Domain.min;
    }
    RVectorr getDomainMaxBoundary() const {
      return Domain.max;
    }
    RVectorr getDomainSize()        const {
      return Domain.max-Domain.min;
    }

    real getUnitVolume() const {
      return UnitVolume;
    }
    real getVolume()     const {
      return Volume;
    }

    real getTimeInc()       const {
      return TimeInc;
    }
    real getCurrentTime()   const {
      return CurrentTime;
    }
    int   getCurrentIter()   const {
      return CurrentIter;
    }
    int   getMaxIter()       const {
      return MaxIter;
    }

    int   getForceSetUp()    const {
      return ForceSetUp;
    }

    bool  isEndOfSimulation() const {
      return (CurrentIter==MaxIter+1);
    }


    void setInternElectricField(Field &Phi, Boundary<int,DIMR> &BC);
    void getInternElectricField(VectorField &E) const;

    VVectorr getExternElectricField() const;
    VVectorr getExternMagneticField() const;
    VVectorr getExternGravitationField() const;

  protected:

    typedef parser::Parser Parser;
    typedef Parser::LUT LUT;
    typedef Parser::LUTPair Pair;

  private:

    enum parser_enum {
      timeinc=1, maxiter, refmass, reflambda, reftemp,
      gridsize, rmin, rmax,
      eforce, emode, evector, bforce, bvector, gforce, gvector
    };

    real TimeInc, CurrentTime;
    int MaxIter, CurrentIter;
    real RefMass, RefLambda, RefTemperature;
    RVectori GridSize;
    Boundary<real,DIMR> Domain;
    bool Eforce;
    RVectorr Eext;
    int EMode;
    VectorField Ei;
    bool Bforce;
    VVectorr Bext;
    bool Gforce;
    RVectorr Gext;

    timeutils::Timer timer;

    LUT emode_map;

    real          Volume;
    real          UnitVolume;
    int           ForceSetUp;

    void fixBoundaryInternElectricField(Field &Phi, Boundary<int,DIMR> &BC);

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };


  inline
  void Scheduler::getInternElectricField(VectorField &E) const
  {
    for (int d=0; d<DIMR; ++d) {
      E(d).resize(Ei(d).shape());
      E(d) = Ei(d);
    }
  }

  inline
  VVectorr Scheduler::getExternElectricField() const
  {
#if (DIMR==2)
    return VVectorr(Eext(0), Eext(1), 0.0);
#elif (DIMR==3)

    return Eext;
#endif
  }

  inline
  VVectorr Scheduler::getExternMagneticField() const
  {
    return Bext;
  }

  inline
  VVectorr Scheduler::getExternGravitationField() const
  {
#if (DIMR==2)
    return VVectorr(Gext(0), Gext(1), 0.0);
#elif (DIMR==3)

    return Gext;
#endif
  }

}

#endif // SCHEDULER_H

