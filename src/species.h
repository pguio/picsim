/**************************************************************************
 *
 * $Id: species.h,v 1.129 2011/03/26 15:36:08 patrick Exp $
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
 *  Abstract Base Class Species
 *
 * 	The kinetic is described by methods in the ABC
 *
 * 	The Derived Classes should provide at least
 * 	the following member function
 *
 * 	  injectParticles()
 *
 */

#ifndef SPECIES_H
#define SPECIES_H

#include <scheduler.h>
#include <pdf.h>

namespace picsim {

  class Species : public parser::Parser {
  public:

    typedef std::fstream fstream;
    typedef std::ostream ostream;
    typedef std::ostringstream ostringstream;
    typedef std::string string;

    typedef blitz::Range Range;

    typedef mudfas::Field Field;
    typedef mudfas::iField iField;
    typedef mudfas::VectorField VectorField;

    typedef mudfas::Array2dr Array2dr;
    typedef mudfas::Array3dr Array3dr;

    struct Coordinate {
      RVectorr R; // Space variable
      VVectorr V; // Velocity variable
      Coordinate () : R(0.0), V(0.0) {}
    };

    typedef std::list<Coordinate>        ParticleList;
    typedef ParticleList::iterator       ParticleListIter;
    typedef ParticleList::const_iterator ParticleListConstIter;

    typedef pdf::NumberGenerator<real> NumberGenerator;

    enum bc_enum { periodic=0, insulated, free_space };
    enum axis_enum { x=0, y, z };

    friend ostream& operator<<(ostream& os, const Species &s);

    Species(int nargs=0, char *args[]=0, const string name="");
    virtual ~Species();

    virtual void initialise(NumberGenerator &draw, Scheduler &scheduler);

    void depositCharge(Scheduler &scheduler, Field &Rho) const;
    void moveParticles(NumberGenerator &draw, Scheduler &scheduler);
    void sync(int numLost, real v2Lost);

    virtual void injectParticles(NumberGenerator &draw, Scheduler &scheduler)=0;

    void dump() const;
    void restore();

    const string &name() const;

    int getParticleNumber() const;
    real getKineticEnergy() const;

    ParticleListIter      begin();
    ParticleListIter      end();

    ParticleListConstIter begin() const;
    ParticleListConstIter end()   const;
    unsigned              size()  const;

    ParticleListIter erase(ParticleListIter &i);

    real getChargeNumber() const;
    real getmu_s()         const;

    const Boundary<real,DIMR> getAlfa() const;

  protected:

    typedef parser::Parser Parser;
    typedef Parser::LUT LUT;
    typedef Parser::LUTPair Pair;

    enum parser_enum {
      _mass=1, _charge,
      _lambda, _drift, _temp,
      _eqStart, _eqTemp,
      _bcmin, _bcmax, _rmin, _rmax,
#if defined(HAVE_MPI)
      _cpuIds,
#endif
      next
    };

#if defined(HAVE_MPI)

    VecInt cpuIds;
#endif

    const string speciesName;

    // Particle properties
    real mass, charge;

    // Species properties
    real lambda;
    VVectorr drift, temp;

    bool  eqStart;
    real eqTemp;

    Boundary<int,DIMR>   BC;
    Boundary<real,DIMR> domain;
    Boundary<real,DIMR> alfa;

    parser::Parser::LUT bcMap;

    // Normalised variable
    real         m_s, Z_s;       //  mass, charge
    VVectorr     u_s, T_s, v_s;  //  drift, temp, thermal speed
    real         Teq;            //  hydro start temperature
    VVectorr     T, S;           //  decomposition of rotation B field
    real         dt, dtOver2;    //  time increment

    real         Z_sOverM_s;

    int          numParticles; 		// current number of particles
    int          numToInject;     // at each iteration
    ParticleList particles;	      // List of particles

    Boundary<int,DIMR> Loss;

    real         sumv2Before, sumv2After;
    real         kineticEnergy;

    virtual void printOn(ostream& os) const;

    virtual void fixBoundaryParticles(Scheduler &scheduler);

    void fixBoundaryDensity(Scheduler &scheduler, Field &Rho) const;

  private:

#if 0
#if (DIMR==2)

    template<class RVector>
    Array2dr depositLinear(const RVector &dr) const;
#elif (DIMR==3)
    template<class RVector>
    Array3dr depositLinear(const RVector &dr) const;
#endif
#else
    template<class RVector>
    Field depositLinear(const RVector &dr) const;
#endif

    template<class RVector>
    real interpolateLinear(const Field A, const RVector &dr) const;

    void move0(Scheduler &scheduler);
    void move1(Scheduler &scheduler);
    void move2(Scheduler &scheduler);
    void move3(Scheduler &scheduler);
    void move4(Scheduler &scheduler);
    void move5(Scheduler &scheduler);
    void move6(Scheduler &scheduler);
    void move7(Scheduler &scheduler);

    void (Species::*move[8])(Scheduler &scheduler);

    void deccel0(Scheduler &scheduler);
    void deccel1(Scheduler &scheduler);
    void deccel2(Scheduler &scheduler);
    void deccel3(Scheduler &scheduler);
    void deccel4(Scheduler &scheduler);
    void deccel5(Scheduler &scheduler);
    void deccel6(Scheduler &scheduler);
    void deccel7(Scheduler &scheduler);

    void (Species::*deccel[8])(Scheduler &scheduler);

    VVectorr interpolateField(VectorField &field, RVectorr &pos,
                              Scheduler &scheduler) const;

    template<class RVector>
    void checkBoundary(int &particle_index, ParticleListConstIter &iter,
                       RVector &rindex, RVectori &index,
                       RVectori &GridSize) const;

    void initParsing (int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };

}

#endif // SPECIES_H

