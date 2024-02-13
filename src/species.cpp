/**************************************************************************
 *
 * $Id: species.cpp,v 1.186 2020/04/22 13:54:13 patrick Exp $
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

#include <species.h>

namespace picsim {

  using blitz::cross;
  using blitz::dot;
  using blitz::epsilon;
  using blitz::floor;
  using blitz::sqrt;

  using blitz::TinyVector;

  using blitz::firstDim;
  using blitz::secondDim;
#if (DIMR==3)

  using blitz::thirdDim;
#endif

  using parser::map_elt;
  using parser::yesno;

  using std::ios;
  using std::ostream;

#define ID "$Id: species.cpp,v 1.186 2020/04/22 13:54:13 patrick Exp $"

  void Species::printOn(ostream& os) const
  {
    os << "mass                    = " << mass << '\n'
       << "charge                  = " << charge << '\n'

       << "lambda                  = " << lambda << '\n'
       << "drift velocity          = " << drift << '\n'
       << "temperature             = " << temp << '\n'

       << "Equilibrium start       = " << yesno(eqStart) << '\n';
    if (eqStart)
      os
          << "Equilibrium temperature = " << eqTemp << '\n';

    os << "rmin                    = " << domain.min << '\n'
       << "rmax                    = " << domain.max << '\n'
       << "bcmin                   = " << map_elt(bcMap,BC.min) << '\n'
       << "bcmax                   = " << map_elt(bcMap,BC.max) << '\n'

       << "Initial part. numb.     = " << numParticles << '\n'
       << "Injected part. numb.    = " << numToInject << '\n';
  }

  ostream& operator<<(ostream& os, const Species &s)
  {
    s.printOn(os);
    return os;
  }

  Species::Species(int nargs, char *args[], const string name)
    : Parser(nargs, args),
#if defined(HAVE_MPI)
      cpuIds(0),
#endif
      speciesName(name), mass(1.0), charge(1.0),
      lambda(20.0), drift(0.0), temp(1.0), eqStart(false), eqTemp(0.0),
      BC(free_space,free_space), domain(0.0,mudfas::DEFAULT_GRIDSIZE-1), alfa(0.0,0.0),
      numParticles(0), numToInject(0), particles(0)
  {
    bcMap.insert(Pair(periodic  , "  periodic"));
    bcMap.insert(Pair(insulated , " insulated"));
    bcMap.insert(Pair(free_space, "free space"));

    initParsing(nargs, args);
    paramParsing();

    checkParam();

    move[0] = &Species::move0;
    move[1] = &Species::move1;
    move[2] = &Species::move2;
    move[3] = &Species::move3;
    move[4] = &Species::move4;
    move[5] = &Species::move5;
    move[6] = &Species::move6;
    move[7] = &Species::move7;

    deccel[0] = &Species::deccel0;
    deccel[1] = &Species::deccel1;
    deccel[2] = &Species::deccel2;
    deccel[3] = &Species::deccel3;
    deccel[4] = &Species::deccel4;
    deccel[5] = &Species::deccel5;
    deccel[6] = &Species::deccel6;
    deccel[7] = &Species::deccel7;
  }

  Species::~Species()
  {}


  void Species::initialise(NumberGenerator &draw, Scheduler &scheduler)
  {
    m_s = mass/scheduler.getRefMass();
    Z_s = charge;

    Z_sOverM_s = Z_s/m_s;

    u_s = drift/std::sqrt(scheduler.getRefTemp()/scheduler.getRefMass());
    T_s = temp/scheduler.getRefTemp();
    v_s = blitz::sqrt(T_s/m_s);

    Teq = eqTemp/scheduler.getRefTemp();

    // Rotation parameters
    // From Plasma physics via computer simulation
    // Birsdall and Langdon, p. 62 eq. 11 and 13
    VVectorr B(scheduler.getExternMagneticField());
    dt = scheduler.getTimeInc();
    dtOver2 = dt/2.0;
    T = -Z_s/m_s*dtOver2*B;
#if !defined(_AIX)

    S = 2.0*T/(1.0+dot(T,T));
#else

    S = T/(1.0+dot(T,T));
    S *= 2.0;
#endif

  }

  template <class RVector>
  void Species::checkBoundary(int &particleIndex, ParticleListConstIter &i,
                              RVector &r, RVectori &ir,
                              RVectori &gridSize) const
  {

    if (any(ir > gridSize-2))  {
      ostringstream os;
      os << "\nany(ir > gridSize-2) = " << any(ir > gridSize-2)
         << "\nparticle ir = " << particleIndex
         << "\nR = " << i->R
         << "\nr = " << r
         << "\nir = " << ir;
      throw ClassException("Species", os.str());
    } else if (any(ir < 0))  {
      ostringstream os;
      os << "\nany(ir < 0) = " << any(ir < 0)
         << "\nparticle ir = " << particleIndex
         << "\nR = " << i->R
         << "\nr = " << r
         << "\nir = " << ir;
      throw ClassException("Species", os.str());
    } else {
      ++particleIndex;
    }
  }

#if (DIMR==2)
  template <class RVector>
  Species::Field Species::depositLinear(const RVector &dr) const
  {
    Field xy(2,2);
    const RVector cdr(1.0-dr);
    xy =
      cdr(0)*cdr(1), cdr(0)*dr(1),
      dr(0)*cdr(1),  dr(0)*dr(1);
    return xy;
  }

  template <class RVector>
  real Species::interpolateLinear(const Field A, const RVector &dr) const
  {
    const RVector cdr(1.0-dr);
    TinyVector<double,4> a(A(0,0), A(0,1),
                           A(1,0), A(1,1));
    TinyVector<double,4> b(cdr(0)*cdr(1), cdr(0)*dr(1),
                           dr(0)*cdr(1),  dr(0)*dr(1));
    return static_cast<real>(dot(a,b));
  }
#elif (DIMR==3)
  template <class RVector>
  Species::Field Species::depositLinear(const RVector &dr) const
  {
    Field xyz(2,2,2);
    const RVector cdr(1.0-dr);
    xyz =
      cdr(0)*cdr(1)*cdr(2), cdr(0)*cdr(1)*dr(2),
      cdr(0)* dr(1)*cdr(2), cdr(0)* dr(1)*dr(2),
      dr(0)*cdr(1)*cdr(2),  dr(0)*cdr(1)*dr(2),
      dr(0)* dr(1)*cdr(2),  dr(0)* dr(1)*dr(2);
    return xyz;
  }

  template <class RVector>
  real Species::interpolateLinear(const Field A, const RVector &dr) const
  {
    const RVector cdr(1.0-dr);
    TinyVector<double,8> a(A(0,0,0), A(0,0,1),
                           A(0,1,0), A(0,1,1),
                           A(1,0,0), A(1,0,1),
                           A(1,1,0), A(1,1,1));
    TinyVector<double,8> b(cdr(0)*cdr(1)*cdr(2), cdr(0)*cdr(1)*dr(2),
                           cdr(0)* dr(1)*cdr(2), cdr(0)* dr(1)*dr(2),
                           dr(0)*cdr(1)*cdr(2),  dr(0)*cdr(1)*dr(2),
                           dr(0)* dr(1)*cdr(2),  dr(0)* dr(1)*dr(2));
    return static_cast<real>(dot(a,b));
  }
#endif

#if 1

  void Species::depositCharge(Scheduler &scheduler, Field &Rho) const
  {
#if defined(BZ_DEBUG)
    int particleIndex = 0;
#endif

    RVectori gridSize(Rho.shape());
    RVectord a((gridSize-1)/scheduler.getDomainSize());
    RVectord domainMin(scheduler.getDomainMinBoundary());

    Rho = 0.0;

    ParticleListConstIter i=particles.begin(), e=particles.end();
    for ( ; i!=e; ++i) {
      RVectord r(a*(i->R-domainMin));
      RVectori ir(floor(r));
      RVectord weight(r-ir);
      for (int d=0; d<DIMR; ++d) {
        // precision loss due to affine transform [a,b[ -> [0,gridSize-1]
        if (ir(d) > gridSize(d)-2) {
          ir(d) = gridSize(d)-2;
          weight(d) = 1.0;
        }
      }
#if defined(BZ_DEBUG)
      checkBoundary(particleIndex, i, r, ir, gridSize);
#endif

      Range I(ir(0),ir(0)+1);
      Range J(ir(1),ir(1)+1);
#if (DIMR==2)

      Rho(I,J) += depositLinear(weight);
#elif (DIMR==3)

      Range K(ir(2),ir(2)+1);
      Rho(I,J,K) += depositLinear(weight);
#endif

    }
    Rho *= Z_s;
    fixBoundaryDensity(scheduler, Rho);
  }

#else

  void Species::depositCharge(Scheduler &scheduler, Field &Rho) const
  {
#if defined(BZ_DEBUG)
    int particleIndex = 0;
#endif

    RVectori gridSize(Rho.shape());
    RVectord a((gridSize-1)/scheduler.getDomainSize());
    RVectord domainMin(scheduler.getDomainMinBoundary());

    Field rho1(gridSize+1);
    rho1 = 0.0;

    ParticleListConstIter i=particles.begin(), e=particles.end();
    for ( ; i!=e; ++i) {
      RVectord r(a*(i->R-domainMin));
      RVectori ir(floor(r));
      RVectord weight(r-ir);
#if defined(BZ_DEBUG)

      checkBoundary(particleIndex, i, r, ir, gridSize);
#endif

      Range I(ir(0),ir(0)+1);
      Range J(ir(1),ir(1)+1);
#if (DIMR==2)

      rho1(I,J) += depositLinear(weight);
#elif (DIMR==3)

      Range K(ir(2),ir(2)+1);
      rho1(I,J,K) += depositLinear(weight);
#endif

    }
#if (DIMR==2)
    Rho = Z_s*rho1(Range(0,gridSize(0)-1),Range(0,gridSize(1)-1));
#elif (DIMR==3)

    Rho = Z_s*rho1(Range(0,gridSize(0)-1),Range(0,gridSize(1)-1),
                   Range(0,gridSize(2)-1));
#endif

    fixBoundaryDensity(scheduler, Rho);
  }

#endif

#if (DIMR==2)
  void Species::fixBoundaryDensity(Scheduler &scheduler, Field &Rho) const
  {
    Range all(Range::all());
    RVectori gridsize(Rho.shape());
    RVectorr dr(scheduler.getDomainSize()/(gridsize-1));
    for (int d=0; d<DIMR; ++d) {
      // first dimension is always accessed thanks to the
      // permutation trick
      int end(Rho.ubound(0));
      switch (BC.min(d)) {
      case periodic:
        Rho(0,all) = Rho(0,all)+Rho(end,all);
        Rho(end,all) = Rho(0,all);
        break;
      case insulated:
        Rho(0,all) = 1.0;
        break;
      case free_space:
        Rho(0,all) *= (2.0-alfa.min(d)*dr(d));
        break;
      }
      switch (BC.max(d)) {
      case periodic:
        break;
      case insulated:
        Rho(end,all) = 1.0;
        break;
      case free_space:
        Rho(end,all) *= (2.0+alfa.max(d)*dr(d));
        break;
      }
      // trick to permute on the left the dimensions
      // which let access always the same indices
      Rho.transposeSelf(secondDim, firstDim);
    }
  }
#elif (DIMR==3)
  void Species::fixBoundaryDensity(Scheduler &scheduler, Field &Rho) const
  {
    Range all(Range::all());
    RVectori gridsize(Rho.shape());
    RVectorr dr(scheduler.getDomainSize()/(gridsize-1));
    for (int d=0; d<DIMR; ++d) {
      // first dimension is always accessed thanks to the
      // permutation trick
      int end(Rho.ubound(0));
      switch (BC.min(d)) {
      case periodic:
        Rho(0,all,all) = Rho(0,all,all)+Rho(end,all,all);
        Rho(end,all,all) = Rho(0,all,all);
        break;
      case insulated:
        Rho(0,all,all) = 1.0;
        break;
      case free_space:
        Rho(0,all,all) *= (2.0-alfa.min(d)*dr(d));
        break;
      }
      switch (BC.max(d)) {
      case periodic:
        break;
      case insulated:
        Rho(end,all,all) = 1.0;
        break;
      case free_space:
        Rho(end,all,all) *= (2.0+alfa.max(d)*dr(d));
        break;
      }
      // trick to permute on the left the dimensions
      // which let access always the same indices
      Rho.transposeSelf(secondDim, thirdDim, firstDim);
    }
  }
#endif

  VVectorr Species::interpolateField(VectorField &field, RVectorr &pos,
                                     Scheduler &scheduler) const
  {
    if (!scheduler.isEforce())
      return VVectorr(0.0,0.0,0.0);

    RVectord delta_r(scheduler.getDomainSize());
    RVectord domainMin(scheduler.getDomainMinBoundary());
    RVectord domainMax(scheduler.getDomainMaxBoundary());
    RVectori gridSize(field(0).shape());
    RVectord r((pos-domainMin)*(gridSize-1)/delta_r);
    RVectori ir;
#if (DIMR==2)
    RVectord weight(0.0,0.0);
    Range I, J;
#elif (DIMR==3)
    RVectord weight(0.0,0.0,0.0);
    Range I, J, K;
#endif

    switch (scheduler.getEMode()) {
    case Scheduler::momentum_conserving:
      ir = floor(r);
      weight = r-ir;
      for (int d=0; d<DIMR; ++d) {
        // precision loss due to affine transform [a,b[ -> [0,gridSize-1]
        if (ir(d) > gridSize(d)-2) {
          ir(d) = gridSize(d)-2;
          weight(d) = 1.0;
        }
      }
#if defined(BZ_DEBUG)
      if (any(pos == domainMax) || any(ir == gridSize-1)) {
        ostringstream os;
        os << "\nany(pos == domainMax) = " << RVectord(any(pos == domainMax))
           << "\npos = " << pos
           << "\ndomainMax = " << domainMax
           << "\nany(ir == gridSize-1) = " << RVectori(any(ir == gridSize-1))
           << "\nir = " << ir
           << "\nr = " << r
           << "\ngridSize-1 = " << RVectori(gridSize-1);
        throw ClassException("Species", os.str());
      }
#endif
      I.setRange(ir(0),ir(0)+1);
      J.setRange(ir(1),ir(1)+1);
#if (DIMR==3)

      K.setRange(ir(2),ir(2)+1);
#endif

      break;
    case Scheduler::energy_conserving:
      // well-behaved for periodic boundary systems
      r -= 0.5;
      ir = floor(r);
      weight = r-ir;

#define ADJUST_RANGE(dim,var)                        \
if (r(dim) < 0.0)                                    \
	var.setRange(gridSize(dim)-2,0,-gridSize(dim)+2);  \
else if (r(dim) >= gridSize(dim)-2)                  \
	var.setRange(0,gridSize(dim)-2,gridSize(dim)-2);   \
else                                                 \
	var.setRange(ir(dim),ir(dim)+1);

      ADJUST_RANGE(0,I)
      ADJUST_RANGE(1,J)
#if (DIMR==3)
      ADJUST_RANGE(2,K)
#endif

#undef ADJUST_RANGE
      break;
    }
    // initialise to zero!
    VVectorr fieldpos(0.0,0.0,0.0);
    for (int d=0; d<DIMR; ++d) {
#if (DIMR==2)
      Field f(field(d)(I,J));
#elif (DIMR==3)

      Field f(field(d)(I,J,K));
#endif

      fieldpos(d) = interpolateLinear(f,weight);
    }
    return fieldpos;
  }

  void Species::move0(Scheduler &scheduler)
  // No forces
  {
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // accumulate old velocity
      real v2 = dot(i->V, i->V);
      sumv2Before += v2;
      // accumulate new velocity
      sumv2After += v2;
      // Velocity to Position
      for (int d=0; d<DIMR; ++d)
        i->R(d) += i->V(d) * dt;
    }
  }

  void Species::move1(Scheduler &scheduler)
  // Electrostatic force
  {
    VVectorr E_extern(scheduler.getExternElectricField());
    VectorField E_intern;
    scheduler.getInternElectricField(E_intern);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // accumulate old velocity
      sumv2Before += dot(i->V, i->V);
      // Translation
      VVectorr Edt = (E_extern +
                      interpolateField(E_intern, i->R, scheduler))*dt;
      Edt *= Z_sOverM_s;
      i->V += Edt;
      // accumulate new velocity
      sumv2After += dot(i->V, i->V);
      // Velocity to Position
      for (int d=0; d<DIMR; ++d)
        i->R(d) += i->V(d) * dt;
    }
  }

  void Species::move2(Scheduler &scheduler)
  // Magnetic force
  {
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // accumulate old velocity
      sumv2Before += dot(i->V, i->V);
      // Rotation
      VVectorr v = i->V + cross(i->V, T);
      i->V += cross(v, S);
      // accumulate new velocity
      sumv2After  += dot(i->V, i->V);
      // Velocity to Position
      for (int d=0; d<DIMR; ++d)
        i->R(d) += i->V(d) * dt;
    }
  }

  void Species::move3(Scheduler &scheduler)
  // Electrostatic + Magnetic forces
  {
    VVectorr E_extern(scheduler.getExternElectricField());
    VectorField E_intern;
    scheduler.getInternElectricField(E_intern);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // accumulate old velocity
      sumv2Before += dot(i->V, i->V);
      // Half-translation
      VVectorr Edt_2 = (E_extern +
                        interpolateField(E_intern, i->R, scheduler))*dtOver2;
      Edt_2 *= Z_sOverM_s;
      i->V += Edt_2;
      // Rotation
      VVectorr v = i->V + cross(i->V, T);
      i->V += cross(v, S);
      // Half-translation
      i->V += Edt_2;
      // accumulate new velocity
      sumv2After += dot(i->V, i->V);
      // Velocity to Position
      for (int d=0; d<DIMR; ++d)
        i->R(d) += i->V(d) * dt;
    }
  }

  void Species::move4(Scheduler &scheduler)
  // Gravitational force
  {
    VVectorr Gdt(scheduler.getExternGravitationField()*dt);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // accumulate old velocity
      sumv2Before += dot(i->V, i->V);
      // Translation
      i->V += Gdt;
      // accumulate new velocity
      sumv2After += dot(i->V, i->V);
      // Velocity to Position
      for (int d=0; d<DIMR; ++d)
        i->R(d) += i->V(d) * dt;
    }
  }

  void Species::move5(Scheduler &scheduler)
  // Electrostatic + Gravitational forces
  {
    VVectorr E_extern(scheduler.getExternElectricField());
    VectorField E_intern;
    scheduler.getInternElectricField(E_intern);
    VVectorr Gdt(scheduler.getExternGravitationField()*dt);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // accumulate old velocity
      sumv2Before += dot(i->V, i->V);
      // Translation
      VVectorr Edt = (E_extern +
                      interpolateField(E_intern, i->R, scheduler))*dt;
      Edt *= Z_sOverM_s;
      i->V += Edt;
      i->V += Gdt;
      // accumulate new velocity
      sumv2After += dot(i->V, i->V);
      // Velocity to Position
      for (int d=0; d<DIMR; ++d)
        i->R(d) += i->V(d) * dt;
    }
  }

  void Species::move6(Scheduler &scheduler)
  // Magnetic + Gravitational forces
  {
    VVectorr Gdt_2(scheduler.getExternGravitationField()*dtOver2);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // accumulate old velocity
      sumv2Before += dot(i->V, i->V);
      // Half-translation
      i->V += Gdt_2;
      // Rotation
      VVectorr v = i->V + cross(i->V, T);
      i->V += cross(v, S);
      // Half-translation
      i->V += Gdt_2;
      // accumulate new velocity
      sumv2After += dot(i->V, i->V);
      // Velocity to Position
      for (int d=0; d<DIMR; ++d)
        i->R(d) += i->V(d) * dt;
    }
  }

  void Species::move7(Scheduler &scheduler)
  // Electrostatic + Magnetic + Gravitational forces
  {
    VVectorr E_extern(scheduler.getExternElectricField());
    VectorField E_intern;
    scheduler.getInternElectricField(E_intern);
    VVectorr Gdt_2(scheduler.getExternGravitationField()*dtOver2);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // accumulate old velocity
      sumv2Before += dot(i->V, i->V);
      // Half-translation
      VVectorr Edt_2 = (E_extern +
                        interpolateField(E_intern, i->R, scheduler))*dtOver2;
      Edt_2 *= Z_sOverM_s;
      i->V += Edt_2;
      i->V += Gdt_2;
      // Rotation
      VVectorr v = i->V + cross(i->V, T);
      i->V += cross(v, S);
      // Half-translation
      i->V += Edt_2;
      i->V += Gdt_2;
      // accumulate new velocity
      sumv2After += dot(i->V, i->V);
      // Velocity to Position
      for (int d=0; d<DIMR; ++d)
        i->R(d) += i->V(d) * dt;
    }
  }

  void Species::moveParticles(NumberGenerator &draw, Scheduler &scheduler)
  {
    if (scheduler.getCurrentIter() == 0) {
      (this->*deccel[scheduler.getForceSetUp()])(scheduler);
    }
    // Kinetic energy calculation
    // from Plasma physics via computer simulation
    // Birsdall and Langdon, p.74 eq. 6
    sumv2Before   = 0.0;
    sumv2After    = 0.0;
    kineticEnergy = 0.0;

    (this->*move[scheduler.getForceSetUp()])(scheduler);

    fixBoundaryParticles(scheduler);
    injectParticles(draw, scheduler);

    kineticEnergy  = 0.5*m_s*0.5*(sumv2Before+sumv2After);
    kineticEnergy *= 2.0/3.0; // To convert in temperature
  }


  void Species::sync(int numLost, real v2Lost)
  {
    numParticles -= numLost;
    sumv2After   -= v2Lost;

    kineticEnergy  = 0.5*m_s*0.5*(sumv2Before+sumv2After);
    kineticEnergy *= 2.0/3.0; // To convert in temperature
  }


  void Species::deccel0(Scheduler &scheduler)
  // No forces
  {}


  void Species::deccel1(Scheduler &scheduler)
  // Electrostatic force
  {
    VVectorr E_extern(scheduler.getExternElectricField());
    VectorField E_intern;
    scheduler.getInternElectricField(E_intern);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // Half-translation
      VVectorr Edt_2 = (E_extern +
                        interpolateField(E_intern, i->R, scheduler))*dtOver2;
      Edt_2  *= Z_sOverM_s;
      i->V -= Edt_2;
    }
  }

  void Species::deccel2(Scheduler &scheduler)
  // Magnetic force
  {
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // Half rotation back
      real v1 = dot(i->V, i->V);
      i->V  -= cross(i->V, T);
      real v2 = dot(i->V, i->V);
      // normalise to keep same velocity!
      i->V  *= std::sqrt(v1/v2);
    }
  }

  void Species::deccel3(Scheduler &scheduler)
  // Electrostatic + Magnetic forces
  {
    VVectorr E_extern(scheduler.getExternElectricField());
    VectorField E_intern;
    scheduler.getInternElectricField(E_intern);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // Half-translation
      VVectorr Edt_2 = (E_extern +
                        interpolateField(E_intern, i->R, scheduler))*dtOver2;
      Edt_2  *= Z_sOverM_s;
      i->V -= Edt_2;
      // Half rotation back
      real v1 = dot(i->V, i->V);
      i->V  -= cross(i->V, T);
      real v2 = dot(i->V, i->V);
      // normalise to keep same velocity!
      i->V  *= std::sqrt(v1/v2);
    }
  }

  void Species::deccel4(Scheduler &scheduler)
  // Gravitational force
  {
    VVectorr Gdt_2(scheduler.getExternGravitationField()*dtOver2);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // Half-translation
      i->V -= Gdt_2;
    }
  }

  void Species::deccel5(Scheduler &scheduler)
  // Electrostatic + Gravitational forces
  {
    VVectorr E_extern(scheduler.getExternElectricField());
    VectorField E_intern;
    scheduler.getInternElectricField(E_intern);
    VVectorr Gdt_2(scheduler.getExternGravitationField()*dtOver2);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // Half-translation
      VVectorr Edt_2 = (E_extern +
                        interpolateField(E_intern, i->R, scheduler))*dtOver2;
      Edt_2  *= Z_sOverM_s;
      i->V -= Edt_2;
      i->V -= Gdt_2;
    }
  }

  void Species::deccel6(Scheduler &scheduler)
  // Magnetic + Gravitational forces
  {
    VVectorr Gdt_2(scheduler.getExternGravitationField()*dtOver2);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // Half-translation
      i->V -= Gdt_2;
      // Half rotation back
      real v1 = dot(i->V, i->V);
      i->V  -= cross(i->V, T);
      real v2 = dot(i->V, i->V);
      // normalise to keep same velocity!
      i->V  *= std::sqrt(v1/v2);
    }
  }

  void Species::deccel7(Scheduler &scheduler)
  // Electrostatic + Magnetic + Gravitational forces
  {
    VVectorr E_extern(scheduler.getExternElectricField());
    VectorField E_intern;
    scheduler.getInternElectricField(E_intern);
    VVectorr Gdt_2(scheduler.getExternGravitationField()*dtOver2);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      // Half-translation
      VVectorr Edt_2 = (E_extern +
                        interpolateField(E_intern, i->R, scheduler))*dtOver2;
      Edt_2 *= Z_sOverM_s;
      i->V -= Edt_2;
      i->V -= Gdt_2;
      // Half rotation back
      real v1 = dot(i->V, i->V);
      i->V  -= cross(i->V, T);
      real v2 = dot(i->V, i->V);
      // normalise to keep same velocity!
      i->V  *= std::sqrt(v1/v2);
    }
  }

  void Species::fixBoundaryParticles(Scheduler &scheduler)
  {
    Loss = 0;
    RVectorr delta_r(scheduler.getDomainSize());
    RVectorr domainMin(scheduler.getDomainMinBoundary());
    RVectorr domainMax(scheduler.getDomainMaxBoundary());
    ParticleListIter i = particles.begin(), e = particles.end();
    do {
      bool lost = false;
      for (int d=0; d<DIMR; ++d) {
        if (i->R(d) < domainMin(d)) {
          switch (BC.min(d)) {
          case periodic:
            i->R(d) += delta_r(d);
            // if exactly on the upper boundary
            // move the particle an `epsilon' inside the upper boundary
            if (i->R(d) == domainMax(d)) {
              i->R(d) -= i->R(d)*epsilon(static_cast<real>(1.0));
            }
            break;
          case insulated:
            i->R(d) = 2.0 * domainMin(d) - i->R(d);
            i->V(d) = -i->V(d);
            break;
          case free_space:
            if (!lost) { // do not count lost several times in several directions
              lost = true;
              Loss.min(d) += 1;
              numParticles--;
            }
            break;
          }
        }
        if (i->R(d) >= domainMax(d)) {
          switch (BC.max(d)) {
          case periodic:
            i->R(d) -= delta_r(d);
            break;
          case insulated:
            i->R(d) = 2.0 * domainMax(d) - i->R(d);
            i->V(d) = -i->V(d);
            // if exactly on the upper boundary
            // move an `epsilon' inside the upper boundary
            if (i->R(d) == domainMax(d)) {
              i->R(d) -= i->R(d)*epsilon(static_cast<real>(1.0));
            }
            break;
          case free_space:
            if (!lost) { // do not count lost several times in several directions
              lost = true;
              Loss.max(d) += 1;
              numParticles--;
            }
            break;
          }
        }
      }
      if (lost) {
        sumv2After -= dot(i->V, i->V);
        i = particles.erase(i);
      } else {
        ++i;
      }
    } while (i != e);

    BEGIN_DEBUG_OUTPUT(5)
    HEADER_DEBUG_OUTPUT2(speciesName,"fixBoundaryParticles")
    VAR_DEBUG_OUTPUT1(Loss)
    END_DEBUG_OUTPUT

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"fixBoundaryParticles")
    VAR_DEBUG_OUTPUT2(particles.size(),numParticles)
    END_DEBUG_OUTPUT
  }

  void Species::dump() const
  {
    string fname(speciesName);
#if defined(HAVE_MPI)

    ostringstream rank;
    rank << rankProc;
    fname += rank.str();
#endif

    fname += ".data";

    fstream file(fname.c_str(), ios::out | ios::binary);

    unsigned len = particles.size();
    file.write(reinterpret_cast<char*>(&len), sizeof(unsigned));

    ParticleListConstIter i = particles.begin(), e = particles.end();
    for ( ; i!=e; ++i) {
      file.write(const_cast<char*>
                 (reinterpret_cast<const char*>(i->R.data())),
                 DIMR * sizeof(real));
      file.write(const_cast<char*>
                 (reinterpret_cast<const char*>(i->V.data())),
                 DIMV * sizeof(real));
    }

  }

  void Species::restore()
  {}

  const std::string & Species::name() const
  {
    return speciesName;
  }

  int Species::getParticleNumber() const
  {
    return numParticles;
  }

  real Species::getKineticEnergy() const
  {
    return kineticEnergy;
  }

  Species::ParticleListIter Species::begin()
  {
    return particles.begin();
  }

  Species::ParticleListIter Species::end()
  {
    return particles.end();
  }

  Species::ParticleListConstIter Species::begin() const
  {
    return particles.begin();
  }

  Species::ParticleListConstIter Species::end() const
  {
    return particles.end();
  }

  unsigned Species::size() const
  {
    return particles.size();
  }

  Species::ParticleListIter Species::erase(ParticleListIter &i)
  {
    return particles.erase(i);
  }

  real Species::getChargeNumber() const
  {
    return Z_s;
  }

  real Species::getmu_s() const
  {
    return m_s;
  }

  const Boundary<real,DIMR> Species::getAlfa() const
  {
    return alfa;
  }

  void Species::initParsing(int nargs, char *args[])
  {
    registerClass(speciesName);
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    string prefix(speciesName);
    prefix += "::";
    setPrefix(prefix.c_str());

    parseLevelDebugOption("dl");

#if 0
    const char rveci[] =  "int[" DIMRSTR "]";
    const char rvecr[] = "real[" DIMRSTR "]";
    const char vvecr[] = "real[" DIMVSTR "]";
#endif

    using parser::types::boolean;
    using parser::types::real;
    using parser::types::intVect;
    using parser::types::realVect;

    insertOption(_mass   , "mass"       , real    , "Particles mass"        , Any(mass));
    insertOption(_charge , "charge"     , real    , "Particles charge"      , Any(charge));

    insertOption(_lambda , "lambda"     , real    , "Particles/Debye volume", Any(lambda));
    insertOption(_drift  , "drift"      , realVect, "Drift velocity"        , Any(drift));
    insertOption(_temp   , "temp"       , realVect, "Temperature"           , Any(temp));

    insertOption(_eqStart, "hydrostatic", boolean , "Initial cond"          , Any(eqStart));
    insertOption(_eqTemp , "Thydro"     , real    , "Hydrostatic Temp"      , Any(eqTemp));

    insertOption(_bcmin  , "bcmin"      , intVect , "LB condition"          , Any(BC.min));
    insertOption(_bcmax  , "bcmax"      , intVect , "UB condition"          , Any(BC.max));
    insertOption(_rmin   , "rmin"       , realVect, "LB coordinate"         , Any(domain.min));
    insertOption(_rmax   , "rmax"       , realVect, "UB coordinate"         , Any(domain.max));

#if defined(HAVE_MPI)

    insertOption(_cpuIds , "cpuid"      , intVect, "Vector of cpu's id"    , Any(cpuIds));
#endif

  }

  void Species::paramParsing()
  {
    parseOption(_mass   , mass);
    parseOption(_charge , charge);

    parseOption(_lambda , lambda);
    parseOption(_drift  , drift);
    parseOption(_temp   , temp);

    parseOption(_eqStart, eqStart);
    parseOption(_eqTemp , eqTemp);

    parseOption(_bcmin  , BC.min);
    parseOption(_bcmax  , BC.max);
    parseOption(_rmin   , domain.min);
    parseOption(_rmax   , domain.max);
  }

  void Species::checkParam() const
  {
    checkMap(_bcmin, bcMap, BC.min);
    checkMap(_bcmax, bcMap, BC.max);
  }

#undef ID

}
