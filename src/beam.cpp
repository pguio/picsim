/**************************************************************************
 *
 * $Id: beam.cpp,v 1.124 2016/06/02 17:19:03 patrick Exp $
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

#include <beam.h>
#include <species-factory.h>

namespace picsim {

#define ID "$Id: beam.cpp,v 1.124 2016/06/02 17:19:03 patrick Exp $"


  const int   DEFAULT_BEAMDIR     =  1;
#if (DIMR==2)

  const real DEFAULT_BEAM_LAMBDA =   5;
  const real DEFAULT_BEAM_DRIFT  =  -4;
  const real DEFAULT_BEAM_TEMP   = 0.5;
  const real DEFAULT_WIDTH       =  20;
  const real DEFAULT_DEPTH       =  40;
#elif (DIMR==3)

  const real DEFAULT_BEAM_LAMBDA =   2;
  const real DEFAULT_BEAM_DRIFT  =  -4;
  const real DEFAULT_BEAM_TEMP   = 0.5;
  const real DEFAULT_WIDTH       =   5;
  const real DEFAULT_DEPTH       =  10;
#endif

  const RVectorr DEFAULT_BEAM_POS = (mudfas::DEFAULT_GRIDSIZE-1)/2;

  const double M_1_SQRT2PI = 0.5*M_2_SQRTPI*M_SQRT1_2;

  namespace {
    RegisterInFactory<Species, Beam> registerMe("Beam");
  }

  namespace factory {
    void dummyBeam()
    {}
  }

  template <typename T_numtype>
  T_numtype leftInt(T_numtype x, T_numtype x0,
                    T_numtype m, T_numtype s, T_numtype dx)
  {
    using blitz::pow2;
    // \int 1/dx*( (x-x0)/dx+1)/sqrt(2*pi)/s.*exp(-(x-m).^2/2/s^2)
    return 1.0/pow2(dx)*(-s*M_1_SQRT2PI*exp(-0.5*pow2((x-m)/s))+
                         0.5*(m-x0+dx)*erf((x-m)*M_SQRT1_2/s));
  }

  template <typename T_numtype>
  T_numtype rightInt(T_numtype x, T_numtype x0,
                     T_numtype m, T_numtype s, T_numtype dx)
  {
    using blitz::pow2;
    // \int 1/dx*(-(x-x0)/dx+1)/sqrt(2*pi)/s.*exp(-(x-m).^2/2/s^2)
    return 1.0/pow2(dx)*( s*M_1_SQRT2PI*exp(-0.5*pow2((x-m)/s))-
                          0.5*(m-x0-dx)*erf((x-m)*M_SQRT1_2/s));
  }


  std::ostream& operator<<(std::ostream& os, const Beam &b)
  {
    b.printOn(os);
    return os;
  }


  Beam::Beam(int nargs, char *args[], const string name)
    : Species(nargs, args, name),
      beamDir(DEFAULT_BEAMDIR), beamPos(DEFAULT_BEAM_POS),
      beamWidth(DEFAULT_WIDTH), beamInitDepth(DEFAULT_DEPTH),
      timeSpec(0)
  {
    initParsing(nargs, args);
    paramParsing();
    checkParam();
  }

  Beam::~Beam()
  {}


  void Beam::initialise(NumberGenerator &draw, Scheduler &scheduler)
  {
    Species::initialise(draw, scheduler);

    // Calculate mean drift velocity of flux pdf
    // f(v) = v Gauss(v,u0,v0) /beta
    //      = v / (2 pi)^1/2 / v0 \exp(-(v-u0)^2/2/v0^2) / beta
    // meanDrift = \int_0^v v f(v) dv
    real v0 = sqrt(temp(beamDir)/mass);
    real u0 = drift(beamDir);
    real nu = std::abs(u0)/v0;
    real sgn = (u0 != 0.0 ? std::abs(u0)/u0 : 1.0);
    real beta = M_1_SQRT2PI*v0*exp(-0.5*nu*nu)+0.5*u0*(1.0+erf(nu/M_SQRT2));
    real meanDrift = sgn/beta*(0.5*(u0*u0+v0*v0)*(1.0+erf(nu/M_SQRT2))+
                               u0*v0*M_1_SQRT2PI*exp(-0.5*nu*nu));

    real flux = lambda * meanDrift;

    beamInitDepthTime = (beamDir <= DIMR-1 ? beamInitDepth(beamDir) / meanDrift : 0.0);

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
    VAR_DEBUG_OUTPUT3(nu,beta,meanDrift)
    VAR_DEBUG_OUTPUT2(lambda,flux)
    END_DEBUG_OUTPUT

    // Correction of the flux if gravitation and electric fields
    VVectorr G(scheduler.getExternGravitationField());
    if (scheduler.isGforce() && G(beamDir) != 0.0 &&
        scheduler.isEforce() && eqStart) {
      real K = m_s*std::abs(G(beamDir))/Teq;
      RVectorr dr = domain.max-domain.min;
      real nmax = dr(beamDir)*K/(1.0-exp(-dr(beamDir)*K));
      real nmin = nmax*exp(-dr(beamDir)*K);
      flux *= (flux*G(beamDir) > 0.0 ? nmin : nmax);

      BEGIN_DEBUG_OUTPUT(20)
      HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
      VAR_DEBUG_OUTPUT3(eqStart,Teq,K)
      VAR_DEBUG_OUTPUT2(nmax,nmin)
      END_DEBUG_OUTPUT
    }

    real area = getBeamArea(scheduler);

    BEGIN_DEBUG_OUTPUT(20)
    HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
    VAR_DEBUG_OUTPUT2(flux,area)
    VAR_DEBUG_OUTPUT2(drift,u_s)
    END_DEBUG_OUTPUT

    if (beamDir <= DIMR-1) {
      // Number of particles to reinject at each time step
      numToInject = (BC.min(beamDir) == periodic && BC.max(beamDir) == periodic ?
                     0 : static_cast<int>(std::abs(flux*area*dt)+0.5));
      // Initial number of particles
      numParticles = static_cast<int>(std::abs(lambda*area*beamInitDepth(beamDir)+0.5));
    } else {
      numToInject = 0;
      numParticles = static_cast<int>(std::abs(lambda*area+0.5));
    }

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
    VAR_DEBUG_OUTPUT2(numToInject,numParticles)
    END_DEBUG_OUTPUT

    // Check that time spec is an array of size at least the number of iterations
    if (! timeSpec.empty() > 0) {
      if (static_cast<int>(timeSpec.size()) < scheduler.getMaxIter())
        throw ClassException("Beam", "timespec size < maxiter");
    }

    // Allocate and initialise particles
    initialiseParticles(draw);
  }

  real Beam::getBeamArea(Scheduler &scheduler)
  {
    RVectori gridSize(scheduler.getGridSize());
    RVectorr rmin(scheduler.getDomainMinBoundary());
    RVectorr rmax(scheduler.getDomainMaxBoundary());
    RVectorr dr(scheduler.getDomainSize()/(gridSize-1));
    real Area = 1.0;
    for (int d=0; d<DIMR; ++d) {
      if (d != beamDir) {
        real dx = dr(d);
        real xl = rmin(d);
        real xh = xl+dx;
        real m = beamPos(d);
        real s = beamWidth(d);
        // integrate only the right hand part at lower boundary
        real prob = rightInt(xh, xl, m, s, dx)-rightInt(xl, xl, m, s, dx);
        real sum = prob;
        real max = prob;
        for (int i=1; i<gridSize(d)-1; ++i) { // for all grid points except the boundary
          prob = leftInt(xh, xh, m, s, dx)-leftInt(xl, xh, m, s, dx);
          xl = xh;
          xh = xl+dx;
          prob += rightInt(xh, xl, m, s, dx)-rightInt(xl, xl, m, s, dx);
          sum += prob;
          max = std::max(max, prob);
        }
        // integrate only the left hand part at upper boundary
        prob = leftInt(xh, xh, m, s, dx)-leftInt(xl, xh, m, s, dx);
        sum += prob;
        max = std::max(max, prob);

        Area *= sum/max*dr(d);
      }
    }
    return Area;
  }

  void Beam::initialiseParticles(NumberGenerator &draw)
  {
    particles.resize(numParticles);
    real T = beamInitDepthTime;

    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      for (int d=0; d<DIMR; ++d) {
        real a = domain.min(d);
        real b = domain.max(d);
        if (d == beamDir) {
          if (BC.min(beamDir) == periodic && BC.max(beamDir) == periodic) {
            // Periodic condition
            // beam is spatially uniformely distributed along the beam
            // direction
            i->R(beamDir) = draw.uniform(a, b);
            i->V(beamDir) = draw.flux(u_s(d), v_s(d));
          } else {
            // Otherwise
            // Estimate first the velocity and calculate the distance the
            // particle can have moved in the simulation area
            if (u_s(d) >= 0.0) {
              i->V(d) = draw.flux(u_s(d), v_s(d));
              i->R(d) = draw.uniform(a, std::min(a+i->V(d)*T, b));
            } else {
              i->V(d) = -draw.flux(-u_s(d), v_s(d));
              i->R(d) = draw.uniform(std::max(a, b+i->V(d)*T), b);
            }
          }
        } else {
          // Direction perpendicular to the beam
          // normal distribution of the beam in configuration space
          do {
            i->R(d) = draw.normal(beamPos(d), beamWidth(d));
          } while (i->R(d) < a || i->R(d) >= b);
          i->V(d) = draw.normal(u_s(d), v_s(d));
        }
      }
#if (DIMR==2)
      i->V(DIMV-1) = draw.normal(u_s(DIMV-1), v_s(DIMV-1));
#endif

    }
  }

  void Beam::injectParticles(NumberGenerator &draw, Scheduler &scheduler)
  {
    unsigned n = (! timeSpec.empty() ? (timeSpec[scheduler.getCurrentIter()]*numToInject)/1000 :
                  numToInject);

    // If beam direction outside simulation plane
    // inject possibly lost particles with free_space boundary only
    // Otherwise should be zero
    if (beamDir > DIMR-1)
      n = blitz::sum(Loss.min)+blitz::sum(Loss.max);

    ParticleList p(n);

    for (ParticleListIter i=p.begin(), e=p.end(); i!=e; ++i) {
      for (int d=0; d<DIMR; ++d) {
        real a = domain.min(d);
        real b = domain.max(d);
        if (d == beamDir) {
          if (u_s(d) >= 0.0) {
            i->V(d) = draw.flux(u_s(d), v_s(d));
            i->R(d) = draw.uniform(a, std::min(a+i->V(d)*dt, b));
          } else {
            i->V(d) = -draw.flux(-u_s(d), v_s(d));
            i->R(d) = draw.uniform(std::max(a, b+i->V(d)*dt), b);
          }
        } else {
          do {
            i->R(d) = draw.normal(beamPos(d), beamWidth(d));
          } while (i->R(d) < a || i->R(d) >= b);
          i->V(d) = draw.normal(u_s(d), v_s(d));
        }
      }
#if (DIMR==2)
      i->V(DIMV-1) = draw.normal(u_s(DIMV-1), v_s(DIMV-1));
#endif

      sumv2After += blitz::dot(i->V, i->V);
    }

    BEGIN_DEBUG_OUTPUT(5)
    HEADER_DEBUG_OUTPUT2(speciesName,"injectParticles")
    VAR_DEBUG_OUTPUT2(p.size(),numToInject)
    END_DEBUG_OUTPUT

    particles.insert(particles.end(), p.begin(), p.end());
    numParticles += n;

    BEGIN_DEBUG_OUTPUT(5)
    HEADER_DEBUG_OUTPUT2(speciesName,"injectParticles")
    VAR_DEBUG_OUTPUT2(particles.size(),numParticles)
    END_DEBUG_OUTPUT
  }

  void Beam::printOn(std::ostream& os) const
  {

    Species::printOn(os);

    const char* dir[]= {"x", "y", "z"
                       };
    using parser::yesno;

    os
        << "Beam direction          = " << dir[beamDir] << '\n'
        << "Beam position           = " << beamPos << '\n'
        << "Beam Initial depth      = " << (beamDir <= DIMR-1 ? beamInitDepth(beamDir) : 0.0) << '\n'
        << "Beam width              = " << beamWidth << '\n'
        << "Beam time spec          = " << yesno(!timeSpec.empty()) << '\n';
  }

  void Beam::initParsing(int nargs, char *args[])
  {
#if 0
    const char rvecr[] = "real[" DIMRSTR "]";
#endif

    registerClass(speciesName);
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    using parser::types::integer;
    using parser::types::realVect;
    using parser::types::intVect;

    insertOption(_beamDir      , "dir"     , integer , "Beam direction"                       , Any(beamDir));
    insertOption(_beamPos      , "pos"     , realVect, "Beam position"                        , Any(beamPos));

    insertOption(_beamWidth    , "width"   , realVect, "Beam width   "                        , Any(beamWidth));
    insertOption(_beamInitDepth, "depth"   , realVect, "Beam depth @ t=0"                     , Any(beamInitDepth));

    insertOption(_timeSpec     , "timespec", intVect, "Time spec for flux intensity [0,1000]", Any(timeSpec));
  }

  void Beam::paramParsing()
  {
    parseOption(_beamDir      , beamDir);
    parseOption(_beamPos      , beamPos);

    parseOption(_beamWidth    , beamWidth);
    parseOption(_beamInitDepth, beamInitDepth);

    parseOption(_timeSpec     , timeSpec);
  }

  void Beam::checkParam() const
  {
    if (beamDir < 0 || beamDir > DIMV)
      throw ClassException("Beam", "Beam direction is  0, 1 or 2");
#if 0

    if (drift(beamDir) == 0.0 && beamInitDepth(beamDir) == 0.0)
      throw ClassException("Beam", "Beam drift velocity and depth are both 0.0!");
#endif

  }

#undef ID

}
