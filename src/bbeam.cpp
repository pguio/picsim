/**************************************************************************
 *
 * $Id: bbeam.cpp,v 1.7 2011/03/26 15:36:08 patrick Exp $
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

#include <bbeam.h>
#include <species-factory.h>

namespace picsim {

#define ID "$Id: bbeam.cpp,v 1.7 2011/03/26 15:36:08 patrick Exp $"

  namespace {
    RegisterInFactory<Species, BBeam> registerMe("BBeam");
  }

  namespace factory {
    void dummyBBeam()
    {}
  }

  std::ostream& operator<<(std::ostream& os, const BBeam &b)
  {
    b.printOn(os);
    return os;
  }


  BBeam::BBeam(int nargs, char *args[], const string name) : Beam(nargs, args, name)
  {
    initParsing(nargs, args);
    paramParsing();
    checkParam();
  }

  BBeam::~BBeam()
  {}


  void BBeam::initialise(NumberGenerator &draw, Scheduler &scheduler)
  {
    Species::initialise(draw, scheduler);

    // Gamma  = n u0
    real flux = lambda * drift(beamDir);

    // Correction of the flux if gravitation and electric fields
    VVectorr G(scheduler.getExternGravitationField());
    if (scheduler.isGforce() && G(beamDir) != 0.0 && scheduler.isEforce()) {
      real K = (eqStart ? m_s*std::abs(G(beamDir))/Teq : 0.0);
      RVectorr dr = domain.max-domain.min;
      real nmax = dr(beamDir)*K/(1.0-exp(-dr(beamDir)*K));
      real nmin = nmax*exp(-dr(beamDir)*K);
      flux *= (flux*G(beamDir) > 0.0 ? nmin : nmax);
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
        throw ClassException("BBeam", "timespec size < maxiter");
    }

    // Allocate and initialise particles
    initialiseParticles(draw);
  }


  void BBeam::initialiseParticles(NumberGenerator &draw)
  {
    particles.resize(numParticles);

    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      for (int d=0; d<DIMR; ++d) {
        real a = domain.min(d);
        real b = domain.max(d);
        if (d == beamDir) {
          if (BC.min(beamDir) == periodic && BC.max(beamDir) == periodic) {
            // Periodic condition
            // beam is spatially uniformely distributed along the beam direction
            i->R(beamDir) = draw.uniform(a, b);
            i->V(beamDir) = draw.normal(u_s(d), v_s(d));
          } else {
            // Otherwise uniform over beamDepthInit
            if (u_s(d) >= 0.0) {
              i->R(d) = draw.uniform(a, std::min(a+beamInitDepth(beamDir), b));
              i->V(d) = draw.normal(u_s(d), v_s(d));
            } else {
              i->R(d) = draw.uniform(std::max(a, b-beamInitDepth(beamDir)), b);
              i->V(d) = draw.flux(u_s(d), v_s(d));
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

  void BBeam::injectParticles(NumberGenerator &draw, Scheduler &scheduler)
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

  void BBeam::printOn(std::ostream& os) const
  {
    Beam::printOn(os);
  }

  void BBeam::initParsing(int nargs, char *args[])
  {
    registerClass(speciesName);
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);
  }

  void BBeam::paramParsing()
  {}

  void BBeam::checkParam() const
  {}

#undef ID

}
