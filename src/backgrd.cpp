/**************************************************************************
 *
 * $Id: backgrd.cpp,v 1.106 2011/03/26 15:36:07 patrick Exp $
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

#include <backgrd.h>
#include <species-factory.h>

namespace picsim {

  using blitz::dot;
  using blitz::product;
  using blitz::sum;

  using parser::map_elt;

  using std::ostream;

#define ID "$Id: backgrd.cpp,v 1.106 2011/03/26 15:36:07 patrick Exp $"

  const double M_1_SQRT2PI = 0.5*M_2_SQRTPI*M_SQRT1_2;

  namespace {
    RegisterInFactory<Species, Backgrd> registerMe("Background");
  }

  namespace factory {
    void dummyBgrd()
    {}
  }


  ostream& operator<<(ostream& os, const Backgrd &b)
  {
    b.printOn(os);
    return os;
  }


  Backgrd::Backgrd(int nargs, char *args[], const string name)
    : Species(nargs, args, name),
      injectMode(conserving,conserving), injectCoef(0.0, 0.0),
      numInjectStat(0,0)
  {
    imodeMap.insert(Pair(conserving,"PN Conserved"));
    imodeMap.insert(Pair(statistic ,"   Statistic"));
    imodeMap.insert(Pair(noinject  ,"No Injection"));

    initParsing(nargs, args);
    paramParsing();
    checkParam();
  }

  Backgrd::~Backgrd()
  {}

  void Backgrd::initialise(NumberGenerator &draw, Scheduler &scheduler)
  {
    Species::initialise(draw, scheduler);

    RVectorr deltar = domain.max-domain.min;
    real V = product(deltar);
    VVectorr nu(u_s/v_s);
    RVectorr Vmin, Vmax;

    for (int d=0; d<DIMR; ++d) {
      Vmin(d)=M_1_SQRT2PI*v_s(d)*exp(double(-0.5*nu(d)*nu(d)))+
              0.5*u_s(d)*(1.0+erf(double(nu(d)/M_SQRT2)));
      Vmax(d)=M_1_SQRT2PI*v_s(d)*exp(double(-0.5*nu(d)*nu(d)))-
              0.5*u_s(d)*(1.0-erf(double(nu(d)/M_SQRT2)));
    }

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
    VAR_DEBUG_OUTPUT2(Vmin,Vmax)
    END_DEBUG_OUTPUT

    numInjectStat = Boundary<int,DIMR>(lambda*V/deltar*Vmin*dt,
                                       lambda*V/deltar*Vmax*dt);
    BEGIN_DEBUG_OUTPUT(5)
    HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
    VAR_DEBUG_OUTPUT1(numInjectStat)
    END_DEBUG_OUTPUT

    if (scheduler.isGforce() && scheduler.isEforce()) {
      VVectorr G(scheduler.getExternGravitationField());
      RVectori gridsize(scheduler.getGridSize());
      RVectorr dr(deltar/(gridsize-1));
      for (int d=0; d<DIMR; ++d) {
        if (G(d) != 0.0) {
          alfa.max(d) = alfa.min(d) = m_s*G(d)/Teq;
          if (G(d) > 0.0) {
            real nmax = deltar(d)*alfa.max(d)/(1.0-exp(-deltar(d)*alfa.max(d)));
            real nmin = nmax*exp(-deltar(d)*alfa.max(d));

            BEGIN_DEBUG_OUTPUT(20)
            HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
            VAR_DEBUG_OUTPUT2(nmax,nmin)
            END_DEBUG_OUTPUT

            numInjectStat.min(d) = static_cast<int>(lambda*nmin*V/
                                                    deltar(d)*Vmin(d)*dt+0.5);
            numInjectStat.max(d) = static_cast<int>(lambda*nmax*V/
                                                    deltar(d)*Vmax(d)*dt+0.5);
          } else {
            real nmax = -deltar(d)*alfa.max(d)/(1.0-exp(deltar(d)*alfa.max(d)));
            real nmin = nmax*exp(deltar(d)*alfa.max(d));

            BEGIN_DEBUG_OUTPUT(20)
            HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
            VAR_DEBUG_OUTPUT2(nmax,nmin)
            END_DEBUG_OUTPUT

            numInjectStat.min(d) = static_cast<int>(lambda*nmax*V/
                                                    deltar(d)*Vmin(d)*dt+0.5);
            numInjectStat.max(d) = static_cast<int>(lambda*nmin*V/
                                                    deltar(d)*Vmax(d)*dt+0.5);
          }
        }
      }
      BEGIN_DEBUG_OUTPUT(5)
      HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
      VAR_DEBUG_OUTPUT2(numInjectStat,alfa)
      END_DEBUG_OUTPUT
    }

    numParticles = static_cast<int>(lambda*product(deltar)+0.5);
    particles.resize(numParticles);

    BEGIN_DEBUG_OUTPUT(20)
    HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
    VAR_DEBUG_OUTPUT2(numParticles,particles.size())
    END_DEBUG_OUTPUT

    if (scheduler.isGforce() && scheduler.isEforce() && eqStart)
      initialiseExponential(draw);
    else
      initialiseUniform(draw);
  }

  void Backgrd::printOn(ostream& os) const
  {
    Species::printOn(os);

    os
        << "Inject mode min         = " << map_elt(imodeMap, injectMode.min) << '\n'
        << "Inject mode max         = " << map_elt(imodeMap, injectMode.max) << '\n'
        << "Inject coef min         = " << injectCoef.min << '\n'
        << "Inject coef max         = " << injectCoef.max << '\n';
  }

  void Backgrd::initialiseUniform(NumberGenerator &draw)
  {
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      for (int d=0; d<DIMR; ++d) {
        i->R(d) = draw.uniform(domain.min(d), domain.max(d));
      }
      for (int d=0; d<DIMV; ++d) {
        i->V(d) = draw.normal(u_s(d), v_s(d));
      }
    }
  }

  void Backgrd::initialiseExponential(NumberGenerator &draw)
  {
    RVectorr s(alfa.max);
    for (ParticleListIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      for (int d=0; d<DIMR; ++d) {
        real a = domain.min(d);
        real b = domain.max(d);
        if (s(d) > 0.0) {
          do {
            i->R(d) = b - draw.exponential(s(d));
          } while (i->R(d) < a);
        } else if (s(d) < 0.0) {
          do {
            i->R(d) = a + draw.exponential(-s(d));
          } while (i->R(d) >= b);
        } else {
          i->R(d) = draw.uniform(a, b);
        }
      }
      for (int d=0; d<DIMV; ++d) {
        i->V(d) = draw.normal(u_s(d), v_s(d));
      }
    }
  }


  void Backgrd::injectParticles(NumberGenerator &draw, Scheduler &scheduler)
  {
    // Trick
    // if bulk velocity != 0 and
    // number of particles is conserved then
    // reinject what is lost at one side into the opposite side
    for (int d=0; d<DIMR; ++d) {
      int loss_min(Loss.min(d));
      int loss_max(Loss.max(d));
      switch (injectMode.min(d)) {
      case statistic:
        Loss.min(d) = numInjectStat.min(d);
        break;
      case conserving:
        if (u_s(d)!=0.0)
          Loss.min(d) = loss_max;
        break;
      case noinject:
        Loss.min(d) = 0;
        break;
      }
      switch (injectMode.max(d)) {
      case statistic:
        Loss.max(d) = numInjectStat.max(d);
        break;
      case conserving:
        if (u_s(d)!=0.0)
          Loss.max(d) = loss_min;
        break;
      case noinject:
        Loss.max(d) = 0;
        break;
      }
    }

    for (int d=0; d<DIMR; ++d) {
      Loss.min(d) = static_cast<int>((1.0+injectCoef.min(d))*Loss.min(d)+0.5);
      Loss.max(d) = static_cast<int>((1.0+injectCoef.max(d))*Loss.max(d)+0.5);
    }

    numToInject = sum(Loss.min+Loss.max);

    ParticleList p(numToInject);
    ParticleListIter i = p.begin();
    for (int d=0; d<DIMR; ++d) { // For all spatial dimensions
      for (int n=0; n<Loss.min(d); ++n) { // Lost at minimum boundary
        for (int pd=0; pd<DIMR; ++pd) {
          real a = domain.min(pd);
          real b = domain.max(pd);
          if (pd == d) { // perpendicular direction to the face
            i->V(pd) = draw.flux(u_s(pd), v_s(pd));
            i->R(pd) = draw.uniform(a, std::min(a+i->V(pd)*dt, b));
          } else { // parallel directions to the face
            i->R(pd) = draw.uniform(a, b);
            i->V(pd) = draw.normal(u_s(pd), v_s(pd));
          }
        }
#if (DIMR==2)
        i->V(DIMV-1) = draw.normal(u_s(DIMV-1), v_s(DIMV-1));
#endif

        sumv2After += dot(i->V, i->V);
        ++i;
      }
      for (int n=0; n<Loss.max(d); ++n) {
        for (int pd=0; pd<DIMR; ++pd) {
          real a = domain.min(pd);
          real b = domain.max(pd);
          if (pd == d) { // perpendicular direction to the face
            i->V(pd) = -draw.flux(-u_s(pd), v_s(pd));
            i->R(pd) = draw.uniform(std::max(a, b+i->V(pd)*dt), b);
          } else { // parallel directions
            i->R(pd) = draw.uniform(a, b);
            i->V(pd) = draw.normal(u_s(pd), v_s(pd));
          }
        }
#if (DIMR==2)
        i->V(DIMV-1) = draw.normal(u_s(DIMV-1), v_s(DIMV-1));
#endif

        sumv2After += dot(i->V, i->V);
        ++i;
      }
    }

    BEGIN_DEBUG_OUTPUT(5)
    HEADER_DEBUG_OUTPUT2(speciesName,"injectParticles")
    VAR_DEBUG_OUTPUT2(p.size(),numToInject)
    END_DEBUG_OUTPUT

    particles.insert(particles.end(), p.begin(), p.end());
    numParticles += numToInject;

    BEGIN_DEBUG_OUTPUT(5)
    HEADER_DEBUG_OUTPUT2(speciesName,"injectParticles")
    VAR_DEBUG_OUTPUT2(particles.size(),numParticles)
    END_DEBUG_OUTPUT
  }

  void Backgrd::initParsing(int nargs, char *args[])
  {
    registerClass(speciesName);
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

#if 0
    const char rveci[] =  "int[" DIMRSTR "]";
    const char rvecr[] = "real[" DIMRSTR "]";
#endif

    using parser::types::intVect;
    using parser::types::realVect;

    insertOption(_imodemin, "imodemin", intVect ,"LB inject mode", Any(injectMode.min));
    insertOption(_imodemax, "imodemax", intVect ,"UB inject mode", Any(injectMode.max));
    insertOption(_icoefmin, "icoefmin", realVect,"LB inject coef", Any(injectCoef.min));
    insertOption(_icoefmax, "icoefmax", realVect,"UB inject coef", Any(injectCoef.max));
  }

  void Backgrd::paramParsing()
  {
    parseOption(_imodemin, injectMode.min);
    parseOption(_imodemax, injectMode.max);
    parseOption(_icoefmin, injectCoef.min);
    parseOption(_icoefmax, injectCoef.max);
  }

  void Backgrd::checkParam() const
  {
    checkMap(_imodemin, imodeMap, injectMode.min);
    checkMap(_imodemax, imodeMap, injectMode.max);
  }

#undef ID

}
