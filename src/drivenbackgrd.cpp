/**************************************************************************
 *
 * $Id: drivenbackgrd.cpp,v 1.81 2017/11/26 18:22:11 patrick Exp $
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

#include <drivenbackgrd.h>
#include <species-factory.h>

namespace picsim {

  using blitz::dot;
  using blitz::pow2;
  using blitz::rint;
  using blitz::sum;

  using parser::map_elt;

  using std::ostream;

#define ID "$Id: drivenbackgrd.cpp,v 1.81 2017/11/26 18:22:11 patrick Exp $"

  const double M_1_SQRT2PI(0.5*M_2_SQRTPI/M_SQRT2);

  namespace {
    RegisterInFactory<Species, DrivenBackgrd> registerMe("DrivenBackground");
  }

  namespace factory {
    void dummyDrivenBgrd()
    {}
  }

  ostream& operator<<(ostream& os, const DrivenBackgrd &b)
  {
    b.printOn(os);
    return os;
  }


  DrivenBackgrd::DrivenBackgrd(int nargs, char *args[], const string name)
    : Backgrd(nargs, args, name),
      Frqcy(0.0,0.0), A(0.0,0.0), Pos(0.5,0.5),
      Start(0.0,0.0), Stop(0.0,0.0)
  {
    perturbMap.insert(Pair(nope    , "nothing"));
    perturbMap.insert(Pair(harmonic, "harmonic"));
    perturbMap.insert(Pair(pulse   , "pulse"));

    initParsing(nargs, args);
    paramParsing();
    checkParam();
  }

  DrivenBackgrd::~DrivenBackgrd()
  {}

  void DrivenBackgrd::initialise(NumberGenerator &draw, Scheduler &scheduler)
	{
	  Backgrd::initialise(draw, scheduler);

    //injectParticles(draw, scheduler);
	}

  void DrivenBackgrd::injectParticles(NumberGenerator &draw, Scheduler &scheduler)
  {

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"injectParticles")
    VAR_DEBUG_OUTPUT1(Loss)
    END_DEBUG_OUTPUT

    // Trick
    // f bulk velocity != 0 and
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
        if (u_s(d) != 0.0)
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
        if (u_s(d) != 0.0)
          Loss.max(d) = loss_min;
        break;
      case noinject:
        Loss.max(d) = 0;
        break;
      }
    }

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"injectParticles")
    VAR_DEBUG_OUTPUT1(Loss)
    END_DEBUG_OUTPUT

    real t(scheduler.getCurrentTime());
    Boundary<real,DIMR> alpha(0.0,0.0);
    Boundary<real,DIMR> omega(2.0*M_PI*Frqcy.min,2.0*M_PI*Frqcy.max);
    for (int d=0; d<DIMR; ++d) {
      if (Perturbation.min(d) == harmonic) {
        // Harmonic perturbation
        // mean value of sin(\omega t) in [t,t+dt]
        // i.e. 1/t\int_t^{t+dt}sin(\omega t) dt
        if (A.min(d) != 0.0 && omega.min(d) != 0.0) {
          alpha.min(d) = -A.min(d)/(omega.min(d)*dt)*
                         (cos(omega.min(d)*(t+dt))-cos(omega.min(d)*t));
        }
      } else if (Perturbation.min(d) == pulse &&
                 t >= Start.min(d) && t <= Stop.min(d)) {
        // Pulsed perturbation only at first iteration
        alpha.min(d) = A.min(d);
      }
      if (Perturbation.max(d) == harmonic) {
        // Harmonic perturbation
        if (A.max(d) != 0.0 && omega.max(d) != 0.0) {
          alpha.max(d) = -A.max(d)/(omega.max(d)*dt)*
                         (cos(omega.max(d)*(t+dt))-cos(omega.max(d)*t));
        }
      } else if (Perturbation.max(d) == pulse &&
                 t >= Start.max(d) && t <= Stop.max(d)) {
        // Pulsed perturbation only at first iteration
        alpha.max(d) = A.max(d);
      }
    }
    Boundary<int,DIMR> Inject(rint((1.0+alpha.min)*Loss.min),
                              rint((1.0+alpha.max)*Loss.max));

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"injectParticles")
    VAR_DEBUG_OUTPUT2(alpha,Inject)
    END_DEBUG_OUTPUT

    numToInject = sum(Inject.min+Inject.max);

    ParticleList p(numToInject);
    ParticleListIter i(p.begin());
    for (int d=0; d<DIMR; ++d) { // For all spatial dimensions
      for (int n=0; n<Inject.min(d); ++n) { // Lost at minimum boundary
        real a1(1.0); // estimate alpha in the first direction for specified r
        for (int pd=0; pd<DIMR; ++pd) {
          real a = domain.min(pd);
          real b = domain.max(pd);
          if (pd == d) { // perpendicular direction to the face
            i->V(pd) = draw.flux(u_s(pd), v_s(pd));
            i->R(pd) = draw.uniform(a, std::min(a+i->V(pd)*dt, b));
          } else { // parallel directions to the face
            if (alpha.min(d) == 0.0) { // No perturbation
              i->R(pd) = draw.uniform(a, b);
            } else { // Perturbation
              real mu(Pos.min(pd));
              real s(Dev.min(pd));
              i->R(pd) = draw.perturb(a, b, mu, s, a1*alpha.min(d));
              real x((i->R(pd)-a)/(b-a));
              a1 = static_cast<real>(M_1_SQRT2PI/s
                                     *exp(-static_cast<double>(pow2((x-mu)/s)/2.0))
                                     *2.0/(erf((1.0-mu)/M_SQRT2/s)+erf(mu/M_SQRT2/s)));
            }
            i->V(pd) = draw.normal(u_s(pd), v_s(pd));
          }
        }
#if (DIMR==2)
        i->V(DIMV-1) = draw.normal(u_s(DIMV-1), v_s(DIMV-1));
#endif

        sumv2After += dot(i->V, i->V);
        ++i;
      }
      for (int n=0; n<Inject.max(d); ++n) {
        real a1(1.0);
        for (int pd=0; pd<DIMR; ++pd) {
          real a = domain.min(pd);
          real b = domain.max(pd);
          if (pd == d) { // perpendicular direction to the face
            i->V(pd) = -draw.flux(u_s(pd), v_s(pd));
            i->R(pd) = draw.uniform(std::max(a, b+i->V(pd)*dt), b);
          } else { // parallel  directions to the face
            if (alpha.max(d) == 0.0) { // No perturbation
              i->R(pd) = draw.uniform(a, b);
            } else { // Perturbation
              real mu(Pos.max(pd));
              real s(Dev.max(pd));
              i->R(pd) = draw.perturb(a, b, mu, s, a1*alpha.max(d));
              real x((i->R(pd)-a)/(b-a));
              a1 = static_cast<real>(M_1_SQRT2PI/s
                                     *exp(-static_cast<double>(pow2((x-mu)/s)/2.0))
                                     *2.0/(erf((1.0-mu)/M_SQRT2/s)+erf(mu/M_SQRT2/s)));
            }
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

  void DrivenBackgrd::printOn(ostream& os) const
  {
    Backgrd::printOn(os);

    os
        << "Perturbation frequency  = " << Frqcy << '\n'
        << "Perturbation amplitude  = " << A << '\n'
        << "Perturbation position   = " << Pos << '\n'
        << "Perturbation std        = " << Dev << '\n'
        << "Perturbation type min   = "
        << map_elt(perturbMap, Perturbation.min) << '\n'
        << "Perturbation type max   = "
        << map_elt(perturbMap, Perturbation.max) << '\n'
        << "Pulse start time        = " << Start << '\n'
        << "Pulse stop  time        = " << Stop;
  }

  void DrivenBackgrd::initParsing(int nargs, char *args[])
  {
    registerClass(speciesName);
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

#if 0
    const char rvecr[] = "real[" DIMRSTR "]";
#endif

    using parser::types::realVect;

    insertOption(fmin      , "fmin"      , realVect, "LB perturbation frequency", Any(Frqcy.min));
    insertOption(fmax      , "fmax"      , realVect, "UB perturbation frequency", Any(Frqcy.max));
    insertOption(amin      , "amin"      , realVect, "LB perturbation amplitude", Any(A.min));
    insertOption(amax      , "amax"      , realVect, "UB perturbation amplitude", Any(A.max));

    insertOption(posmin    , "posmin"    , realVect, "LB perturbation position" , Any(Pos.min));
    insertOption(posmax    , "posmax"    , realVect, "UB perturbation position" , Any(Pos.max));
    insertOption(devmin    , "devmin"    , realVect, "LB perturbation deviation", Any(Dev.min));
    insertOption(devmax    , "devmax"    , realVect, "UB perturbation deviation", Any(Dev.max));

    insertOption(perturbmin, "perturbmin", realVect, "LB perturbation type"     , Any(Perturbation.min));
    insertOption(perturbmax, "perturbmax", realVect, "UB perturbation type"     , Any(Perturbation.max));

    insertOption(startmin  , "startmin"  , realVect, "LB pulse start"           , Any(Start.min));
    insertOption(startmax  , "startmax"  , realVect, "UB pulse start"           , Any(Start.max));
    insertOption(stopmin   , "stopmin"   , realVect, "LB pulse stop"            , Any(Stop.min));
    insertOption(stopmax   , "stopmax"   , realVect, "UB pulse stop"            , Any(Stop.max));
  }

  void DrivenBackgrd::paramParsing()
  {
    parseOption(fmin      , Frqcy.min);
    parseOption(fmax      , Frqcy.max);
    parseOption(amin      , A.min);
    parseOption(amax      , A.max);
    parseOption(posmin    , Pos.min);
    parseOption(posmax    , Pos.max);
    parseOption(devmin    , Dev.min);
    parseOption(devmax    , Dev.max);
    parseOption(perturbmin, Perturbation.min);
    parseOption(perturbmax, Perturbation.max);
    parseOption(startmin  , Start.min);
    parseOption(startmax  , Start.max);
    parseOption(stopmin   , Stop.min);
    parseOption(stopmax   , Stop.max);
  }

  void DrivenBackgrd::checkParam() const
  {
    checkMap(perturbmin, perturbMap, Perturbation.min);
    checkMap(perturbmax, perturbMap, Perturbation.max);

    // Try out parameter alpha for perturbation
    // Check whether it is not larger than maximum allowed
    NumberGenerator draw;
    for (int d=0; d<DIMR; ++d) {
      if (Perturbation.min(d) == harmonic && A.min(d) != 0.0) {
        real a1(1.0);
        for (int pd=0; pd<DIMR; ++pd) {
          if (pd != d) {

            BEGIN_DEBUG_OUTPUT(20)
            HEADER_DEBUG_OUTPUT2(speciesName,"checkParam")
            VAR_DEBUG_OUTPUT1(a1)
            END_DEBUG_OUTPUT

            draw.perturb(domain.min(pd), domain.max(pd),
                         Pos.min(pd), Dev.min(pd), -a1*A.min(d));
            real mu = Pos.min(pd);
            real s  = Dev.min(pd);
            // In the case DIMR==3, p(x,y) = p(y|x) p(x)
            // maximum a1=p(x) when x=\mu
            a1 = (real)(M_1_SQRT2PI/s*2.0/
                        (erf((1.0-mu)/M_SQRT2/s)+erf(mu/M_SQRT2/s)));
          }
        }
      }
      if (Perturbation.max(d) == harmonic && A.max(d) != 0.0) {
        real a1(1.0);
        for (int pd=0; pd<DIMR; ++pd) {
          if (pd != d) {

            BEGIN_DEBUG_OUTPUT(20)
            HEADER_DEBUG_OUTPUT2(speciesName,"checkParam")
            VAR_DEBUG_OUTPUT1(a1)
            END_DEBUG_OUTPUT

            draw.perturb(domain.min(pd), domain.max(pd),
                         Pos.max(pd), Dev.max(pd), -a1*A.max(d));
            real mu(Pos.max(pd));
            real s(Dev.max(pd));
            a1 = (real)(M_1_SQRT2PI/s*2.0/
                        (erf((1.0-mu)/M_SQRT2/s)+erf(mu/M_SQRT2/s)));
          }
        }
      }
    }
  }

#undef ID

}
