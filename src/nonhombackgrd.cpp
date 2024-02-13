/**************************************************************************
 *
 * $Id: nonhombackgrd.cpp,v 1.16 2015/11/06 18:24:54 patrick Exp $
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

#include <fstream>
#include <nonhombackgrd.h>
#include <species-factory.h>

namespace picsim {

  using blitz::cast;
  using blitz::dot;
  using blitz::epsilon;
  using blitz::product;
  using blitz::sum;

  using blitz::firstDim;
  using blitz::secondDim;
#if (DIMR==3)
  using blitz::thirdDim;
#endif

  using parser::map_elt;

  using std::cout;
  using std::endl;
  using std::ostream;

#define ID "$Id: nonhombackgrd.cpp,v 1.16 2015/11/06 18:24:54 patrick Exp $"

  const double M_1_SQRT2PI = 0.5*M_2_SQRTPI*M_SQRT1_2;

  namespace {
    RegisterInFactory<Species, NonHomBackgrd> registerMe("NonHomBackground");
  }

  namespace factory {
    void dummyNonHomBgrd()
    {}
  }


  ostream& operator<<(ostream& os, const NonHomBackgrd &b)
  {
    b.printOn(os);
    return os;
  }


  NonHomBackgrd::NonHomBackgrd(int nargs, char *args[], const string name)
    : Species(nargs, args, name),
      injectMode(conserving,conserving),
      injectCoef(0.0, 0.0), injectNHCoef(0.5, 0.5),
      numInjectStat(0,0), numInjectHomPerCellStat(0,0)
  {
    imodeMap.insert(Pair(conserving,"PN Conserved"));
    imodeMap.insert(Pair(statistic ,"   Statistic"));
    imodeMap.insert(Pair(noinject  ,"No Injection"));

    initParsing(nargs, args);
    paramParsing();
    checkParam();
  }

  NonHomBackgrd::~NonHomBackgrd()
  {}

  void NonHomBackgrd::initialise(NumberGenerator &draw, Scheduler &scheduler)
  {
    // First initialise abstract base class
    Species::initialise(draw, scheduler);

    // then load non homogeneous model from file
    // and init number of particles/cell to inject
    loadNonHomModel();

    particles.resize(numParticles);

    BEGIN_DEBUG_OUTPUT(20)
    HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
    VAR_DEBUG_OUTPUT2(numParticles,particles.size())
    END_DEBUG_OUTPUT

    ParticleListIter p=particles.begin();
    RVectorr dr = (domain.max-domain.min)/numCells;
    int ntot = 0;
#if (DIMR==3)
    for (int k=0; k<numCells(2); ++k)
#endif
      for (int j=0; j<numCells(1); ++j)
        for (int i=0; i<numCells(0); ++i) {
#if (DIMR==2)
          int npic = numParticlesPerCell(i,j);
          RVectorr rmin = domain.min+RVectorr(dr(0)*i,dr(1)*j);
#elif (DIMR==3)
          int npic = numParticlesPerCell(i,j,k);
          RVectorr rmin = domain.min+RVectorr(dr(0)*i,dr(1)*j,dr(2)*k);
#endif
          RVectorr rmax = rmin+dr;
//cout<<"i="<<i<<",j="<<j<<",rmin="<<rmin<<",rmax="<<rmax<<",npic="<<npic<<endl;
// drho/dx = 1/2((rho(x1,y0)-rho(x0,y0))/(x1-x0)+
//               (rho(x1,y1)-rho(x0,y1))/(x1-x0))
// drho/dy = 1/2((rho(x0,y1)-rho(x0,y0))/(y1-y0)+
//               (rho(x1,y1)-rho(x1,y0))/(y1-y0))
          for (int n=1; n<=npic; ++n, ++p) {
            ntot ++ ;
            for (int d=0; d<DIMR; ++d) {
              p->R(d) = draw.uniform(rmin(d), rmax(d));
            }
            for (int d=0; d<DIMV; ++d) {
              p->V(d) = draw.normal(u_s(d), v_s(d));
            }
          }
        }
    BEGIN_DEBUG_OUTPUT(20)
    HEADER_DEBUG_OUTPUT2(speciesName,"initialise")
    VAR_DEBUG_OUTPUT2(numParticles,ntot)
    END_DEBUG_OUTPUT
  }

  void NonHomBackgrd::printOn(ostream& os) const
  {
    Species::printOn(os);

    os
        << "Inject mode min         = " << map_elt(imodeMap, injectMode.min) << '\n'
        << "Inject mode max         = " << map_elt(imodeMap, injectMode.max) << '\n'
        << "Inject coef min (hom)   = " << injectCoef.min << '\n'
        << "Inject coef max (hom)   = " << injectCoef.max << '\n'
        << "NH Model Filename       = " << nhmFilename << '\n'
        << "Inject coef min (nh)    = " << injectNHCoef.min << '\n'
        << "Inject coef max (nh)    = " << injectNHCoef.max << '\n';
  }

  void NonHomBackgrd::fixBoundaryParticles(Scheduler &scheduler)
  {
    Loss = 0;
    for (int d=0; d<DIMR; ++d) {
      LossAtLB(d) = 0;
      LossAtUB(d) = 0;
    }
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
              RVectori idimr((i->R-domainMin)/(domainMax-domainMin)*numCells);
              BVectori ibndy;
              BVectori::iterator iter(ibndy.begin());
              for (int id=0; id<DIMR; ++id) {
                if (id != d) { // direction parallel boundary
                  *iter++ = (idimr(id) < LossAtLB(d).lbound(id) ?
                             LossAtLB(d).lbound(id) :
                             (idimr(id) > LossAtLB(d).ubound(id) ?
                              LossAtLB(d).ubound(id) : idimr(id) ) );
#if defined(BZ_DEBUG)
                  // Check whether ibndy is within limits
                  if (idimr(id) < LossAtLB(d).lbound(id) ||
                      idimr(id) > LossAtLB(d).ubound(id)) {
                    cout << "ERROR " << "idimr(" << id << ")=" << idimr(id)
                         << " but LossAtLB(" << d << ").ubound(" << id << ")="
                         << LossAtLB(d).ubound(id) << endl;
                  }
#endif
                }
              }
              LossAtLB(d)(ibndy) += 1;
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
              RVectori idimr((i->R-domainMin)/(domainMax-domainMin)*numCells);
              BVectori ibndy;
              BVectori::iterator iter(ibndy.begin());
              for (int id=0; id<DIMR; ++id) {
                if (id != d) { // direction parallel to boundary
                  *iter++ = (idimr(id) < LossAtUB(d).lbound(id) ?
                             LossAtLB(d).lbound(id) :
                             (idimr(id) > LossAtUB(d).ubound(id) ?
                              LossAtLB(d).ubound(id) : idimr(id) ) );
#if defined(BZ_DEBUG)
                  // Check whether ibndy is within limits
                  if (idimr(id) < LossAtUB(d).lbound(id) ||
                      idimr(id) > LossAtUB(d).ubound(id)) {
                    cout << "ERROR! " << "idimr(" << id << ")=" << idimr(id)
                         << " but LossAtUB(" << d << ").ubound(" << id << ")="
                         << LossAtLB(d).ubound(id) << endl;
                  }
#endif
                }
              }
              LossAtUB(d)(ibndy) += 1;
            }
            break;
          }
        }
        if (!lost && (i->R(d) < domainMin(d) || i->R(d) >= domainMax(d))) {
          cout << "ERROR"<< " d=" << d << " i->R="<< i->R(d)
               <<" domainMin="<<domainMin(d)<<" domainMax="<<domainMax(d)<<endl;
        }
      }
      if (lost) {
        sumv2After -= dot(i->V, i->V);
        i = particles.erase(i);
      } else {
        ++i;
      }
    } while (i != e);

    BEGIN_DEBUG_OUTPUT(30)
    HEADER_DEBUG_OUTPUT2(speciesName,"fixBoundaryParticles");
    for (int d=0; d<DIMR; ++d) {
      if (BC.min(d) == free_space)
        cout << "LossAtLB(" << d << ")=" << LossAtLB(d) << endl;
      if (BC.max(d) == free_space)
        cout << "LossAtUB(" << d << ")=" << LossAtUB(d) << endl;
    }
    for (int d=0; d<DIMR; ++d) {
      cout << "sum(LossAtLB(" << d << "))=" << sum(LossAtLB(d))
           << " "
           << "sum(LossAtUB(" << d << "))=" << sum(LossAtUB(d))
           << endl;
    }
    END_DEBUG_OUTPUT0

    BEGIN_DEBUG_OUTPUT(5)
    HEADER_DEBUG_OUTPUT2(speciesName,"fixBoundaryParticles")
    VAR_DEBUG_OUTPUT1(Loss)
    END_DEBUG_OUTPUT

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"fixBoundaryParticles")
    VAR_DEBUG_OUTPUT2(particles.size(),numParticles)
    END_DEBUG_OUTPUT

  }

  void NonHomBackgrd::injectParticles(NumberGenerator &draw, Scheduler &scheduler)
  {
    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"injectParticles")
    VAR_DEBUG_OUTPUT1(Loss)
    END_DEBUG_OUTPUT

    // Trick
    // if bulk velocity != 0 and
    // number of particles is conserved then
    // reinject what is lost at one side into the opposite side
    for (int d=0; d<DIMR; ++d) {
      blitz::Array<int,DIMR-1> loss_min(LossAtLB(d));
      blitz::Array<int,DIMR-1> loss_max(LossAtUB(d));
      switch (injectMode.min(d)) {
      case statistic:
        LossAtLB(d) = numInjectPerCellAtLB(d);
        break;
      case conserving:
        if (u_s(d)!=0.0)
          LossAtLB(d) = loss_max;
        break;
      case noinject:
        LossAtLB(d) = 0;
        break;
      }
      switch (injectMode.max(d)) {
      case statistic:
        LossAtUB(d) = numInjectPerCellAtUB(d);
        break;
      case conserving:
        if (u_s(d)!=0.0)
          LossAtUB(d) = loss_min;
        break;
      case noinject:
        LossAtUB(d) = 0;
        break;
      }
    }

    for (int d=0; d<DIMR; ++d) {
      if (injectMode.min(d) != noinject) {
        // (1+injectcoef) * constant_density + injectNHcoef * perturb_density
        LossAtLB(d) = cast<int>((1.0+injectCoef.min(d))*
                                numInjectHomPerCellStat.min(d)+
                                injectNHCoef.min(d)*
                                (LossAtLB(d)-numInjectHomPerCellStat.min(d)) +
                                0.5);
      }
      Loss.min(d) = sum(LossAtLB(d));
      if (injectMode.max(d) != noinject) {
        LossAtUB(d) = cast<int>((1.0+injectCoef.max(d))*
                                numInjectHomPerCellStat.max(d)+
                                injectNHCoef.max(d)*
                                (LossAtUB(d)-numInjectHomPerCellStat.max(d)) +
                                0.5);
      }
      Loss.max(d) = sum(LossAtUB(d));
    }

    BEGIN_DEBUG_OUTPUT(30)
    HEADER_DEBUG_OUTPUT2(speciesName,"injectParticles");
    for (int d=0; d<DIMR; ++d) {
      if (BC.min(d) == free_space)
        cout << "LossAtLB(" << d << ")=" << LossAtLB(d) << endl;
      if (BC.max(d) == free_space)
        cout << "LossAtUB(" << d << ")=" << LossAtUB(d) << endl;
    }
    for (int d=0; d<DIMR; ++d) {
      cout << "sum(LossAtLB(" << d << "))=" << sum(LossAtLB(d))
           << " "
           << "sum(LossAtUB(" << d << "))=" << sum(LossAtUB(d))
           << endl;
    }
    END_DEBUG_OUTPUT0

    RVectorr dr = (domain.max-domain.min)/numCells;

    numToInject = sum(Loss.min+Loss.max);

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"injectParticles")
    VAR_DEBUG_OUTPUT1(Loss)
    VAR_DEBUG_OUTPUT1(numToInject)
    END_DEBUG_OUTPUT

    ParticleList p(numToInject);
    ParticleListIter ip = p.begin();
    int numInserted = 0;
    for (int d=0; d<DIMR; ++d) { // For all spatial dimensions
#if (DIMR==3)
      for (int k=0; k<(d==2? 1 : numCells(2)); ++k) {
#endif
        for (int j=0; j<(d==1? 1 : numCells(1)); ++j) {
          for (int i=0; i<(d==0? 1 : numCells(0)); ++i) {
#if (DIMR==2)
            RVectori idimr(i,j);
#elif (DIMR==3)
            RVectori idimr(i,j,k);
#endif
            BVectori ibndy;
            BVectori::iterator iter(ibndy.begin());
            for (int id=0; id<DIMR; ++id) {
              if (id != d) { // direction in boundary
                *iter++ = idimr(id);
              }
            }
            //cout << "idimr=" << idimr << " ibndy=" << ibndy << endl;
#if (DIMR==2)
            RVectorr rmin = domain.min+RVectorr(dr(0)*i,dr(1)*j);
#elif (DIMR==3)
            RVectorr rmin = domain.min+RVectorr(dr(0)*i,dr(1)*j,dr(2)*k);
#endif
            RVectorr rmax = rmin+dr;
            //cout << "rmin=" << rmin << " rmax=" << rmax << endl;
            for (int n=0; n<LossAtLB(d)(ibndy); ++n) { // Lower boundary
              for (int pd=0; pd<DIMR; ++pd) {
                real a = (pd == d ? domain.min(pd) : rmin(pd));
                real b = (pd == d ? domain.max(pd) : rmax(pd));
                if (pd == d) { // perpendicular direction to the face
                  ip->V(pd) = draw.flux(u_s(pd), v_s(pd));
                  ip->R(pd) = draw.uniform(a, std::min(a+ip->V(pd)*dt, b));
                } else { // parallel directions to the face
                  ip->R(pd) = draw.uniform(a, b);
                  ip->V(pd) = draw.normal(u_s(pd), v_s(pd));
                }
                if (ip->R(pd) < domain.min(pd) || ip->R(pd) >= domain.max(d)) {
                  cout << "ERROR LB"<< " d,pd="<<d<<","<<pd<<" i->R="<< ip->R(pd)
                       <<" domain.min="<<domain.min(pd)<<" domain.max="<<domain.max(pd)<<endl;
                }
              }
#if (DIMR==2)
              ip->V(DIMV-1) = draw.normal(u_s(DIMV-1), v_s(DIMV-1));
#endif
              sumv2After += dot(ip->V, ip->V);
              ++ip;
              numInserted++;
            }
            //cout << "numToInject=" << numToInject << "numInserted=" << numInserted << endl;
#if 1
            for (int n=0; n<LossAtUB(d)(ibndy); ++n) { // Upper boundary
              for (int pd=0; pd<DIMR; ++pd) {
                real a = (pd == d ? domain.min(pd) : rmin(pd));
                real b = (pd == d ? domain.max(pd) : rmax(pd));
                if (pd == d) { // perpendicular direction to the face
                  ip->V(pd) = -draw.flux(-u_s(pd), v_s(pd));
                  //ip->R(pd) = draw.uniform(std::max(a, b+ip->V(pd)*dt), b);
                  real l=std::max(a, b+ip->V(pd)*dt);
                  real u=b;
                  ip->R(pd) = draw.uniform(l,u);
                  if (ip->R(pd) >= domain.max(pd))
                    cout<<"draw.uniform1("<<std::max(a, b+ip->V(pd)*dt)
                        << ","<<b<<")="<<ip->R(pd)
                        <<" (ip->R(pd)-b)="<< b-ip->R(pd) <<endl;
                  if (ip->R(pd) >= u)
                    cout<<"draw.uniform2("<<l<<","<<u<<")="<<ip->R(pd)<<" x-u="<<ip->R(pd)-u<<endl;
                } else { // parallel directions
                  ip->R(pd) = draw.uniform(a, b);
                  ip->V(pd) = draw.normal(u_s(pd), v_s(pd));
                }
              }
              for (int pd=0; pd<DIMR; ++pd) {
                if (ip->R(pd) < domain.min(pd) || ip->R(pd) >= domain.max(d)) {
                  cout<<"ERROR UB"<<" ibndy="<<ibndy<<" n="<<n<<" d,pd="<<d<<","<<pd<<" i->R="<< ip->R
                      <<" domain.min="<<domain.min<<" domain.max="<<domain.max<<endl;
                }
              }
#if (DIMR==2)
              ip->V(DIMV-1) = draw.normal(u_s(DIMV-1), v_s(DIMV-1));
#endif

              sumv2After += dot(ip->V, ip->V);
              ++ip;
              numInserted++;
            }
            //cout << "numToInject=" << numToInject << "numInserted=" << numInserted << endl;
#endif
          }
        }
#if (DIMR==3)
      }
#endif
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

  void NonHomBackgrd::loadNonHomModel()
  {
    // Open non-homogeneous model file
    std::ifstream ifs(nhmFilename.c_str());
    if (nhmFilename.empty() || ifs.bad()) {
      ostringstream os;
      os << "Unable to open file:  " << nhmFilename;
      throw ClassException("NonHomBackgrd", os.str());
    }

    // file generated using matlab function write_nhbgrd_model
    for (int d=0; d<DIMR; d++) {
      ifs >> gridSize(d) ;
    }
    relDens.resize(gridSize);

#if (DIMR==3)
    for (int k=0; k<gridSize(2); ++k)
#endif
      for (int j=0; j<gridSize(1); ++j)
        for (int i=0; i<gridSize(0); ++i) { // Fastest changing variable saved
          if (ifs.eof() || ifs.bad()) {
            ostringstream os;
            os << "Premature end of file on " << nhmFilename;
            throw ClassException("NonHomBackgrd", os.str());
          }
#if (DIMR==2)
          ifs >> relDens(i,j);
#elif (DIMR==3)
          ifs >> relDens(i,j,k);
#endif
        }

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"loadNonHomModel")
    VAR_DEBUG_OUTPUT1(nhmFilename)
    VAR_DEBUG_OUTPUT1(gridSize)
    VAR_DEBUG_OUTPUT2(min(relDens),max(relDens))
    END_DEBUG_OUTPUT

    // n-1 intervals for n-points grid
    numCells = gridSize-1;
    numParticlesPerCell.resize(numCells);
    cellRelDens.resize(numCells);
    RVectorr deltar = (domain.max-domain.min)/numCells;
    real dV = product(deltar);

    Range I(0,numCells(0)-1);
    Range J(0,numCells(1)-1);
#if (DIMR==2)
    cellRelDens(I,J) = 0.25*(relDens(I,J)  + relDens(I+1,J) +
                             relDens(I,J+1)+ relDens(I+1,J+1));
    numParticlesPerCell(I,J) = cast<int>(lambda*cellRelDens(I,J)*dV + 0.5);
#elif (DIMR==3)
    Range K(0,numCells(2)-1);
    cellRelDens(I,J,K) = 0.125*(relDens(I,J,K)+ relDens(I+1,J,K)+
                                relDens(I,J+1,K)+ relDens(I+1,J+1,K)+
                                relDens(I,J,K+1)+ relDens(I+1,J,K+1)+
                                relDens(I,J+1,K+1)+relDens(I+1,J+1,K+1));
    numParticlesPerCell(I,J,K) = cast<int>(lambda*cellRelDens(I,J,K)*dV + 0.5);
#endif

    numParticles = blitz::sum(numParticlesPerCell);

    BEGIN_DEBUG_OUTPUT(10)
    int numParticlesHom = int(lambda*product(domain.max-domain.min)+0.5);
    HEADER_DEBUG_OUTPUT2(speciesName,"loadNonHomModel")
    VAR_DEBUG_OUTPUT2(numParticles,numParticlesHom)
    END_DEBUG_OUTPUT

    VVectorr nu(u_s/v_s);
    RVectorr Vmin, Vmax;
    for (int d=0; d<DIMR; ++d) {
      Vmin(d)=M_1_SQRT2PI*v_s(d)*exp(double(-0.5*nu(d)*nu(d)))+
              0.5*u_s(d)*(1.0+erf(double(nu(d)/M_SQRT2)));
      Vmax(d)=M_1_SQRT2PI*v_s(d)*exp(double(-0.5*nu(d)*nu(d)))-
              0.5*u_s(d)*(1.0-erf(double(nu(d)/M_SQRT2)));
    }


    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(speciesName,"loadNonHomModel")
    VAR_DEBUG_OUTPUT2(Vmin,Vmax)
    END_DEBUG_OUTPUT

    numInjectHomPerCellStat = Boundary<int,DIMR>(lambda*dV/deltar*Vmin*dt,
                              lambda*dV/deltar*Vmax*dt);

    BEGIN_DEBUG_OUTPUT(5)
    HEADER_DEBUG_OUTPUT2(speciesName,"loadNonHomModel")
    VAR_DEBUG_OUTPUT1(numInjectHomPerCellStat)
    END_DEBUG_OUTPUT

    Range all(Range::all());
    RVectori lbnd(cellRelDens.lbound());
    RVectori ubnd(cellRelDens.ubound());
    for (int d=0; d<DIMR; ++d) { // for all physical dimensions
      int perp1 = numCells((d+1) % DIMR);
#if (DIMR==2)
      numInjectPerCellAtLB(d).resize(perp1);
      numInjectPerCellAtLB(d) = cast<int>(
                                  lambda*cellRelDens(lbnd(d), all)*
                                  dV/deltar(d)*Vmin(d)*dt+0.5);
      numInjectPerCellAtUB(d).resize(perp1);
      numInjectPerCellAtUB(d) = cast<int>(
                                  lambda*cellRelDens(ubnd(d), all)*
                                  dV/deltar(d)*Vmax(d)*dt+0.5);
      LossAtLB(d).resize(perp1);
      LossAtUB(d).resize(perp1);
      // trick to (circular) permute on the left the dimensions
      // which let access always the same indices
      cellRelDens.transposeSelf(secondDim, firstDim);
#elif (DIMR==3)
      int perp2 = numCells((d+2) % DIMR);
      numInjectPerCellAtLB(d).resize(perp1, perp2);
      numInjectPerCellAtLB(d) = cast<int>(
                                  lambda*cellRelDens(lbnd(d), all, all)*
                                  dV/deltar(d)*Vmin(d)*dt+0.5);
      numInjectPerCellAtUB(d).resize(perp1, perp2);
      numInjectPerCellAtUB(d) = cast<int>(
                                  lambda*cellRelDens(ubnd(d), all, all)*
                                  dV/deltar(d)*Vmin(d)*dt+0.5);
      LossAtLB(d).resize(perp1, perp2);
      LossAtUB(d).resize(perp1, perp2);
      cellRelDens.transposeSelf(secondDim, thirdDim, firstDim);
#endif
    }

    RVectorr injectMin, injectMax;
    for (int d=0; d<DIMR; ++d) {
      injectMin(d) = sum(numInjectPerCellAtLB(d));
      injectMax(d) = sum(numInjectPerCellAtUB(d));
    }
    numInjectStat = Boundary<int,DIMR>(injectMin,injectMax);

    BEGIN_DEBUG_OUTPUT(5)
    HEADER_DEBUG_OUTPUT2(speciesName,"loadNonHomModel")
    VAR_DEBUG_OUTPUT1(numInjectStat)
    END_DEBUG_OUTPUT
  }

  void NonHomBackgrd::initParsing(int nargs, char *args[])
  {
    registerClass(speciesName);
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

#if 0
    const char rveci[] =  "int[" DIMRSTR "]";
    const char rvecr[] = "real[" DIMRSTR "]";
#endif

    using parser::types::intVect;
    using parser::types::realVect;
    using parser::types::charStr;

    insertOption(_imodemin, "imodemin", intVect ,
                 "Lower boundary inject mode", Any(injectMode.min));
    insertOption(_imodemax, "imodemax", intVect ,
                 "Upper boundary inject mode", Any(injectMode.max));
    insertOption(_icoefmin, "icoefmin", realVect,
                 "Homogeneous lower boundary inject coefficient",
                 Any(injectCoef.min));
    insertOption(_icoefmax, "icoefmax", realVect,
                 "Homogeneous upper boundary inject coefficient",
                 Any(injectCoef.max));

    insertOption(_nhmfilename, "nhmfilename", charStr,
                 "Filename for non-homogeneous model", Any(nhmFilename));

    insertOption(_inhcoefmin, "inhcoefmin", realVect,
                 "Non homogeneous lower boundary inject coefficient",
                 Any(injectNHCoef.min));
    insertOption(_inhcoefmax, "inhcoefmax", realVect,
                 "Non homogeneous upper boundary inject coefficient",
                 Any(injectNHCoef.max));
  }

  void NonHomBackgrd::paramParsing()
  {
    parseOption(_imodemin, injectMode.min);
    parseOption(_imodemax, injectMode.max);
    parseOption(_icoefmin, injectCoef.min);
    parseOption(_icoefmax, injectCoef.max);

    parseOption(_nhmfilename, nhmFilename);

    parseOption(_inhcoefmin, injectNHCoef.min);
    parseOption(_inhcoefmax, injectNHCoef.max);
  }

  void NonHomBackgrd::checkParam() const
  {
    checkMap(_imodemin, imodeMap, injectMode.min);
    checkMap(_imodemax, imodeMap, injectMode.max);
  }

#undef ID

}
