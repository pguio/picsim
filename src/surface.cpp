/**************************************************************************
 *
 * $Id: surface.cpp,v 1.45 2011/03/26 15:36:08 patrick Exp $
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

#include <surface.h>

namespace picsim {

  using blitz::all;
  using blitz::dot;
  using blitz::pow2;
  using blitz::where;

  using parser::map_elt;

  using std::ostream;

#define ID "$Id: surface.cpp,v 1.45 2011/03/26 15:36:08 patrick Exp $"

  using mudfas::DEFAULT_GRIDSIZE;

  const RVectorr DEFAULT_SURF_MIN = (DEFAULT_GRIDSIZE-1)/2-(DEFAULT_GRIDSIZE-1)/4;
  const RVectorr DEFAULT_SURF_MAX = (DEFAULT_GRIDSIZE-1)/2+(DEFAULT_GRIDSIZE-1)/4;


  ostream& operator<<(ostream& os, const Surface &s)
  {
#if (DIMR==3)
    const char *axis[]= { "x", "y", "z"
                        };
#endif

#if defined(HAVE_MPI)

    if (s.rankProc == masterProc)
#endif

      return os
             << "type                 = " << map_elt(s.typeMap, s.type) << '\n'
#if (DIMR==3)
             << (s.type == Surface::cylindrical ?  "Cylindre direction   = " : "")
             << (s.type == Surface::cylindrical ?  axis[s.cylDir] : "")
             << (s.type == Surface::cylindrical ?  "\n" : "")
#endif
             << "property             = " << map_elt(s.propertyMap, s.property) << '\n'
             << "surf rho             = " << s.surfRho << '\n'
             << "rmin                 = " << s.limits.min << '\n'
             << "rmax                 = " << s.limits.max << '\n';
#if defined(HAVE_MPI)

    else
      return os;
#endif
  }


  Surface::Surface(int nargs, char *args[], const string name)
    : Parser(nargs, args),
      surfaceName(name), type(cartesian), cylDir(zDir), property(absorbing),
      surfRho(1.0), limits(DEFAULT_SURF_MIN,DEFAULT_SURF_MAX)
  {
    typeMap.insert(Pair(cartesian  ,   "  cartesian"));
    typeMap.insert(Pair(spherical  ,   "  spherical"));
#if (DIMR==3)

    typeMap.insert(Pair(cylindrical,   "cylindrical"));
#endif

    dirMap.insert(Pair(xDir  ,   "x"));
    dirMap.insert(Pair(yDir  ,   "y"));
    dirMap.insert(Pair(zDir  ,   "z"));

    propertyMap.insert(Pair(absorbing, "  absorbing"));
    propertyMap.insert(Pair(emitting , "   emitting"));

    initParsing(nargs, args);
    paramParsing();

    checkParam();
  }

  Surface::~Surface()
  {}

  void Surface::initialise(const Scheduler &scheduler)
  {
    using blitz::tensor::i;
    using blitz::tensor::j;
#if (DIM==3)
    using blitz::tensor::k;
#endif

    Boundary<real,DIMR> domain(scheduler.getDomainBoundary());
    RVectori gridSize(scheduler.getGridSize());

    if (type == spherical || type == cylindrical) {
#if !defined(_AIX)
      center = 0.5*(limits.min+limits.max);
      scale  = 2.0/(limits.max-limits.min);
#else

      for (int d=0; d<DIMR; ++d) {
        center(d) = 0.5*(limits.min(d)+limits.max(d));
        scale(d)  = 2.0/(limits.max(d)-limits.min(d));
      }
#endif

    }

    rhos.resize(gridSize);
    rhos = 0;

    mudfas::real mn, mx;
    int n1;
#if (DIM==2)

    Field x(gridSize), y(gridSize);
    mn = domain.min(0);
    mx = domain.max(0);
    n1 = gridSize(0)-1;
    x = (mx-mn)*i/n1+mn+0.0*j;

    mn = domain.min(1);
    mx = domain.max(1);
    n1 = gridSize(1)-1;
    y = 0.0*i+(mx-mn)*j/n1+mn;
#elif (DIM==3)

    Field x(gridSize), y(gridSize), z(gridSize);

    mn = domain.min(0);
    mx = domain.max(0);
    n1 = gridSize(0)-1;
    x = (mx-mn)*i/n1+mn+0.0*j+0.0*k;

    mn = domain.min(1);
    mx = domain.max(1);
    n1 = gridSize(1)-1;
    y = 0.0*i+(mx-mn)*j/n1+mn+0.0*k;

    mn = domain.min(2);
    mx = domain.max(2);
    n1 = gridSize(2)-1;
    z = 0.0*i+0.0*j+(mx-mn)*k/n1+mn;
#endif

    Field rho0(gridSize);
    rho0 = surfRho;

    switch (type) {
    case  cartesian:
#if (DIM==2)

      rhos = where( x >= limits.min(0) && x <= limits.max(0) &&
                    y >= limits.min(1) && y <= limits.max(1), rho0, rhos);
#elif (DIM==3)

      rhos = where( x >= limits.min(0) && x <= limits.max(0) &&
                    y >= limits.min(1) && y <= limits.max(1) &&
                    z >= limits.min(2) && z <= limits.max(2), rho0, rhos);
#endif

      break;

    case  spherical:
#if (DIM==2)

      rhos = where( pow2((x-center(0))*scale(0)) +
                    pow2((y-center(1))*scale(1)) <= 1.0, rho0, rhos);
#elif (DIM==3)

      rhos = where( pow2((x-center(0))*scale(0)) +
                    pow2((y-center(1))*scale(1)) +
                    pow2((z-center(2))*scale(2)) <= 1.0, rho0, rhos);
#endif

      break;

#if (DIM==3)

    case cylindrical:
      switch (cylDir) {
      case xDir:
        rhos = where( x >= limits.min(0) && x <= limits.max(0) &&
                      pow2((y-center(1))*scale(1)) +
                      pow2((z-center(2))*scale(2)) <= 1.0, rho0, rhos);
        break;
      case yDir:
        rhos = where( y >= limits.min(1) && y <= limits.max(1) &&
                      pow2((x-center(0))*scale(0)) +
                      pow2((z-center(2))*scale(2)) <= 1.0, rho0, rhos);
        break;
      case zDir:
        rhos = where( z >= limits.min(2) && z <= limits.max(2) &&
                      pow2((x-center(0))*scale(0)) +
                      pow2((y-center(1))*scale(1)) <= 1.0, rho0, rhos);
        break;
      }
      break;
#endif

    }
  }

  void Surface::processSpecies(SpeciesHandler &species)
  {
    SpeciesListIter i=species.begin(), e=species.end();
    for ( ; i!=e; ++i) {
      processParticles(**i);
    }
  }

  void Surface::processDensity(Field &rho)
  {
    rho += rhos;
  }

  void Surface::processParticles(Species &species)
  {

    ParticleListIter i=species.begin(), e=species.end();
    int numLost  = 0;
    real v2Lost = 0;

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT1("Surface::processParticles")
    VAR_DEBUG_OUTPUT2(species.size(),species.getParticleNumber())
    END_DEBUG_OUTPUT

    switch (type) {
    case cartesian:
      do {
#if defined(__GNUC__) && (__GNUC__ < 3)
        if ( all( blitz::operator>= ((*i).R, limits.min) &&
                  blitz::operator<= ((*i).R, limits.max ) ) ) {
#else
        if ( all( (*i).R >= limits.min && (*i).R <= limits.max ) ) {
#endif
          ++numLost;
          v2Lost += dot((*i).V, (*i).V);
          i = species.erase(i);
        } else {
          ++i;
        }
      }
      while (i != e)
        ;
      break;

    case spherical:
      do {
        RVectorr r(((*i).R-center)*scale);
        if ( dot(r,r) <= 1.0 ) {
          ++numLost;
          v2Lost += dot((*i).V, (*i).V);
          i = species.erase(i);
        } else {
          ++i;
        }
      } while (i != e);
      break;

#if (DIMR==3)

    case cylindrical:
      do {
        RVectorr r(((*i).R-center)*scale);
        r(cylDir) = 0.0;
        if ( dot(r,r) <= 1.0 &&
             (*i).R(cylDir) >= limits.min(cylDir) &&
             (*i).R(cylDir) <= limits.max(cylDir) ) {
          ++numLost;
          v2Lost += dot((*i).V, (*i).V);
          i = species.erase(i);
        } else {
          ++i;
        }
      } while (i != e);
      break;
#endif

    }

    species.sync(numLost, v2Lost);

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT2(surfaceName,"processParticles")
    VAR_DEBUG_OUTPUT2(species.size(),species.getParticleNumber())
    VAR_DEBUG_OUTPUT1(numLost)
    END_DEBUG_OUTPUT
  }

  void Surface::initParsing(int nargs, char *args[])
  {
    registerClass(surfaceName);
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    string prefix(surfaceName);
    prefix += "::";
    setPrefix(prefix.c_str());

    parseLevelDebugOption("dl");

    using parser::types::integer;
    using parser::types::real;
    using parser::types::realVect;

#if 0
    const char rvecr[] = "real[" DIMRSTR "]";
#endif

    insertOption(_type    , "type"    , integer , "Type of the surface"      , Any(type));
    insertOption(_cylDir  , "cylDir"  , integer , "Direction of the cylinder", Any(cylDir));

    insertOption(_property, "property", integer , "Property of the surface"  , Any(property));
    insertOption(_surfRho , "rho"     , real    , "Density of the surface"   , Any(surfRho));

    insertOption(_rmin    , "rmin"    , realVect, "Minimum of the surface"   , Any(limits.min));
    insertOption(_rmax    , "rmax"    , realVect, "Minimum of the surface"   , Any(limits.max));
  }

  void Surface::paramParsing()
  {
    parseOption(_type    , type);
    parseOption(_cylDir  , cylDir);
    parseOption(_property, property);
    parseOption(_surfRho , surfRho);
    parseOption(_rmin    , limits.min);
    parseOption(_rmax    , limits.max);
  }

  void Surface::checkParam() const
  {
    checkMap(_type    , typeMap    , type);
    checkMap(_cylDir  , dirMap     , cylDir);
    checkMap(_property, propertyMap, property);
  }

#undef ID

}
