/**************************************************************************
 *
 * $Id: moments.cpp,v 1.29 2011/03/26 15:36:08 patrick Exp $
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

#include <moments.h>

namespace picsim {

  using blitz::all;
  using blitz::floor;
  using blitz::pow2;

  using std::ostream;

#define ID "$Id: moments.cpp,v 1.29 2011/03/26 15:36:08 patrick Exp $"

  ostream& operator<<(ostream& os, const Moments &m)
  {
    return os
           << "Moments    name    = " << m.momentName << '\n'
           << "  Species names    = " << m.speciesNames << '\n'
           << "  Start            = " << m.start << '\n'
           << "  End              = " << m.end << '\n'
           << "  Stride           = " << m.stride << '\n'
           << "  Velocity flag    = " << m.uFlag << '\n'
           << "  Temperature flag = " << m.TFlag << '\n';
  }

  Moments::Moments(int nargs, char *args[], const string name)
    : Parser(nargs, args), momentName(name),
      start(0.0), end(0.0), stride(0.0), uFlag(false), TFlag(false)
  {
    initParsing(nargs, args);
    paramParsing();
  }

  Moments::~Moments()
  {}

  void Moments::initialise(const Scheduler & scheduler)
  {
    for (int d=0; d<DIMR; ++d) {
      gridSize(d) = static_cast<int>((end(d)-start(d))/stride(d)+0.5);
      if ( std::abs(end(d)-start(d)-gridSize(d)*stride(d)) >
           blitz::epsilon(static_cast<real>(1.0))*end(d) ) {
        std::ostringstream s;
        s.precision(15);
        s << "Improper start,stride,end spec: (end-start)/stride=("
          << end(d) << "-" << start(d) << ")/" << stride(d) << " not integral.";
        throw ClassException("Moments", s.str());
      }
    }
    grid.resize(gridSize);

    rhoFactor = 1.0/scheduler.getRefLambda()/blitz::product(stride);
    for (int d=0; d<DIMV; ++d) {
      n.resize(gridSize);
      if ( uFlag(d) || TFlag(d) )
        u(d).resize(gridSize);
      if ( TFlag(d) )
        T(d).resize(gridSize);
    }
  }

  void Moments::computeMoments(SpeciesHandler &species, Field &dens,
                               VectorField &vel, VectorField &Temp)
  {
    for (int d=0; d<DIMV; ++d) {
      if ( uFlag(d) || TFlag(d) )
        u(d) = 0.0;
      if ( TFlag(d) )
        T(d) = 0.0;
    }
#if defined(HAVE_MPI)
    SpeciesListIter i = species.begin();
    if ( isSpeciesinMoments((**i).name())) {
      computeMoments(**i);
    }
    MPI_Reduce(n.data(), dens.data(), n.size(), MPI_COORDINATE,
               MPI_SUM, mpiRoot, MPI_COMM_WORLD);
    for (int d=0; d<DIMV; ++d) {
      if ( uFlag(d) )
        MPI_Reduce(u(d).data(), vel(d).data(), u(d).size(), MPI_COORDINATE,
                   MPI_SUM, mpiRoot, MPI_COMM_WORLD);
      if ( TFlag(d) )
        MPI_Reduce(T(d).data(), Temp(d).data(), T(d).size(), MPI_COORDINATE,
                   MPI_SUM, mpiRoot, MPI_COMM_WORLD);
    }
#else
    for (SpeciesListIter i=species.begin(), e=species.end(); i!=e; ++i) {
      if ( isSpeciesinMoments((**i).name()) ) {
        computeMoments(**i);
      }
    }
    dens = n;
    for (int d=0; d<DIMV; ++d) {
      if ( uFlag(d) )
        vel(d) = u(d);
      if ( TFlag(d) )
        Temp(d) = T(d);
    }
#endif

  }

  const std::string &Moments::name() const
  {
    return momentName;
  }

  bool Moments::isVel(int d) const
  {
    return uFlag(d);
  }

  bool Moments::isTemp(int d) const
  {
    return TFlag(d);
  }

  mudfas::real Moments::getStart(int d) const
  {
    return start(d);
  }

  mudfas::real Moments::getEnd(int d) const
  {
    return end(d);
  }

  mudfas::real Moments::getStride(int d) const
  {
    return stride(d);
  }

  RVectori Moments::dim() const
  {
    return gridSize;
  }

  int Moments::dim(int d) const
  {
    return gridSize(d);
  }

  bool Moments::isSpeciesinMoments(const string &name) const
  {
    NameListConstIter i=speciesNames.begin(), e=speciesNames.end();
    for ( ; i!=e; ++i) {
      if (name == *i)
        return true;
    }
    return false;
  }

  void Moments::computeMoments(Species &species)
  {
    // Sort particles on grid
    RVectord a(gridSize/(end-start));
    for (ParticleListConstIter i=species.begin(), e=species.end(); i!=e; ++i) {
#if defined(__GNUC__) && (__GNUC__ < 3)
      if ( all( blitz::operator>= ((*i).R, start) &&
                blitz::operator< ((*i).R, end) ) ) {
#else
      if (all((*i).R >= start) && all((*i).R < end)) {
#endif
        RVectord rindex(a*((*i).R-start));
        RVectori index(floor(rindex));
        for (int d=0; d<DIMR; ++d) {
          index(d) = (index(d) > gridSize(d)-1 ? gridSize(d)-1 : index(d));
        }
        grid(index).push_back(&(*i));
      }
    }

    real mu_s = species.getmu_s();

    for (int d = 0; d<DIMV; ++d) { // For each velocity dimension
      if ( uFlag(d) || TFlag(d) ) {
        for (int i=0; i<gridSize(0); ++i) {
          for (int j=0; j<gridSize(1); ++j) {
#if (DIMR==2)
            RVectori index(i, j);
#elif (DIMR==3)

            for (int k=0; k<gridSize(2); ++k) {
              RVectori index(i, j, k);
#endif

            int numParticles = grid(index).size();
            ParticlePListConstIter pi, pe;

            // density
            n(index) = numParticles*rhoFactor;

            // mean drift velocity along the d-th velocity dimension
            if ( uFlag(d) || TFlag(d) ) {
              for (pi=grid(index).begin(), pe=grid(index).end(); pi!=pe; ++pi)
                u(d)(index) += (**pi).V(d);
              if (numParticles > 0)
                u(d)(index) /= numParticles;
            }

            // temperature along the d-th velocity dimension
            if ( TFlag(d) ) {
              for (pi=grid(index).begin(), pe=grid(index).end(); pi!=pe; ++pi)
                T(d)(index) += pow2((**pi).V(d)-u(d)(index));
              if (numParticles > 1)
                T(d)(index) /= (numParticles-1)/mu_s;
            }

#if (DIMR==3)

          }
#endif

        }
      }
    }
  }
  for (int i=0; i<gridSize(0); ++i)
    for (int j=0; j<gridSize(1); ++j)
#if (DIMR==2)

      grid(i, j).clear();
#elif (DIMR==3)

      for (int k=0; k<gridSize(2); ++k)
        grid(i, j, k).clear();
#endif

}

void Moments::initParsing(int nargs, char *args[])
{
  registerClass("Moments");
  registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

  string prefix(momentName);
  prefix += "::";
  setPrefix(prefix.c_str());

  parseLevelDebugOption("dl");

  using parser::types::boolVect;
  using parser::types::stringVect;
  using parser::types::realVect;

#if 0
  const char rvecr[] = "real[" DIMRSTR "]";
  const char vvecb[] = "bool[" DIMVSTR "]";
#endif

  insertOption(_speciesNames, "species", stringVect, "List of the species names", Any(speciesNames));

  insertOption(_start       , "start"  , realVect  , "Subvolume start"          , Any(start));
  insertOption(_end         , "end"    , realVect  , "Subvolume end"            , Any(end));
  insertOption(_stride      , "stride" , realVect  , "Subvolume stride"         , Any(stride));

  insertOption(_vel         , "vel"    , boolVect  , "u=<v> flag"               , Any(uFlag));
  insertOption(_Temp        , "Temp"   , boolVect  , "T=<|v-u|^2> flag"         , Any(TFlag));
}

void Moments::paramParsing()
{
  parseOption(_speciesNames, speciesNames);

  parseOption(_start       , start);
  parseOption(_end         , end);
  parseOption(_stride      , stride);

  parseOption(_vel         , uFlag);
  parseOption(_Temp        , TFlag);
}

#undef ID

}
