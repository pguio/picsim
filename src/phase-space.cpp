/**************************************************************************
 *
 * $Id: phase-space.cpp,v 1.51 2016/06/02 17:05:13 patrick Exp $
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

#include <phase-space.h>

namespace picsim {

  using std::ostream;

#define ID "$Id: phase-space.cpp,v 1.51 2016/06/02 17:05:13 patrick Exp $"

  ostream& operator<<(ostream& os, const PhaseSpace &p)
  {
    return os
           << "PhaseSpace name = " << p.spaceName << '\n'
           << "  Species names = " << p.speciesNames << '\n'
           << "  Start         = " << p.start << '\n'
           << "  End           = " << p.end << '\n'
           << "  Stride        = " << p.stride << '\n'
           << "  Integrated    = " << p.integrate << '\n'
           << "  Speed         = " << p.speed << '\n';
  }

  PhaseSpace::PhaseSpace(int nargs, char *args[], const string name)
    : Parser(nargs, args), spaceName(name),
      start(0.0), end(0.0), stride(0.0), integrate(false), speed(0.0),
      spaceNumDim(0), spaceDim(1)
  {
    initParsing(nargs, args);
    paramParsing();
  }

  PhaseSpace::~PhaseSpace()
  {}


  void PhaseSpace::initialise()
  {
    for (int d=0; d<DIMRV; ++d) {
      if (! integrate(d)) {
        spaceDim(d) = static_cast<int>((end(d)-start(d))/stride(d)+0.5);
        if ( std::abs(end(d)-start(d)-spaceDim(d)*stride(d)) >
             blitz::epsilon(static_cast<real>(1.0))*end(d) ) {
          std::ostringstream s;
          s.precision(8);
          s << "Improper start,stride,end spec:\n\t --> (end-start)/stride=('"
            << end(d) << "'-'" << start(d) << "')/" << stride(d)
            << "=" << (end(d)-start(d))/stride(d)
            << "\n\t --> not integral.";
          throw ClassException("PhaseSpace", s.str());
        }
        ++spaceNumDim;
      }
    }

    BEGIN_DEBUG_MASTER_OUTPUT(10)
    HEADER_DEBUG_OUTPUT1("PhaseSpace::initialise")
    VAR_DEBUG_OUTPUT1(spaceDim)
    VAR_DEBUG_OUTPUT1(spaceNumDim)
    END_DEBUG_MASTER_OUTPUT
  }


  void PhaseSpace::computePhaseSpace(const Scheduler &scheduler,
                                     const SpeciesHandler &species,
                                     RVArrayr &space)
  {
    real t = scheduler.getCurrentTime();
#if defined(HAVE_MPI)
    space.resize(spaceDim);
    RVArrayr spaces(space.shape());
    spaces = 0.0;
    SpeciesListConstIter i = species.begin();
    if ( isSpeciesinSpace((**i).name())) {
      integrateParticles(t, **i, spaces);
    }
    MPI_Reduce(spaces.data(), space.data(), spaces.size(), MPI_COORDINATE,
               MPI_SUM, mpiRoot, MPI_COMM_WORLD);
#else

    space.resize(spaceDim);
    space = 0.0;
    for (SpeciesListConstIter i=species.begin(), e=species.end(); i!=e; ++i) {
      if ( isSpeciesinSpace((*i)->name()) ) {
        integrateParticles(t, **i, space);
      }
    }
#endif

  }

  const std::string & PhaseSpace::name() const
  {
    return spaceName;
  }

  unsigned PhaseSpace::numDim() const
  {
    return spaceNumDim;
  }

  RVVectori PhaseSpace::dim() const
  {
    return RVVectori(spaceDim);
  }

  int PhaseSpace::dim(int d) const
  {
    return spaceDim(d);
  }

  bool PhaseSpace::isIntegrated(int d) const
  {
    return integrate(d);
  }

  mudfas::real PhaseSpace::getStart(int d) const
  {
    return start(d);
  }

  mudfas::real PhaseSpace::getEnd(int d) const
  {
    return end(d);
  }

  mudfas::real PhaseSpace::getStride(int d) const
  {
    return stride(d);
  }

  bool PhaseSpace::isSpeciesinSpace(const string & name) const
  {
    NameListConstIter i=speciesNames.begin(), e=speciesNames.end();
    for ( ; i!=e; ++i) {
      if (name == *i)
        return true;
    }
    return false;
  }

  void PhaseSpace::integrateParticles(real t, const Species &species, RVArrayr &space)
  {
    RVectorr rstart, rend;
    for (int d=0; d<DIMR; ++d) {
      rstart(d) = start(d) + t*speed(d);
      rend(d)   = end(d)   + t*speed(d);
    }
    for (ParticleListConstIter i=species.begin(), e=species.end(); i!=e; ++i) {
      bool inSpace = true;
      RVVectori I(-1);
      for (int d=0, j=0; d<DIMR; ++d, ++j) { // Loop over configuration space
        if ( (i->R(d) >= rstart(j) && i->R(d) < rend(j)) ||
             rstart(j) == rend(j) ) {
          if ( ! integrate(j)) {
            I(j) = int(std::floor(double(i->R(d)-rstart(j))/stride(j)));
          } else {
            I(d) = 0;
          }
        } else {
          inSpace = false;
          break;
        }
      }
      if (inSpace) {
        for (int d=0, j=DIMR; d<DIMV; ++d, ++j) { // Loop over velocity space
          if ( (i->V(d) >= start(j) && i->V(d) < end(j) ) || start(j) == end(j)) {
            if ( ! integrate(j)) {
              I(j) = int(std::floor(double(i->V(d)-start(j))/stride(j)));
            } else {
              I(j) = 0;
            }
          }	else {
            inSpace = false;
            break;
          }
        }
      }
      if (inSpace) {
        space(I) += 1;
      }
    }
  }

  void PhaseSpace::initParsing(int nargs, char *args[])
  {
    registerClass("PhaseSpace");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    string prefix(spaceName);
    prefix += "::";
    setPrefix(prefix.c_str());

    parseLevelDebugOption("dl");

    using parser::types::boolVect;
    using parser::types::stringVect;
    using parser::types::realVect;

#if 0

    const char rvvecb[] = "bool[" DIMRVSTR "]";
    const char rvvecr[] = "real[" DIMRVSTR "]";
    const char rvecr[]  = "real[" DIMRSTR "]";
#endif

    insertOption(_speciesNames, "species"  , stringVect, "List of the species names", Any(speciesNames));

    insertOption(_start       , "start"    , realVect  , "Subvolume start"          , Any(start));
    insertOption(_end         , "end"      , realVect  , "Subvolume end"            , Any(end));
    insertOption(_stride      , "stride"   , realVect  , "Subvolume stride"         , Any(stride));
    insertOption(_integrate   , "integrate", boolVect  , "Subvolume integrate flag" , Any(integrate));

    insertOption(_speed       , "speed"    , realVect  , "Subvolume speed"          , Any(speed));
  }

  void PhaseSpace::paramParsing()
  {
    parseOption(_speciesNames, speciesNames);

    parseOption(_start       , start);
    parseOption(_end         , end);
    parseOption(_stride      , stride);
    parseOption(_integrate   , integrate);
    parseOption(_speed       , speed);
  }


#undef ID

}
