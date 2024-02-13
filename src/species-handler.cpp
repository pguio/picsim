/**************************************************************************
 *
 * $Id: species-handler.cpp,v 1.52 2024/01/26 17:47:20 patrick Exp $
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

#include <species-handler.h>

namespace picsim {

  using blitz::mean;
  using blitz::min;
  using blitz::max;

  using parser::header;

  using std::cout;
  using std::endl;
  using std::ios;
  using std::ostream;

#define ID "$Id: species-handler.cpp,v 1.52 2024/01/26 17:47:20 patrick Exp $"

  ostream& operator<<(ostream& os, const SpeciesHandler &s)
  {
#if defined(HAVE_MPI)
    if (s.rankProc == masterProc)
#endif
    {
      os << header("Species Handler setup");
      SpeciesHandler::SpeciesTypes::const_iterator i = s.types.begin();
      SpeciesHandler::SpeciesTypes::const_iterator e = s.types.end();

#if defined(HAVE_MPI)
      SpeciesHandler::SpeciesNumCpu::const_iterator cpus = s.numCpu.begin();
      for (; i!= e; ++i, ++cpus) {
        if (! i->second.empty()) {
          os << i->first << " species        = " << i->second << '\n'
             << "Number of CPUs            = " << cpus->second << endl;
        }
      }
#else
      for (; i!= e; ++i) {
        if (! i->second.empty()) {
          os << i->first << " species        = " << i->second << endl;
        }
      }
#endif
    }


#if defined(HAVE_MPI)

#if 0
    SpeciesHandler::SpeciesListConstIter  i=s.begin();
    for (int ip=0; ip<s.nbProc; ++ip) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (ip == s.rankProc) {
        os << "\n## species " << (**i).name() << '\n' << **i << endl;
      }
    }
#else
    MPI_Barrier(MPI_COMM_WORLD);
    int recv, sent=s.rankProc;
    SpeciesHandler::SpeciesListConstIter  i=s.begin();
    MPI_Status status;
    if (s.rankProc > 0)
      MPI_Recv(&recv, 1, MPI_INT,
               s.rankProc-1, 100, MPI_COMM_WORLD, &status);
    os << "\n## species " << (**i).name() << '\n' << **i << endl;
    if (s.rankProc < s.nbProc-1)
      MPI_Send(&sent, 1, MPI_INT,
               s.rankProc+1, 100, MPI_COMM_WORLD);
#endif

#else

    SpeciesHandler::SpeciesListConstIter i=s.begin(), e=s.end();
    for ( ; i!=e; ++i) {
      os << "\n## species " << (**i).name() << '\n' << **i << '\n';
    }

#endif
    return os;
  }


  SpeciesHandler::SpeciesHandler(int nargs, char *args[]) : Parser(nargs, args),
    speciesStartKey(0), numSpecies(0)
  {
    initParsing(nargs, args);
    paramParsing();

    Factory::instance().init();

#if defined(HAVE_MPI)

    SpeciesTypes::const_iterator i = types.begin();
    SpeciesTypes::const_iterator e = types.end();
    SpeciesNumCpu::const_iterator n = numCpu.begin();
    for (; i!= e; ++i, ++n) { // For each species types
      NameList::const_iterator name = i->second.begin();
      NameList::const_iterator last = i->second.end();
      VecInt::const_iterator ncpu = n->second.begin();
      for (; name != last; ++name, ++ncpu) { // For each object of species types
        for (int j=0; j<*ncpu; ++j) { // For *ncpu processors
          if (rankProc == numSpecies) {
            species.insert(species.end(),
                           Factory::instance()
                           .create(i->first, nargs, args, *name));
          }
          ++numSpecies;
        }
      }
    }

#else

    SpeciesTypes::const_iterator i = types.begin();
    SpeciesTypes::const_iterator e = types.end();
    for (; i!= e; ++i) {
      NameList::const_iterator name = i->second.begin();
      NameList::const_iterator last = i->second.end();
      for (; name != last; ++name) {
        species.insert(species.end(),
                       Factory::instance()
                       .create(i->first, nargs, args, *name));
        ++numSpecies;
      }
    }

#endif

    checkParam();
  }

  SpeciesHandler::~SpeciesHandler()
  {
#if defined(HAVE_MPI)

    SpeciesListIter i=species.begin();

    BEGIN_DEBUG_OUTPUT(30)
    HEADER_DEBUG_OUTPUT1("SpeciesHandler::~SpeciesHandler")
    VAR_DEBUG_OUTPUT1((**i).name())
    END_DEBUG_OUTPUT

    delete *i;

#else

    for (SpeciesListIter i=species.begin(), e=species.end(); i != e; ++i) {

      BEGIN_DEBUG_OUTPUT(30)
      HEADER_DEBUG_OUTPUT1("SpeciesHandler::~SpeciesHandler")
      VAR_DEBUG_OUTPUT1((**i).name())
      END_DEBUG_OUTPUT

      delete *i;
    }

#endif

  }


#if defined(HAVE_MPI)

#define PARSE(Fun)                                                  \
bool SpeciesHandler::Fun() const                                    \
{                                                                   \
	bool parsed = false;                                              \
	if (rankProc == masterProc) {                                     \
		parsed = Parser::Fun();                                         \
	}                                                                 \
	MPI_Bcast(static_cast<void *>(&parsed), 1, MPI_BYTE,              \
			      mpiRoot, MPI_COMM_WORLD);                               \
	if (parsed) {                                                     \
		for (int ip=0; ip<nbProc; ++ip) {                               \
			if (ip == rankProc) {                                         \
				SpeciesListConstIter  i = species.begin();                  \
				(**i).Fun();                                                \
			}                                                             \
			MPI_Barrier(MPI_COMM_WORLD);                                  \
		}                                                               \
		return true;                                                    \
	}                                                                 \
	return false;                                                     \
}

#else

#define PARSE(Fun)                                                  \
bool SpeciesHandler::Fun() const                                    \
{                                                                   \
	if ( Parser::Fun() ) {                                            \
		SpeciesListConstIter i=species.begin(), e=species.end();        \
		for ( ; i!=e; ++i) {                                            \
			(**i).Fun();                                                  \
		}                                                               \
		return true;                                                    \
	}                                                                 \
	return false;                                                     \
}

#endif

  PARSE(parseHelp)
  PARSE(parseVersion)
  PARSE(parseTemplate)

#undef PARSE


  void SpeciesHandler::initialise(NumberGenerator &draw, Scheduler &scheduler)
  {
    if (species.empty()) {
      Factory::instance().content();
      throw ClassException("SpeciesHandler",
                           "You must specify at least one species from the factory content.\n"
                           "Try to run command with the --help option.");
    }

    rhoFactor = 1.0/scheduler.getRefLambda()/scheduler.getUnitVolume();

#if defined(HAVE_MPI)

    SpeciesListIter i=species.begin();
    (**i).initialise(draw, scheduler);
#else

    SpeciesListIter i=species.begin(), e=species.end();
    for ( ; i!=e; ++i) {
      (**i).initialise(draw, scheduler);
    }
#endif

  }

#if defined(HAVE_MPI)

  void SpeciesHandler::depositCharge(Scheduler &scheduler, Field &Rho)
  {
    Field Rhos(Rho.shape());
    SpeciesListConstIter i = species.begin();
    (**i).depositCharge(scheduler, Rhos);

#if defined(SOL1)
    MPI_Reduce(Rhos.data(), Rho.data(), Rhos.size(), MPI_FIELD, MPI_SUM,
               mpiRoot, MPI_COMM_WORLD);
#else // SOL2
    MPI_Allreduce(Rhos.data(), Rho.data(), Rhos.size(), MPI_FIELD, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT1("SpeciesHandler::depositCharge")
    VAR_DEBUG_OUTPUT1(mean(Rhos))
    END_DEBUG_OUTPUT

    BEGIN_DEBUG_MASTER_OUTPUT(10)
    HEADER_DEBUG_OUTPUT1("SpeciesHandler::depositCharge")
    VAR_DEBUG_OUTPUT2(mean(Rho),mean(Rho)*rhoFactor)
    VAR_DEBUG_OUTPUT2(min(Rho),max(Rho))
    END_DEBUG_MASTER_OUTPUT

#if defined(SOL1)
    if (rankProc == masterProc) {
      Rho *= rhoFactor;
    }
#else // SOL2
    Rho *= rhoFactor;
#endif
  }

#else // !defined(HAVE_MPI)

  void SpeciesHandler::depositCharge(Scheduler &scheduler, Field &Rho)
  {
    Rho = 0.0;
    Field Rhos(Rho.shape());
    for (SpeciesListConstIter i=species.begin(), e=species.end(); i!=e; ++i) {
      (**i).depositCharge(scheduler, Rhos);
      Rho += Rhos;

      BEGIN_DEBUG_OUTPUT(10)
      HEADER_DEBUG_OUTPUT1("SpeciesHandler::depositCharge")
      VAR_DEBUG_OUTPUT1(mean(Rhos))
      END_DEBUG_OUTPUT

    }

    BEGIN_DEBUG_MASTER_OUTPUT(10)
    HEADER_DEBUG_OUTPUT1("SpeciesHandler::depositCharge")
    VAR_DEBUG_OUTPUT2(mean(Rho),mean(Rho)*rhoFactor)
    END_DEBUG_MASTER_OUTPUT

    Rho *= rhoFactor;
  }

#endif

  void SpeciesHandler::moveParticles(NumberGenerator &draw, Scheduler &scheduler)
  {
#if defined(HAVE_MPI)
    SpeciesListIter i = species.begin();
    (**i).moveParticles(draw, scheduler);
#else

    for (SpeciesListIter i=species.begin(), e=species.end(); i!=e; ++i) {
      (**i).moveParticles(draw, scheduler);
    }
#endif

  }

#if defined(HAVE_MPI)

  void SpeciesHandler::getState(VecReal &numParticles, VecReal &kineticEnergy,
                                VecReal &kineticTemp)
  {
    SpeciesListConstIter i = species.begin();
    real send[2];
    send[0] = (**i).getKineticEnergy();
    send[1] = (**i).getParticleNumber();
    real * recv = new real [2*numSpecies];
    // Carlo Cavazzoni @ CINECA suggested to try
    // MPI_Reduce instead of MPI_Gather
    // as it is most likely best optimised routine
    MPI_Gather(send, 2, MPI_COORDINATE, recv, 2, MPI_COORDINATE, mpiRoot, MPI_COMM_WORLD);
    if (rankProc == masterProc) {
      for (int is=0; is<numSpecies; ++is) {
        kineticEnergy[is] = recv[2*is];
        numParticles[is]  = recv[2*is+1];
        kineticTemp[is]   = (numParticles[is] != 0 ?
                             kineticEnergy[is]/numParticles[is] :
                             0.0);
      }
      std::ios::fmtflags f = cout.flags() & std::ios::floatfield;
      int p = cout.precision();
      cout << std::setprecision(0) << std::fixed
           << "numParticles      = " << numParticles << '\n'
           << std::resetiosflags(std::ios::fixed)
           << std::setprecision(2) << std::scientific
           << "kineticEnergy     = " << kineticEnergy << '\n'
           << "kineticTemp       = " << kineticTemp
           << std::resetiosflags(std::ios::scientific)
           << std::setiosflags(f) << std::setprecision(p) << endl;
    }
    delete[] recv;
  }

#else // !defined(HAVE_MPI)

  void SpeciesHandler::getState(VecReal &numParticles, VecReal &kineticEnergy,
                                VecReal &kineticTemp)
  {
    int is = 0;
    for (SpeciesListConstIter i=species.begin(), e=species.end(); i!=e; ++i) {
      kineticEnergy[is] = (**i).getKineticEnergy();
      numParticles[is]  = (**i).getParticleNumber();
      kineticTemp[is]   = (numParticles[is] != 0 ?
                           kineticEnergy[is]/numParticles[is] :
                           0.0);
      ++is;
    }
    std::ios::fmtflags f = cout.flags() & std::ios::floatfield;
    int p = cout.precision();
    cout << std::setprecision(0) << std::fixed
         << "numParticles      = " << numParticles << '\n'
         << std::resetiosflags(std::ios::fixed)
         << std::setprecision(2) << std::scientific
         << "kineticEnergy     = " << kineticEnergy << '\n'
         << "kineticTemp       = " << kineticTemp
         << std::resetiosflags(std::ios::scientific)
         << std::setiosflags(f) << std::setprecision(p) << endl;
  }

#endif

#if defined(HAVE_MPI)

  Boundary<real,DIMR> SpeciesHandler::getAlfa() const
  {
    Boundary<mudfas::real,DIMR> alfas(0.0);
    Boundary<mudfas::real,DIMR> alfa;
    MPI_Reduce(alfas.min.data(), alfa.min.data(), DIMR, MPI_COORDINATE,
               MPI_SUM, mpiRoot, MPI_COMM_WORLD);
    MPI_Reduce(alfas.max.data(), alfa.max.data(), DIMR, MPI_COORDINATE,
               MPI_SUM, mpiRoot, MPI_COMM_WORLD);
    return alfa;
  }

#else // !defined(HAVE_MPI)

  Boundary<real,DIMR> SpeciesHandler::getAlfa() const
  {
    Boundary<mudfas::real,DIMR> alfa(0.0);
    for (SpeciesListConstIter i=species.begin(), e=species.end(); i!=e; ++i) {
      alfa += (**i).getAlfa();
    }
    return alfa;
  }


#endif

  int SpeciesHandler::getSpeciesNumber() const
  {
    return numSpecies;
  }

  SpeciesHandler::SpeciesListIter      SpeciesHandler::begin()
  {
    return species.begin();
  }

  SpeciesHandler::SpeciesListIter      SpeciesHandler::end()
  {
    return species.end();
  }

  SpeciesHandler::SpeciesListConstIter SpeciesHandler::begin() const
  {
    return species.begin();
  }

  SpeciesHandler::SpeciesListConstIter SpeciesHandler::end() const
  {
    return species.end();
  }

  int SpeciesHandler::size() const
  {
    return species.size();
  }

  void SpeciesHandler::initParsing(int nargs, char *args[])
  {
    registerClass("SpeciesHandler");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    parseLevelDebugOption("SpeciesHandler::dl");

    using parser::types::stringVect;
#if defined(HAVE_MPI)

    using parser::types::intVect;
#endif

    Factory::Species_iterator i = Factory::instance().begin();
    Factory::Species_iterator e = Factory::instance().end();
    int speciesKey = speciesStartKey;
    for (; i!= e ; ++i) {
      const IDKeyType speciesName(Factory::instance().getKey(i));
      types[speciesName] = NameList(0);
      const string speciesDesc("Vector of "+speciesName+" species names");
      insertOption(speciesKey++, speciesName, stringVect, speciesDesc, Any(types[speciesName]));

#if defined(HAVE_MPI)

      const string cpuOpt(speciesName+"::ncpu");
      const string cpuOptDesc("Vector of number of CPU/"+speciesName);
      insertOption(speciesKey++, cpuOpt, intVect, cpuOptDesc, Any(numCpu[speciesName]));
#endif

    }
  }

  void SpeciesHandler::paramParsing()
  {

    SpeciesTypes::const_iterator i = types.begin();
    SpeciesTypes::const_iterator e = types.end();
    int speciesKey = speciesStartKey;
    for (; i!= e ; ++i) {
      parseOption(speciesKey++, types[i->first]);

#if defined(HAVE_MPI)

      parseOption(speciesKey++, numCpu[i->first]);
#endif

    }
  }

  void SpeciesHandler::checkParam() const
  {
#if defined(HAVE_MPI)
#if 0
    if (numSpecies == 0) {
      if (rankProc == masterProc) {
        throw ClassException("SpeciesHandler",
                             "You must specify at least one species");
      } else {
        throw(EXIT_FAILURE);
      }
    }
#endif
    if (numSpecies != nbProc) {
      if (rankProc == masterProc) {
        ostringstream os;
        os << "You must have as many processes as species:\n\t" <<
           "numSpecies=" << numSpecies << " nbProc=" << nbProc;
        throw ClassException("SpeciesHandler", os.str());
      } else {
        throw(EXIT_FAILURE);
      }
    }
#else
#if 0
    if (numSpecies == 0) {
      throw ClassException("SpeciesHandler",
                           "You must specify at least one species");
    }
#endif
#endif

  }

#undef ID

}
