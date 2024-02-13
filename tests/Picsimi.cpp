/**************************************************************************
 *
 * $Id: Picsim.cpp.in,v 1.42 2019/05/10 16:38:06 patrick Exp $
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


#include <new>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <picsimi.h>

#if defined(HAVE_SYS_RESOURCE_H)
#include <sys/resource.h>
#endif

void my_new_handler()
{
  cerr << "Out of memory" << endl;
  abort();
}

#if defined(HAVE_MPI)


int main(int nargs, char *args[])
{
  try {
    std::set_new_handler(my_new_handler);
    cout.precision(4);

    int nbProc;
    int rankProc;
    int nameLen;
    char cpuName[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&nargs, &args);

    MPI_Comm_size(MPI_COMM_WORLD, &nbProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankProc);
    MPI_Get_processor_name(cpuName, &nameLen);


    // Program parser
    parser::Parser parser(nargs, args);
    if (rankProc == masterProc) {
      parser.registerProgram(args[0]);
      parser.registerPackage(PACKAGE, VERSION, PICSIM_COPYRIGHT);
    }

    // Simulator constructor
    picsim::PicSimI simulator(nargs, args);

#define PARSE(Fun)                                                          \
{                                                                           \
  bool parsed = false;                                                      \
  if (rankProc == masterProc) {                                             \
    parsed = parser.Fun();                                                  \
	}                                                                         \
  MPI_Bcast(static_cast<void*>(&parsed),1,MPI_BYTE,mpiRoot,MPI_COMM_WORLD); \
  if (parsed) {                                                             \
    simulator.Fun();                                                        \
		MPI_Finalize();                                                         \
    throw 0;                                                                \
  }                                                                         \
}                                                                           \
 
    PARSE(parseHelp)
    PARSE(parseVersion)
    PARSE(parseTemplate)

#undef PARSE

    simulator.initialise();

    cout << simulator << endl;

    // Output configuration
    {
      string confFile(simulator.getDiagnosticFilename());
      confFile.append(".conf");
      std::ofstream os;
      if (rankProc == masterProc)
        os.open(confFile.c_str(), ios::out|ios::trunc);
      else
        os.open(confFile.c_str(), ios::out|ios::app);
      os << simulator << endl;
    }

    // Time each process
    timeutils::Timer timer;
    timer.start();

    bool last;
    do {
      last = simulator.integrateOneStepForward();
    } while (!last);
    timer.stop();

    if (rankProc == masterProc) cout << endl;

    for (int ip=0; ip<nbProc; ip++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (ip == rankProc)
        cout << "Time (rank " << rankProc << "@" << cpuName << ") ="
             << " real " << timeutils::toHMS(timer.realElapsed())
             << ", user " << timeutils::toHMS(timer.userElapsed())
             << ", sys " << timeutils::toHMS(timer.sysElapsed())
             << ", mpi " << timeutils::toHMS(timer.mpiElapsed())
             << std::endl;
    }

#if defined(HAVE_SYS_RESOURCE_H)
    if (rankProc == masterProc) cout << endl;

    struct rusage usage;
    int status = getrusage(RUSAGE_SELF, &usage);
    if (status == 0) {
      for (int ip=0; ip<nbProc; ip++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (ip == rankProc)
          cout << "Ruse (rank " << rankProc << "@" << cpuName << ") ="
               << " utime " << usage.ru_utime.tv_sec << " s"
               << ", stime " << usage.ru_stime.tv_sec << " s"
               << ", maxrss " << usage.ru_maxrss/1024 << " Mb"
               << ", minflt " << usage.ru_minflt
               << ", majflt " << usage.ru_majflt
               << ", inblock " << usage.ru_inblock
               << ", oublock " << usage.ru_oublock
               << ", nvcsw " << usage.ru_nvcsw
               << ", nivcsw " << usage.ru_nivcsw
               << std::endl;
      }
		}
#endif

    MPI_Finalize();

    return 0;

  } catch(const ClassException& c) {
      cerr << c.what() << endl;
    return !0;
  } catch(const std::exception& e) {
      cerr << e.what() << endl;
    return !0;
  } catch(const int status) {
      cerr << "status=" << status << endl;
    return status;
  }
}

#else // no MPI support


int main(int nargs, char *args[])
{
  try {
    std::set_new_handler(my_new_handler);
    cout.precision(4);

    parser::Parser parser(nargs, args);
    parser.registerProgram(args[0]);
    parser.registerPackage(PACKAGE, VERSION, PICSIM_COPYRIGHT);

    picsim::PicSimI simulator(nargs, args);

#define PARSE(Fun)                 \
if (parser.Fun()) {                \
	simulator.Fun();                 \
	return 0;                        \
}                                  \
 
    PARSE(parseHelp)
    PARSE(parseVersion)
    PARSE(parseTemplate)

#undef PARSE

    simulator.initialise();

    cout << simulator << endl;

    {
      string confFile(simulator.getDiagnosticFilename());
      confFile.append(".conf");
      std::ofstream os(confFile.c_str(), ios::out|ios::trunc);
      os << simulator << endl;
    }


    // Time each process
    timeutils::Timer timer;
    timer.start();

    bool last;
    do {
      last = simulator.integrateOneStepForward();
    } while (!last);
    timer.stop();

    cout << "\nTime ="
         << " real " << timeutils::toHMS(timer.realElapsed())
         << ", user " << timeutils::toHMS(timer.userElapsed())
         << ", sys " << timeutils::toHMS(timer.sysElapsed())
         << endl;

#if defined(HAVE_SYS_RESOURCE_H)
    struct rusage usage;
    int status = getrusage(RUSAGE_SELF, &usage);
		if (status == 0) {
      cout << "\nRuse ="
           << " utime " << usage.ru_utime.tv_sec << " s"
           << ", stime " << usage.ru_stime.tv_sec << " s"
           << ", maxrss " << usage.ru_maxrss/1024 << " Mb"
           << ", minflt " << usage.ru_minflt
           << ", majflt " << usage.ru_majflt
           << ", inblock " << usage.ru_inblock
           << ", oublock " << usage.ru_oublock
           << ", nvcsw " << usage.ru_nvcsw
           << ", nivcsw " << usage.ru_nivcsw
           << std::endl;
    } 

#endif

    return 0;

  } catch(const ClassException& c) {
      cerr << c.what() << endl;
    return !0;
  } catch(const std::exception& e) {
      cerr << e.what() << endl;
    return !0;
  } catch(const int status) {
      cerr << "status=" << status << endl;
    return status;
  }
}

#endif
