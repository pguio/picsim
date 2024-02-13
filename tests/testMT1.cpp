/**************************************************************************
 *
 * $Id: testMT1.cpp,v 1.2 2011/03/26 15:36:09 patrick Exp $
 *
 * Copyright (c) 2008-2011 Patrick Guio <patrick.guio@gmail.com>
 * All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published
 * by the * Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
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

#include <iostream>
#include <random/uniform.h>
#if defined(HAVE_MPI)
#include <mpi.h>
#endif

using namespace ranlib;

int main(int nargs, char *args[])
{
  int seed = 1;

#if defined(HAVE_MPI)
  int nbProc, rankProc;

  MPI_Init(&nargs, &args);
  MPI_Comm_size(MPI_COMM_WORLD, &nbProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankProc);
#endif

  Uniform<double,MersenneTwister,independentState> r;

#if defined(HAVE_MPI)

  r.seed(seed + rankProc);

  for (int i=0; i<nbProc; ++i) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (i==rankProc) {
      std::cout << "Process " << rankProc <<  ":";
      for (int j=0; j<5; ++j) {
        std::cout << " " << r.random();
      }
      std::cout << std::endl;
    }
  }

  MPI_Finalize ();

#else

  r.seed(seed);
  for (int j=0; j<5; ++j) {
    std::cout << " " << r.random();
  }
  std::cout << std::endl;

#endif

  return 0;
}
