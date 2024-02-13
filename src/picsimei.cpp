/**************************************************************************
 *
 * $Id: picsimei.cpp,v 1.17 2011/11/09 16:36:21 patrick Exp $
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

#include <picsimei.h>

namespace picsim {

#define ID "$Id: picsimei.cpp,v 1.17 2011/11/09 16:36:21 patrick Exp $"

  std::ostream& operator<<(std::ostream& os, const PicSimEI &s)
  {
    s.printOn(os);
    return os;
  }

  PicSimEI::PicSimEI(int nargs, char *args[]) :
    GenericEsPicSim<mudfas::PoissonSolver>(nargs, args)
  {
    initParsing(nargs, args);
    paramParsing();
  }

  PicSimEI::~PicSimEI()
  {}

  void PicSimEI::solveFields()
  {
#if defined(HAVE_MPI)
    if (isEforce()) {
#if defined(SOL1)
      if (rankProc == masterProc) {
        // Initialise multigrid solver
        if (! solver.isInitGuess() ) Phi = 0.0;
        solver.setPhiAndRhs(Phi, Rho);
        // Solver Poisson equation
        solver.solve();
        solver.getPhi(Phi);
      }
#else // SOL2
      // Initialise multigrid solver
      if (! solver.isInitGuess() ) Phi = 0.0;
      solver.setPhiAndRhs(Phi, Rho);
      // Solver Poisson equation
      solver.solve();
      solver.getPhi(Phi);
      // Broadcast Phi to all processors
      MPI_Bcast(Phi.data(), Phi.size(), MPI_FIELD, mpiRoot, MPI_COMM_WORLD);
#endif
      // Calculate E = -\nabla\phi with conditions BC
      Boundary<int,DIMR> BC;
      solver.getBoundaryCondition(BC);
      scheduler.setInternElectricField(Phi, BC);
    }
#else // no NPI
    if (isEforce()) {
      // Initialise multigrid solver
      if (! solver.isInitGuess() ) Phi = 0.0;
      solver.setPhiAndRhs(Phi, Rho);
      // Solver Poisson equation
      solver.solve();
      solver.getPhi(Phi);

      // Calculate E = -\nabla\phi with conditions BC
      Boundary<int,DIMR> BC;
      solver.getBoundaryCondition(BC);
      scheduler.setInternElectricField(Phi, BC);
    }
#endif
  }

  void PicSimEI::initDiagnostics()
  {
#if defined(HAVE_MPI)
    if (rankProc == masterProc) {
      diagnostics.startInit();
      diagnostics.initFields(scheduler);
      diagnostics.initEnergies(scheduler, numParticles, kineticEnergy,
                               kineticTemp);
    }
    diagnostics.initPhaseSpaces(scheduler);
    diagnostics.initMoments(scheduler);
    if (rankProc == masterProc) {
      diagnostics.initPotentialProbes(scheduler);
      diagnostics.initSpectra(scheduler);
      diagnostics.endInit();
    }
#else
    diagnostics.startInit();
    diagnostics.initFields(scheduler);
    diagnostics.initEnergies(scheduler, numParticles, kineticEnergy, kineticTemp);
    diagnostics.initPhaseSpaces(scheduler);
    diagnostics.initMoments(scheduler);
    diagnostics.initPotentialProbes(scheduler);
    diagnostics.initSpectra(scheduler);
    diagnostics.endInit();
#endif
  }

  void PicSimEI::processDiagnostics()
  {
#if defined(HAVE_MPI)
    if (rankProc == masterProc) {
      diagnostics.saveFields(scheduler, Rho, Phi);
      diagnostics.saveEnergies(scheduler, numParticles, kineticEnergy);
    }
    diagnostics.savePhaseSpaces(scheduler, species);
    diagnostics.saveMoments(scheduler, species);
    if (rankProc == masterProc) {
      diagnostics.savePotentialProbes(scheduler, Phi);
      // calculate spectra pver Rho = electron charge density
      Field Rho_e(Rho);
      diagnostics.saveSpectra(scheduler, Rho_e);
    }
#else
    diagnostics.saveFields(scheduler, Rho, Phi);
    diagnostics.saveEnergies(scheduler, numParticles, kineticEnergy);
    diagnostics.savePhaseSpaces(scheduler, species);
    diagnostics.saveMoments(scheduler, species);
    diagnostics.savePotentialProbes(scheduler, Phi);
    Field Rho_e(Rho);
    diagnostics.saveSpectra(scheduler, Rho_e);
#endif
  }

  void PicSimEI::initParsing(int nargs, char *args[])
  {
    registerClass("PicSimEI");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);
  }

  void PicSimEI::paramParsing()
  {}


#undef ID

}
