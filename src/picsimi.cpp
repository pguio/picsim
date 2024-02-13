/**************************************************************************
 *
 * $Id: picsimi.cpp,v 1.73 2011/11/09 16:36:21 patrick Exp $
 *
 * Copyright (c) 2000-2011 * Patrick Guio <patrick.guio@gmail.com>
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

#include <picsimi.h>

namespace picsim {

  using std::cout;
  using std::endl;


#define ID "$Id: picsimi.cpp,v 1.73 2011/11/09 16:36:21 patrick Exp $"

  std::ostream& operator<<(std::ostream& os, const PicSimI &s)
  {
    s.printOn(os);
    return os;
  }

  PicSimI::PicSimI(int nargs, char *args[]) :
    GenericEsPicSim<mudfas::PoissonBoltzmannSolver>(nargs, args), ESE(1), KEE(2)
  {
    initParsing(nargs, args);
    paramParsing();
  }

  PicSimI::~PicSimI()
  {}


  void PicSimI::solveFields()
  {
#if defined(HAVE_MPI)
    if (isEforce()) {
      std::ios::fmtflags f = cout.flags() & std::ios::floatfield;
      int p = cout.precision();
#if defined(SOL1)
      if (rankProc == masterProc) {
        // Initialise multigrid solver
        if (! solver.isInitGuess() ) Phi = 0.0;
        solver.setPhiAndRhs(Phi, Rho);
        // Solver Poisson equation
        solver.solve();
        solver.getPhi(Phi);
#if 0
        if (solver.isInitGuess() && scheduler.getCurrentIter() == 0 ) {
          // solve again with Phi different from zero
          // NEEDS work
          solver.setPhiAndRhs(Phi, Rho);
          solver.solve();
        }
#endif
        // Calculate electrons and field energy
        solver.getElectricFieldEnergy(ESE[0]);
        solver.getFluidEnergy(KEE[0], KEE[1]);
        cout << std::setprecision(2) << std::scientific
             << "ESE     = " << ESE << '\n'
             << "KEE     = " << KEE
             << std::resetiosflags(std::ios::scientific)
             << std::setiosflags(f) << std::setprecision(p) << endl;
      }
      // Broadcast Phi to all processors
      MPI_Bcast(Phi.data(), Phi.size(), MPI_FIELD, mpiRoot, MPI_COMM_WORLD);
#else // SOL2
      // Initialise multigrid solver
      if (! solver.isInitGuess() ) Phi = 0.0;
      solver.setPhiAndRhs(Phi, Rho);
      // Solver Poisson equation
      solver.solve();
      solver.getPhi(Phi);
#if 0
      if (solver.isInitGuess() && scheduler.getCurrentIter() == 0 ) {
        // solve again with Phi different from zero
        // NEEDS work
        solver.setPhiAndRhs(Phi, Rho);
        solver.solve();
      }
#endif
      // Calculate electrons and field energy
      solver.getElectricFieldEnergy(ESE[0]);
      solver.getFluidEnergy(KEE[0], KEE[1]);

      if (rankProc == masterProc) {
        cout << std::setprecision(2) << std::scientific
             << "ESE     = " << ESE << '\n'
             << "KEE     = " << KEE
             << std::resetiosflags(std::ios::scientific)
             << std::setiosflags(f) << std::setprecision(p) << endl;
      }
#endif
      // Calculate E = -\nabla\phi with conditions BC
      Boundary<int,DIMR> BC;
      solver.getBoundaryCondition(BC);
      scheduler.setInternElectricField(Phi, BC);
    }
#else // no MPI
    if (isEforce()) {
      std::ios::fmtflags f = cout.flags() & std::ios::floatfield;
      int p = cout.precision();
      // Initialise multigrid solver
      if (! solver.isInitGuess() ) Phi = 0.0;
      solver.setPhiAndRhs(Phi, Rho);
      // Solver Poisson equation
      solver.solve();
      solver.getPhi(Phi);
#if 0
      if (solver.isInitGuess() && scheduler.getCurrentIter() == 0 ) {
        // solve again with Phi different from zero
        // NEEDS work
        solver.setPhiAndRhs(Phi, Rho);
        solver.solve();
      }
#endif
      // Calculate electrons and field energy
      solver.getElectricFieldEnergy(ESE[0]);
      solver.getFluidEnergy(KEE[0], KEE[1]);
      cout << std::setprecision(2) << std::scientific
           << "ESE     = " << ESE << '\n'
           << "KEE     = " << KEE
           << std::resetiosflags(std::ios::scientific)
           << std::setiosflags(f) << std::setprecision(p) << endl;
      // Calculate E = -\nabla\phi with conditions BC
      Boundary<int,DIMR> BC;
      solver.getBoundaryCondition(BC);
      scheduler.setInternElectricField(Phi, BC);
    }
#endif

  }

  void PicSimI::initDiagnostics()
  {
#if defined(HAVE_MPI)
    if (rankProc == masterProc) {
      diagnostics.startInit();
      diagnostics.initFields(scheduler);
      diagnostics.initEnergies(scheduler, numParticles, kineticEnergy,
                               kineticTemp);
      if (isEforce()) {
        diagnostics.initEnergies(scheduler, ESE, KEE);
      }
    }
    diagnostics.initPhaseSpaces(scheduler);
    diagnostics.initMoments(scheduler);
    if (rankProc == masterProc) {
      if (isEforce()) {
        diagnostics.initPotentialProbes(scheduler);
        diagnostics.initSpectra(scheduler);
      }
      diagnostics.endInit();
    }
#else
    diagnostics.startInit();
    diagnostics.initFields(scheduler);
    diagnostics.initEnergies(scheduler, numParticles, kineticEnergy, kineticTemp);
    if (isEforce()) {
      diagnostics.initEnergies(scheduler, ESE, KEE);
    }
    diagnostics.initPhaseSpaces(scheduler);
    diagnostics.initMoments(scheduler);
    if (isEforce()) {
      diagnostics.initPotentialProbes(scheduler);
      diagnostics.initSpectra(scheduler);
    }
    diagnostics.endInit();
#endif

  }


  void PicSimI::processDiagnostics()
  {
#if defined(HAVE_MPI)
    if (rankProc == masterProc) {
      //Field Ne(Rho.shape()); solver.getDens(Ne); diagnostics.saveFields(scheduler, Rho, Ne);
      diagnostics.saveFields(scheduler, Rho, Phi);
      diagnostics.saveEnergies(scheduler, ESE, KEE, numParticles, kineticEnergy);
    }
    diagnostics.savePhaseSpaces(scheduler, species);
    diagnostics.saveMoments(scheduler, species);
    if (rankProc == masterProc && isEforce()) {
      diagnostics.savePotentialProbes(scheduler, Phi);
      Field Rho_e;
      solver.getElectronDensity(Rho_e);
      diagnostics.saveSpectra(scheduler, Rho_e);
    }
#else
    //Field Ne(Rho.shape()); solver.getDens(Ne); diagnostics.saveFields(scheduler, Rho, Ne);
    diagnostics.saveFields(scheduler, Rho, Phi);
    diagnostics.saveEnergies(scheduler, ESE, KEE, numParticles, kineticEnergy);
    diagnostics.savePhaseSpaces(scheduler, species);
    diagnostics.saveMoments(scheduler, species);
    if (isEforce()) {
      diagnostics.savePotentialProbes(scheduler, Phi);
      Field Rho_e;
      solver.getElectronDensity(Rho_e);
      diagnostics.saveSpectra(scheduler, Rho);
    }
#endif

  }

  void PicSimI::initParsing(int nargs, char *args[])
  {
    registerClass("PicSimI");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);
  }

  void PicSimI::paramParsing()
  {}

#undef ID

}
