/**************************************************************************
 *
 * $Id: scheduler.cpp,v 1.99 2012/02/25 09:25:01 patrick Exp $
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

#include <scheduler.h>

namespace picsim {

  using blitz::product;

  using parser::header;
  using parser::map_elt;
  using parser::yesno;

  using std::cout;
  using std::endl;
  using std::ios;
  using std::ostream;
  using std::setiosflags;


#define ID "$Id: scheduler.cpp,v 1.99 2012/02/25 09:25:01 patrick Exp $"

  ostream& operator<<(ostream& os, const Scheduler &s)
  {

    os << header("Dynamic scheduler setup");
    os
        << "timeinc         = " << s.TimeInc << '\n'
        << "maxiter         = " << s.MaxIter << '\n'

        << "refmass         = " << s.RefMass << '\n'
        << "reflambda       = " << s.RefLambda << '\n'
        << "reftemp         = " << s.RefTemperature << '\n';

    os << "electric force  = " << yesno(s.Eforce) << '\n';
    if (s.Eforce) {
      os << "emode           = " << map_elt(s.emode_map, s.EMode) << '\n';
      os << "E               = " << s.Eext << '\n';
    }

    os << "magnetic force  = " << yesno(s.Bforce) << '\n';
    if (s.Bforce) {
      os << "B               = " << s.Bext << '\n';
    }

    os << "gravity force   = " << yesno(s.Gforce) << '\n';
    if (s.Gforce) {
      os << "G               = " << s.Gext << '\n';
    }
    return os;
  }


  Scheduler::Scheduler(int nargs, char *args[]) : parser::Parser(nargs, args),
    TimeInc(0.2), CurrentTime(0.0), MaxIter(50), CurrentIter(0),
    RefMass(1), RefLambda(20), RefTemperature(1),
    GridSize(mudfas::DEFAULT_GRIDSIZE), Domain(0.0, mudfas::DEFAULT_GRIDSIZE-1),
    Eforce(false), Eext(DEFAULT_E), EMode(momentum_conserving),
    Bforce(false), Bext(DEFAULT_B), Gforce(false), Gext(DEFAULT_G)
  {
    emode_map.insert(Pair(momentum_conserving,"Momentum conserving"));
    emode_map.insert(Pair(energy_conserving,  "  Energy conserving"));

    initParsing(nargs, args);
    paramParsing();

    checkParam();
  }

  Scheduler::~Scheduler()
  {}



  void Scheduler::initialise()
  {
    RVectorr dr(Domain.max-Domain.min);
    Volume = product(dr);
    dr = dr/(GridSize-1);
    UnitVolume = product(dr);

    if (Eforce) {
      switch (EMode) {
      case momentum_conserving:
        for (int d=0; d<DIMR; ++d)
          Ei(d).resize(GridSize);
        break;
      case energy_conserving:
        for (int d=0; d<DIMR; ++d)
          Ei(d).resize(GridSize-1);
        break;
      }
    } else {
      Eext = 0.0;
    }


    if (!Bforce) {
      Bext = 0.0;
    }
    if (!Gforce) {
      Gext = 0.0;
    }

    ForceSetUp = (Eforce ? 1 : 0) * 1 +
                 (Bforce ? 1 : 0) * 2 +
                 (Gforce ? 1 : 0) * 4;

    BEGIN_DEBUG_MASTER_OUTPUT(20)
    HEADER_DEBUG_OUTPUT1("Scheduler::initialise")
    VAR_DEBUG_OUTPUT1(ForceSetUp)
    END_DEBUG_MASTER_OUTPUT

    timer.start();
  }


  void Scheduler::update()
  {
    CurrentTime += TimeInc;
    ++CurrentIter;

#if defined(HAVE_MPI)

    if (rankProc == masterProc) {
#endif
      double usedTime = timer.realRunningElapsed();
      double remTime = usedTime*(MaxIter-CurrentIter)/CurrentIter;

      using timeutils::iterStatus;
      using timeutils::toHMS;

      if (CurrentIter <= MaxIter) {
        std::ios::fmtflags f = cout.flags() & std::ios::floatfield;
        int p = cout.precision();
        cout << "Iteration " << iterStatus(CurrentIter, MaxIter)
             << " tau=" << std::setw(4) << std::setprecision(3) << std::fixed
             << CurrentTime << std::resetiosflags(std::ios::fixed)
             << std::setiosflags(f) << std::setprecision(p)
             << " (used time " << toHMS(usedTime)
             << ", remaining=" << toHMS(remTime) << ')' << endl;
      }
#if defined(HAVE_MPI)

    }
#endif

  }


#if (DIMR==2)
  void Scheduler::setInternElectricField(Field &phi, Boundary<int,DIMR> &BC)
  {
    Range all(Range::all());
    Range I, J;
#if !defined(_AIX)

    RVectorr h2(2.0*(Domain.max-Domain.min)/(GridSize-1));
#else

    RVectorr h2((Domain.max-Domain.min)/(GridSize-1));
    h2 *= 2.0;
#endif

    switch (EMode) {
    case momentum_conserving:
      I.setRange(1,GridSize(0)-2);
      Ei(0)(I,all) = -(phi(I+1,all)-phi(I-1,all))/h2(0);

      J.setRange(1,GridSize(1)-2);
      Ei(1)(all,J) = -(phi(all,J+1)-phi(all,J-1))/h2(1);

      fixBoundaryInternElectricField(phi, BC);
      break;
    case energy_conserving:
      I.setRange(0,GridSize(0)-2);
      J.setRange(0,GridSize(1)-2);

      Ei(0) = -(phi(I+1,J)+phi(I+1,J+1)-phi(I,J)-phi(I,J+1))/h2(0);
      Ei(1) = -(phi(I,J+1)+phi(I+1,J+1)-phi(I,J)-phi(I+1,J))/h2(1);
      break;
    }
  }
#elif (DIMR==3)
  void Scheduler::setInternElectricField(Field &phi, Boundary<int,DIMR> &BC)
  {
    Range all(Range::all());
    Range I, J , K;
#if !defined(_AIX)

    RVectorr h2(2.0*(Domain.max-Domain.min)/(GridSize-1));
    RVectorr h4(2.0*h2);
#else

    RVectorr h2((Domain.max-Domain.min)/(GridSize-1));
    h2 *= 2.0;
    RVectorr h4(h2);
    h4 *= 2.0;
#endif

    switch (EMode) {
    case momentum_conserving:
      I.setRange(1,GridSize(0)-2);
      Ei(0)(I,all,all) = -(phi(I+1,all,all)-phi(I-1,all,all))/h2(0);

      J.setRange(1,GridSize(1)-2);
      Ei(1)(all,J,all) = -(phi(all,J+1,all)-phi(all,J-1,all))/h2(1);

      K.setRange(1,GridSize(2)-2);
      Ei(2)(all,all,K) = -(phi(all,all,K+1)-phi(all,all,K-1))/h2(2);

      fixBoundaryInternElectricField(phi,BC);
      break;
    case energy_conserving:
      I.setRange(0,GridSize(0)-2);
      J.setRange(0,GridSize(1)-2);
      K.setRange(0,GridSize(2)-2);

      Ei(0) = -(phi(I+1,J,K)+phi(I+1,J+1,K)+phi(I+1,J+1,K+1)+phi(I+1,J,K+1)-
                phi(I,J,K)-phi(I,J+1,K)-phi(I,J+1,K+1)-phi(I,J,K+1))/h4(0);

      Ei(1) = -(phi(I,J+1,K)+phi(I+1,J+1,K)+phi(I+1,J+1,K+1)+phi(I,J+1,K+1)-
                phi(I,J,K)-phi(I+1,J,K)-phi(I+1,J,K+1)-phi(I,J,K+1))/h4(1);

      Ei(2) = -(phi(I,J,K+1)+phi(I+1,J,K+1)+phi(I+1,J+1,K+1)+phi(I,J+1,K+1)-
                phi(I,J,K)-phi(I+1,J,K)-phi(I+1,J+1,K)-phi(I,J+1,K))/h4(2);
      break;
    }
  }
#endif


#if (DIMR==2)
  void Scheduler::fixBoundaryInternElectricField(
    Field &phi, Boundary<int,DIMR> &BC)
  {
    using mudfas::MultiGridSolver;

    Range all(Range::all());
    RVectorr h((Domain.max-Domain.min)/(GridSize-1));
#if !defined(_AIX)

    RVectorr h2(2.0*h);
#else

    RVectorr h2(h*2.0);
#endif

    real alfamin(0.0);
    real alfamax(0.0);
    for (int d=0; d<DIMR; ++d) {
      int end = GridSize(d)-1;
      switch (BC.min(d)) {
      case MultiGridSolver::periodic:
        switch (d) {
        case 0:
          Ei(d)(0,all) = -(phi(1,all)-phi(end-1,all))/h2(d);
          Ei(d)(end,all) = Ei(d)(0,all);
          break;
        case 1:
          Ei(d)(all,0) = -(phi(all,1)-phi(all,end-1))/h2(d);
          Ei(d)(all,end) = Ei(d)(all,0);
          break;
        }
        break;
      case MultiGridSolver::dirichlet:
        switch (d) {
        case 0:
          Ei(d)(0,all) = -(phi(1,all)-phi(0,all))/h(d);
          break;
        case 1:
          Ei(d)(all,0) = -(phi(all,1)-phi(all,0))/h(d);
          break;
        }
        break;
      case MultiGridSolver::mixed:
        switch (d) {
        case 0:
          Ei(d)(0,all) = -alfamin;
          break;
        case 1:
          Ei(d)(all,0) = -alfamin;
          break;
        }
        break;
      }
      switch (BC.max(d)) {
      case MultiGridSolver::periodic:
        break;
      case MultiGridSolver::dirichlet:
        switch (d) {
        case 0:
          Ei(d)(end,all) = -(phi(end,all)-phi(end-1,all))/h(d);
          break;
        case 1:
          Ei(d)(all,end) = -(phi(all,end)-phi(all,end-1))/h(d);
          break;
        }
        break;
      case MultiGridSolver::mixed:
        switch (d) {
        case 0:
          Ei(d)(end,all) = alfamax;
          break;
        case 1:
          Ei(d)(all,end) = alfamax;
          break;
        }
        break;
      }
    }
  }
#elif (DIMR==3)
  void Scheduler::fixBoundaryInternElectricField(
    Field &phi, Boundary<int,DIMR> &BC)
  {
    using mudfas::MultiGridSolver;

    Range all(Range::all());
    RVectorr h((Domain.max-Domain.min)/(GridSize-1));
#if !defined(_AIX)

    RVectorr h2(2.0*h);
#else

    RVectorr h2(h*2.0);
#endif

    real alfamin(0.0);
    real alfamax(0.0);
    for (int d=0; d<DIMR; ++d) {
      int end=GridSize(d)-1;
      switch (BC.min(d)) {
      case MultiGridSolver::periodic:
        switch (d) {
        case 0:
          Ei(d)(0,all,all) = -(phi(1,all,all)-phi(end-1,all,all))/h2(d);
          Ei(d)(end,all,all) = Ei(d)(0,all,all);
          break;
        case 1:
          Ei(d)(all,0,all) = -(phi(all,1,all)-phi(all,end-1,all))/h2(d);
          Ei(d)(all,end,all) = Ei(d)(all,0,all);
          break;
        case 2:
          Ei(d)(all,all,0) = -(phi(all,all,1)-phi(all,all,end-1))/h2(d);
          Ei(d)(all,all,end) = Ei(d)(all,all,0);
          break;
        }
        break;
      case MultiGridSolver::dirichlet:
        switch (d) {
        case 0:
          Ei(d)(0,all,all) = -(phi(1,all,all)-phi(0,all,all))/h(d);
          break;
        case 1:
          Ei(d)(all,0,all) = -(phi(all,1,all)-phi(all,0,all))/h(d);
          break;
        case 2:
          Ei(d)(all,all,0) = -(phi(all,all,1)-phi(all,all,0))/h(d);
          break;
        }
        break;
      case MultiGridSolver::mixed:
        switch (d) {
        case 0:
          Ei(d)(0,all,all) = -alfamin;
          break;
        case 1:
          Ei(d)(all,0,all) = -alfamin;
          break;
        case 2:
          Ei(d)(all,all,0) = -alfamin;
          break;
        }
        break;
      }
      switch (BC.max(d)) {
      case MultiGridSolver::periodic:
        break;
      case MultiGridSolver::dirichlet:
        switch (d) {
        case 0:
          Ei(d)(end,all,all) = -(phi(end,all,all)-phi(end-1,all,all))/h(d);
          break;
        case 1:
          Ei(d)(all,end,all) = -(phi(all,end,all)-phi(all,end-1,all))/h(d);
          break;
        case 2:
          Ei(d)(all,all,end) = -(phi(all,all,end)-phi(all,all,end-1))/h(d);
          break;
        }
        break;
      case MultiGridSolver::mixed:
        switch (d) {
        case 0:
          Ei(d)(end,all,all) = alfamax;
          break;
        case 1:
          Ei(d)(all,end,all) = alfamax;
          break;
        case 2:
          Ei(d)(all,all,end) = alfamax;
          break;
        }
        break;
      }
    }
  }
#endif

  void Scheduler::initParsing(int nargs, char *args[])
  {
    registerClass("Scheduler");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    parseLevelDebugOption("Scheduler::dl");

    using parser::types::boolean;
    using parser::types::integer;
    using parser::types::real;
    using parser::types::intVect;
    using parser::types::realVect;

#if 0
    const char rveci[] = "int["  DIMRSTR "]";
    const char rvecr[] = "real[" DIMRSTR "]";
#endif

    insertOption(timeinc  , "timeinc"  , real    , "Time increment"                  , Any(TimeInc));
    insertOption(maxiter  , "maxiter"  , integer , "Maximum iterations"              , Any(MaxIter));

    insertOption(refmass  , "refmass"  , real    , "Reference mass"                  , Any(RefMass));
    insertOption(reflambda, "reflambda", real    , "Reference lambda"                , Any(RefLambda));
    insertOption(reftemp  , "reftemp"  , real    , "Reference temperature"           , Any(RefTemperature));

    insertOption(gridsize , "gridsize" , intVect , "Grid size"                       , Any(GridSize));
    insertOption(rmin     , "rmin"     , realVect, "Value of lower boundary"         , Any(Domain.min));
    insertOption(rmax     , "rmax"     , realVect, "Value of upper boundary"         , Any(Domain.max));

    insertOption(eforce   , "eforce"   , boolean , "Enable/disable electric force"   , Any(Eforce));
    insertOption(emode    , "emode"    , integer , "Electric field calculation mode" , Any(EMode));
    insertOption(evector  , "evector"  , realVect, "Electric field"                  , Any(Eext));

    insertOption(bforce   , "bforce"   , boolean , "Enable/disable magnetic force"   , Any(Bforce));
    insertOption(bvector  , "bvector"  , realVect, "Magnetic field"                  , Any(Bext));

    insertOption(gforce   , "gforce"   , boolean , "Enable/disable gravitation force", Any(Gforce));
    insertOption(gvector  , "gvector"  , realVect, "Gravitation field"               , Any(Gext));
  }

  void Scheduler::paramParsing()
  {
    parseOption(timeinc  , TimeInc);
    parseOption(maxiter  , MaxIter);

    parseOption(refmass  , RefMass);
    parseOption(reflambda, RefLambda);
    parseOption(reftemp  , RefTemperature);

    parseOption(gridsize , GridSize);
    parseOption(rmin     , Domain.min);
    parseOption(rmax     , Domain.max);

    parseOption(eforce   , Eforce);
    parseOption(emode    , EMode);
    parseOption(evector  , Eext);

    parseOption(bforce   , Bforce);
    parseOption(bvector  , Bext);

    parseOption(gforce   , Gforce);
    parseOption(gvector  , Gext);
  }

  void Scheduler::checkParam() const
  {
    checkMap(emode, emode_map, EMode);
  }


#undef ID

}
