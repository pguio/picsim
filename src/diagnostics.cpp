/************************************************************************** *
 * $Id: diagnostics.cpp,v 1.160 2011/03/26 15:36:08 patrick Exp $
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

#include <diagnostics.h>

namespace picsim {

  using parser::header;
  using parser::yesno;

  using std::ostream;
  using std::string;

#define ID "$Id: diagnostics.cpp,v 1.160 2011/03/26 15:36:08 patrick Exp $"

  const string numParticleHdfName("nbp");
  const string densityHdfName("density");
  const string potentialHdfName("potential");

  const string keiHdfName("kei");
  const string keeHdfName("kee");
  const string eseHdfName("ese");

  const string speciesHdfName("species");

  const string timeFieldHdfName("time-field");
  const string timeVectorHdfName("time-vector");

#if (DIMR==2)

  const char *rvAxis[]= { "x", "y", "vx", "vy", "vz"
                        };
#elif (DIMR==3)

  const char *rvAxis[]= { "x", "y", "z", "vx", "vy", "vz"
                        };
#endif

  const char *momsAxis[]= { "vx", "vy", "vz", "Tx", "Ty", "Tz"
                          };

  ostream& operator<<(ostream& os, const Diagnostics &d)
  {
    return os << header("Diagnostics setup")

           << d.hdf << '\n'

           << "Iterations          = " << d.iSpec << '\n'
           << "Density             = " << yesno(d.densityDiag) << '\n'
           << "Potential           = " << yesno(d.potentialDiag) << '\n'
           << "Energy              = " << yesno(d.energyDiag) << '\n'

           << d.phaseSpaces << '\n'
           << d.moments << '\n'
           << d.phiProbes << '\n'
           << d.spectra << '\n';
  }

  Diagnostics::Diagnostics(int nargs, char *args[]) : Parser(nargs, args),
    hdf(nargs, args), iSpec(DEFAULT_SaveIter),
    densityDiag(false), potentialDiag(false), energyDiag(false),
    phaseSpaces(nargs, args), moments(nargs, args),
    phiProbes(nargs, args), spectra(nargs, args)
  {
    initParsing(nargs, args);
    paramParsing();
  }

  Diagnostics::~Diagnostics()
  {}


#define PARSE(Fun)                                       \
bool Diagnostics::Fun() const                            \
{                                                        \
	if (Parser::Fun()) {                                   \
		hdf.Fun();                                           \
		phaseSpaces.Fun();                                   \
		moments.Fun();                                       \
		phiProbes.Fun();                                     \
		return true;                                         \
	}                                                      \
	return false;                                          \
}

  PARSE(parseHelp)
  PARSE(parseVersion)
  PARSE(parseTemplate)

#undef PARSE

  void Diagnostics::initParsing(int nargs, char *args[])
  {
    registerClass("Diagnostics");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

    parseLevelDebugOption("diagnostics::dl");

    using parser::types::boolean;
    using parser::types::intVect;

    insertOption(density  ,"density"  , boolean, "Enable/Disable density diag"  , Any(densityDiag));
    insertOption(potential,"potential", boolean, "Enable/Disable potential diag", Any(potentialDiag));
    insertOption(energy   ,"energy"   , boolean, "Enable/Disable energy diag"   , Any(energyDiag));
    insertOption(_iSpec   ,"ispec"    , intVect, "Diag iteration spec"          , Any(iSpec.data));
  }

  void Diagnostics::paramParsing()
  {
    parseOption(density  , densityDiag);
    parseOption(potential, potentialDiag);
    parseOption(energy   , energyDiag);
    parseOption(_iSpec   , iSpec.data);
  }

  void Diagnostics::startInit()
  {
    if (densityDiag ||
        potentialDiag ||
        energyDiag ||
        ! phaseSpaces.empty() ||
        ! moments.empty() ||
        ! phiProbes.empty() ||
        ! spectra.empty()
       )
      hdf.createSD(PACKAGE " " VERSION " " ID);
  }

  void Diagnostics::initFields(Scheduler &scheduler)
  {
    if (densityDiag) {
      prepareField(densityHdfName, scheduler);
      hdf.initSDvar(V);
    }
    if (potentialDiag && scheduler.isEforce()) {
      prepareField(potentialHdfName, scheduler);
      hdf.initSDvar(V);
    }
  }

  void Diagnostics::initEnergies(Scheduler &scheduler, VecReal &NPS,
                                 VecReal &KES, VecReal &KES_NPS)
  {
    if (energyDiag) {
      prepareVector(numParticleHdfName, static_cast<int>(NPS.size()), scheduler);
      hdf.initSDvar(V);
      prepareVector(keiHdfName, static_cast<int>(KES.size()), scheduler);
      hdf.initSDvar(V);
    }
  }

  void Diagnostics::initEnergies(Scheduler &scheduler, VecReal &ESE,
                                 VecReal &KEE)
  {
    if (energyDiag && scheduler.isEforce()) {
      prepareVector(keeHdfName, static_cast<int>(KEE.size()), scheduler);
      hdf.initSDvar(V);
      prepareVector(eseHdfName, static_cast<int>(ESE.size()), scheduler);
      hdf.initSDvar(V);
    }
  }

  void Diagnostics::initPhaseSpaces(Scheduler &scheduler)
  {
    PhaseSpaceListIter i=phaseSpaces.begin(), e=phaseSpaces.end();
    for ( ; i!=e; ++i) {
      (**i).initialise();
      preparePhaseSpace((**i).name(), scheduler, **i);
#if defined(HAVE_MPI)

      if (rankProc == masterProc)
#endif

        hdf.initSDvar(V);
    }
  }


  void Diagnostics::initMoments(Scheduler &scheduler)
  {
    MomentsListIter i=moments.begin(), e=moments.end();
    for ( ; i!=e; ++i) {
      (**i).initialise(scheduler);
      // density
      {
        prepareMoments((**i).name(), "n", scheduler, **i);
#if defined(HAVE_MPI)

        if (rankProc == masterProc)
#endif

          hdf.initSDvar(V);
      }
      // mean drift velocity
      for (int d=0; d<DIMV; ++d) {
        if ( (**i).isVel(d)) {
          prepareMoments((**i).name(), momsAxis[d], scheduler, **i);
#if defined(HAVE_MPI)

          if (rankProc == masterProc)
#endif

            hdf.initSDvar(V);
        }
      }
      // temperature
      for (int d=0; d<DIMV; ++d) {
        if ( (**i).isTemp(d)) {
          prepareMoments((**i).name(), momsAxis[DIMV+d], scheduler, **i);
#if defined(HAVE_MPI)

          if (rankProc == masterProc)
#endif

            hdf.initSDvar(V);
        }
      }
    }
  }


  void Diagnostics::initPotentialProbes(Scheduler &scheduler)
  {
    if ( scheduler.isEforce() && ! phiProbes.empty() ) {
      ProbesListIter i=phiProbes.begin(), e=phiProbes.end();
      for ( ; i!=e; ++i) {
        (**i).initialise(scheduler);
        preparePotentialProbes((**i).name(), scheduler, **i);
        hdf.initSDvar(V);
      }
    }
  }


  void Diagnostics::initSpectra(Scheduler &scheduler)
  {
    SpectraListIter i=spectra.begin(), e=spectra.end();
    for ( ; i!=e; ++i) {
      rvector St((**i).getTimeAxis(scheduler));
      rvector k((**i).getkc());
      rvector Sk((**i).getSk());
      rvector Sfrq((**i).getSfrq());

      if (scheduler.isEforce() && St.size()>0 && k.size()>0) {
        string vname((**i).name());
        vname += "k";
        string tname("t");
        tname += vname;
        string kname("k");
        kname += vname;
        initHdfVar(vname.c_str(),St,tname.c_str(),k,kname.c_str());
        hdf.initSDvar(V);
      }
      if (scheduler.isEforce() && St.size()>0 && Sk.size()>0 && Sfrq.size()>0) {
        string vname((**i).name());
        vname += "kf";
        string tname("t");
        tname += vname;
        string kname("k");
        kname += vname;
        string fname("f");
        fname += fname;
        initHdfVar(vname.c_str(),St,tname.c_str(),Sk,kname.c_str(),
                   Sfrq,fname.c_str());
        hdf.initSDvar(V);
      }
    }
  }

  void Diagnostics::endInit()
  {
#if 0
    if (densityDiag || potentialDiag || energyDiag ||
        ! phaseSpaces.empty() || ! phiProbes.empty() || ! spectra.empty() ) {
      hdf.endInitSD();
    }
#endif
  }


  void Diagnostics::prepareVector(const string &name, int dim,
                                  const Scheduler &scheduler)
  {
    V.resize(2);

    V.varname = name;
    V.vartype =  hdf::type<float>::def;

    V.dimsize[0] = scheduler.getMaxIter()+1;
    V.dimtype[0] = hdf::type<float>::def;
    V.dimname[0] = timeVectorHdfName;
    V.start[0] = 0;
    V.stride[0] = scheduler.getTimeInc();
    V.end[0] = scheduler.getTimeInc()*scheduler.getMaxIter();

    V.dimsize[1] = dim;
    V.dimtype[1] = hdf::type<float>::def;
    if (name ==  keeHdfName) {
      V.dimname[1] = "kee terms";
    } else if (name == eseHdfName) {
      V.dimname[1] = "ese terms";
    } else {
      V.dimname[1] = speciesHdfName;
    }
    V.start[1] = 0;
    V.stride[1] = 1.0;
    V.end[1] = static_cast<real>(dim)-1.0;
  }

  void Diagnostics::initHdfVar(const string &name,
                               const rvector &axis, const string &axisName)
  {
    V.resize(1);

    V.varname = name;
    V.vartype =  hdf::type<float>::def;

    V.dimsize[0] = axis.rows();
    V.dimtype[0] = hdf::type<float>::def;
    V.dimname[0] = axisName;
    V.start[0] = axis(0);
    V.stride[0] = (axis(V.dimsize[0]-1)-axis(0))/(V.dimsize[0]-1);
    V.end[0] = axis(V.dimsize[0]-1);
  }

  void Diagnostics::initHdfVar(const string &name,
                               const rvector &axis1, const string &axisName1,
                               const rvector &axis2, const string &axisName2)
  {
    V.resize(2);

    V.varname = name;
    V.vartype =  hdf::type<float>::def;

    V.dimsize[0] = axis1.rows();
    V.dimtype[0] = hdf::type<float>::def;
    V.dimname[0] = axisName1;
    V.start[0] = axis1(0);
    V.stride[0] = (axis1(V.dimsize[0]-1)-axis1(0))/(V.dimsize[0]-1);
    V.end[0] = axis1(V.dimsize[0]-1);

    V.dimsize[1] = axis2.rows();
    V.dimtype[1] = hdf::type<float>::def;
    V.dimname[1] = axisName2;
    V.start[1] = axis2(0);
    V.stride[1] = (axis2(V.dimsize[1]-1)-axis2(0))/(V.dimsize[1]-1);
    V.end[1] = axis2(V.dimsize[1]-1);
  }

  void Diagnostics::initHdfVar(
    const string &name,
    const rvector &axis1, const string &axisName1,
    const rvector &axis2, const string &axisName2,
    const	rvector &axis3, const string &axisName3)
  {
    V.resize(3);

    V.varname = name;
    V.vartype =  hdf::type<float>::def;

    V.dimsize[0] = axis1.rows();
    V.dimtype[0] = hdf::type<float>::def;
    V.dimname[0] = axisName1;
    V.start[0] = axis1(0);
    V.stride[0] = (axis1(V.dimsize[0]-1)-axis1(0))/(V.dimsize[0]-1);
    V.end[0] = axis1(V.dimsize[0]-1);

    V.dimsize[1] = axis2.rows();
    V.dimtype[1] = hdf::type<float>::def;
    V.dimname[1] = axisName2;
    V.start[1] = axis2(0);
    V.stride[1] = (axis2(V.dimsize[1]-1)-axis2(0))/(V.dimsize[1]-1);
    V.end[1] = axis2(V.dimsize[1]-1);

    V.dimsize[2] = axis3.rows();
    V.dimtype[2] = hdf::type<float>::def;
    V.dimname[2] = axisName3;
    V.start[2] = axis3(0);
    V.stride[2] = (axis3(V.dimsize[2]-1)-axis3(0))/(V.dimsize[2]-1);
    V.end[2] = axis3(V.dimsize[2]-1);
  }

  void Diagnostics::prepareField(const string &name, Scheduler &scheduler)
  {
    V.resize(1+DIMR);

    V.varname = name;
    V.vartype =  hdf::type<float>::def;

    V.dimsize[0] = (iSpec.end()-iSpec.start())/iSpec.stride()+1;
    V.dimtype[0] = hdf::type<float>::def;
    V.dimname[0] = timeFieldHdfName;
    V.start[0] = scheduler.getTimeInc()*iSpec.start();
    V.stride[0] = scheduler.getTimeInc()*iSpec.stride();
    V.end[0] = scheduler.getTimeInc()*iSpec.end();

    mudfas::FieldVeci gridsize(scheduler.getGridSize());
    Boundary<real,DIMR> domain(scheduler.getDomainBoundary());
    for (int d=0; d<DIMR; ++d) {
      int i = d+1;
      V.dimsize[i] = gridsize(d);
      V.dimtype[i] = hdf::type<float>::def;
      V.dimname[i] = rvAxis[d];
      V.start[i]   = domain.min(d);
      V.stride[i]  = (domain.max(d)-domain.min(d))/(gridsize(d)-1);
      V.end[i]     = domain.max(d);
    }
  }

  void Diagnostics::preparePhaseSpace(
    const string &name, Scheduler &scheduler, PhaseSpace &space)
  {
    V.resize(1+space.numDim());

    V.varname = name;
    V.vartype = hdf::type<float>::def;

    V.dimsize[0] = (iSpec.end()-iSpec.start())/iSpec.stride()+1;
    V.dimtype[0] = hdf::type<float>::def;
    V.dimname[0] = timeFieldHdfName;
    V.start[0] = scheduler.getTimeInc()*iSpec.start();
    V.stride[0] = scheduler.getTimeInc()*iSpec.stride();
    V.end[0] = scheduler.getTimeInc()*iSpec.end();

    for (int d=0, i=1; d<DIMRV; ++d) {
      if ( ! space.isIntegrated(d)) {
        V.dimtype[i] = hdf::type<float>::def;
        V.dimname[i] = rvAxis[d];
        V.dimname[i] = V.dimname[i] + "-" + name;
        V.start[i] = space.getStart(d);
        V.stride[i] = space.getStride(d);
        V.end[i] = space.getEnd(d);
        V.dimsize[i] = space.dim(d);
        ++i;
      }
    }
  }

  void Diagnostics::prepareMoments(
    const string &name, const string &moms_name, Scheduler &scheduler,
    Moments &moms)
  {
    V.resize(1+DIMR);

    V.varname = name;
    V.varname = V.varname + '-' + moms_name;
    V.vartype = hdf::type<float>::def;

    V.dimsize[0] = (iSpec.end()-iSpec.start())/iSpec.stride()+1;
    V.dimtype[0] = hdf::type<float>::def;
    V.dimname[0] = timeFieldHdfName;
    V.start[0] = scheduler.getTimeInc()*iSpec.start();
    V.stride[0] = scheduler.getTimeInc()*iSpec.stride();
    V.end[0] = scheduler.getTimeInc()*iSpec.end();

    for (int d=0, i=1; d<DIMR; ++d, ++i) {
      V.dimtype[i] = hdf::type<float>::def;
      V.dimname[i] = rvAxis[d];
      V.dimname[i] = V.dimname[i] + "-" + name;
      V.start[i] = moms.getStart(d);
      V.stride[i] = moms.getStride(d);
      V.end[i] = moms.getEnd(d);
      V.dimsize[i] = moms.dim(d);
    }
  }

  void Diagnostics::preparePotentialProbes(const string &name, Scheduler &scheduler,
      Probes &probes)
  {
    V.resize(1+DIMR);

    V.varname = name;
    V.vartype = hdf::type<float>::def;

    V.dimsize[0] = scheduler.getMaxIter()+1;
    V.dimtype[0] = hdf::type<float>::def;
    V.dimname[0] = timeVectorHdfName;
    V.start  [0] = 0;
    V.stride [0] = scheduler.getTimeInc();
    V.end    [0] = scheduler.getTimeInc()*scheduler.getMaxIter();

    for (int d=0, i=1; d<DIMR; ++d, ++i) {
      V.dimsize[i] = probes.size(d);
      V.dimtype[i] = hdf::type<float>::def;
      V.dimname[i] = rvAxis[d];
      V.dimname[i] = V.dimname[i] + "-" + name;
      V.start  [i] = probes.start(d);
      V.stride [i] = probes.stride(d);
      V.end    [i] = probes.end(d);
    }
  }

  void Diagnostics::saveFields(Scheduler &scheduler, Field &Rho, Field &Phi)
  {
    int cur_iter   = scheduler.getCurrentIter();
    int first_iter = iSpec.start();
    int strid_iter = iSpec.stride();
    int last_iter  = iSpec.end();
    if ((cur_iter-first_iter >= 0)  && (last_iter-cur_iter >= 0)  &&
        ((cur_iter-first_iter) % strid_iter == 0)) {

      mudfas::FieldVeci gridsize(scheduler.getGridSize());
      start.resize(1+DIMR);
      start[0] = (cur_iter-first_iter)/strid_iter;
      edge.resize(1+DIMR);
      edge[0] = 1;
      for (int d=1; d<=DIMR; ++d) {
        start[d] = 0;
        edge[d] = gridsize(d-1);
      }

      if (densityDiag) {
        hdf.writeSDvar(densityHdfName, start, edge, Rho.data());
      }

      if (potentialDiag && scheduler.isEforce()) {
        hdf.writeSDvar(potentialHdfName, start, edge, Phi.data());
      }
    }
  }

  void Diagnostics::saveEnergies(Scheduler &scheduler, VecReal &ESE,
                                 VecReal &KEE, VecReal &NPS, VecReal &KES)
  {
    int cur_iter(scheduler.getCurrentIter());
    if (energyDiag) {
      start.resize(2);
      edge.resize(2);
      start[0] = cur_iter;
      edge[0] = 1;
      start[1] = 0;
      edge[1] = static_cast<int>(NPS.size());

      hdf.writeSDvar(numParticleHdfName, start, edge, &NPS[0]);
      hdf.writeSDvar(keiHdfName, start, edge, &KES[0]);
      if (scheduler.isEforce()) {
        start.resize(2);
        edge.resize(2);
        start[1] = 0;
        edge[1] = 2;
        hdf.writeSDvar(keeHdfName, start, edge, &KEE[0]);
        start[1] = 0;
        edge[1] = 1;
        hdf.writeSDvar(eseHdfName, start, edge, &ESE[0]);
      }
    }
  }

  void Diagnostics::saveEnergies(Scheduler &scheduler, VecReal &NPS,
                                 VecReal &KES)
  {
    int cur_iter(scheduler.getCurrentIter());
    if (energyDiag) {
      start.resize(2);
      edge.resize(2);
      start[0] = cur_iter;
      edge[0] = 1;
      start[1] = 0;
      edge[1] = static_cast<int>(NPS.size());

      hdf.writeSDvar(numParticleHdfName, start, edge, &NPS[0]);
      hdf.writeSDvar(keiHdfName, start, edge, &KES[0]);
    }
  }

  void Diagnostics::savePhaseSpaces(Scheduler &scheduler, SpeciesHandler &species)
  {
    int cur_iter   = scheduler.getCurrentIter();
    int first_iter = iSpec.start();
    int strid_iter = iSpec.stride();
    int last_iter  = iSpec.end();
    if ((cur_iter-first_iter >= 0)  &&
        (last_iter-cur_iter >= 0)  &&
        ((cur_iter-first_iter) % strid_iter == 0)) {

      PhaseSpaceListIter i=phaseSpaces.begin(), e=phaseSpaces.end();
      for ( ; i!=e; ++i) {
        int ndims = (**i).numDim();
        RVVectori dim((**i).dim());
        start.resize(1+ndims);
        edge.resize(1+ndims);
        start[0] = (cur_iter-first_iter)/strid_iter;
        edge[0] = 1;
        for (int d=0, j=1; d<DIMRV; ++d) {
          if (dim(d) != 1) {
            start[j] = 0;
            edge[j] = dim(d);
            ++j;
          }
        }
        RVArrayr space(dim);
        (**i).computePhaseSpace(scheduler, species, space);
#if defined(HAVE_MPI)

        if (rankProc == masterProc)
#endif

          hdf.writeSDvar((**i).name(), start, edge, space.data());
      }
    }
  }

  void Diagnostics::saveMoments(Scheduler &scheduler, SpeciesHandler &species)
  {
    int cur_iter   = scheduler.getCurrentIter();
    int first_iter = iSpec.start();
    int strid_iter = iSpec.stride();
    int last_iter  = iSpec.end();
    if ((cur_iter-first_iter >= 0)  &&
        (last_iter-cur_iter >= 0)  &&
        ((cur_iter-first_iter) % strid_iter == 0)) {

      MomentsListIter i=moments.begin(), e=moments.end();
      for ( ; i!=e; ++i) {
        start.resize(1+DIMR);
        edge.resize(1+DIMR);
        start[0] = (cur_iter-first_iter)/strid_iter;
        edge[0] = 1;
        RVectori dim((**i).dim());
        for (int d=0, j=1; d<DIMR; ++d, ++j) {
          start[j] = 0;
          edge[j] = dim(d);
        }
        Field n(dim);
        VectorField u, T;
        for (int d=0; d<DIMV; ++d) {
          if ( (**i).isVel(d))
            u(d).resize(dim);
          if ( (**i).isTemp(d))
            T(d).resize(dim);
        }
        (**i).computeMoments(species, n, u, T);

        // density
        {
          string name((**i).name());
          name = name + '-' + 'n';
#if defined(HAVE_MPI)

          if (rankProc == masterProc)
#endif

            hdf.writeSDvar(name.c_str(), start, edge, n.data());
        }

        // mean drift velocity
        for (int d=0; d<DIMV; ++d) {
          if ( (**i).isVel(d)) {
            string name((**i).name());
            name = name + '-' + momsAxis[d];
#if defined(HAVE_MPI)

            if (rankProc == masterProc)
#endif

              hdf.writeSDvar(name.c_str(), start, edge, u(d).data());
          }
        }
        // temperature
        for (int d=0; d<DIMV; ++d) {
          if ( (**i).isTemp(d)) {
            string name((**i).name());
            name = name + '-' + momsAxis[DIMV+d];
#if defined(HAVE_MPI)

            if (rankProc == masterProc)
#endif

              hdf.writeSDvar(name.c_str(), start, edge, T(d).data());
          }
        }
      }
    }
  }


  void Diagnostics::savePotentialProbes(Scheduler &scheduler, Field &Phi)
  {
    int cur_iter(scheduler.getCurrentIter());
    if (scheduler.isEforce() && ! phiProbes.empty()) {
      start.resize(1+DIMR);
      edge.resize(1+DIMR);
      start[0] = cur_iter;
      edge [0] = 1;
      ProbesListIter i=phiProbes.begin(), e=phiProbes.end();
      for ( ; i!=e; ++i) {
        for (int d=0, j=1; d<DIMR; ++d, ++j) {
          start[j] = 0;
          edge [j] = (**i).size(d);
        }
#if (DIMR==2)
        Field phi((**i).size(0), (**i).size(1));
        phi = Phi((**i).range(0), (**i).range(1));
#elif (DIMR==3)

        Field phi((**i).size(0), (**i).size(1), (**i).size(2));
        phi = Phi((**i).range(0), (**i).range(1), (**i).range(2));
#endif

        hdf.writeSDvar((**i).name(), start, edge, phi.data());
      }
    }
  }

  void Diagnostics::saveSpectra(Scheduler &scheduler, Field &Rho)
  {
    SpectraListIter i=spectra.begin(), e=spectra.end();
    for ( ; i!=e; ++i) {
      (**i).update(Rho, scheduler);

      int it;
      rvector phik2;
      rmatrix phikw;

      if ( (**i).isSpectraComplete(scheduler) ) {

        (**i).getSpectra(it, phik2, phikw);

        {
          string vname((**i).name());
          vname += "k";

          start.resize(2);
          edge.resize(2);

          start[0] = it;
          edge [0] = 1;
          start[1] = 0;
          edge [1] = phik2.rows();

          hdf.writeSDvar(vname.c_str(), start, edge, phik2.data());
        }

        {
          string vname((**i).name());
          vname += "kf";

          start.resize(3);
          edge.resize(3);

          start[0] = it;
          edge [0] = 1;
          start[1] = 0;
          edge [1] = phikw.rows();
          start[2] = 0;
          edge [2] = phikw.cols();

          hdf.writeSDvar(vname.c_str(), start, edge, phikw.data());
        }
      }
    }
  }

#undef ID

}
