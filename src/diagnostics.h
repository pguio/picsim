/**************************************************************************
 *
 * $Id: diagnostics.h,v 1.104 2011/03/26 15:36:08 patrick Exp $
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

#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include <scheduler.h>
#include <species-handler.h>
#include <phase-space-handler.h>
#include <moments-handler.h>
#include <probes-handler.h>
#include <spectra-handler.h>
#include <hdf-interface.h>

namespace picsim {

  class Diagnostics : public parser::Parser {
  public:

    typedef std::ostream ostream;

    typedef mudfas::Field Field;

    typedef PhaseSpaceHandler::PhaseSpaceList PhaseSpaceList;
    typedef PhaseSpaceHandler::PhaseSpaceListIter PhaseSpaceListIter;
    typedef PhaseSpaceHandler::PhaseSpaceListConstIter PhaseSpaceListConstIter;

    typedef Moments::VectorField VectorField;
    typedef MomentsHandler::MomentsList MomentsList;
    typedef MomentsHandler::MomentsListIter MomentsListIter;
    typedef MomentsHandler::MomentsListConstIter MomentsListConstIter;

    typedef ProbesHandler::ProbesList ProbesList;
    typedef ProbesHandler::ProbesListIter ProbesListIter;
    typedef ProbesHandler::ProbesListConstIter ProbesListConstIter;

    typedef SpectraHandler::SpectraList SpectraList;
    typedef SpectraHandler::SpectraListIter SpectraListIter;
    typedef SpectraHandler::SpectraListConstIter SpectraListConstIter;

    typedef Spectra::rvector rvector;
    typedef Spectra::rmatrix rmatrix;

    friend ostream& operator<<(ostream& os, const Diagnostics &d);

    Diagnostics(int nargs, char *args[]);
    ~Diagnostics();

    virtual bool parseHelp() const;
    virtual bool parseVersion() const;
    virtual bool parseTemplate() const;

    void startInit();
    void initFields(Scheduler &scheduler);
    void initEnergies(Scheduler &scheduler, VecReal &NPS, VecReal &KES,
                      VecReal &KES_NPS);
    void initEnergies(Scheduler &scheduler, VecReal &ESE, VecReal &KEE);
    void initPhaseSpaces(Scheduler &scheduler);
    void initMoments(Scheduler &scheduler);
    void initPotentialProbes(Scheduler &scheduler);
    void initSpectra(Scheduler &scheduler);
    void endInit();

    void saveFields(Scheduler &scheduler, Field &Rho, Field &Phi);
    void saveEnergies(Scheduler &scheduler, VecReal &ESE,
                      VecReal &KEE, VecReal &NPS, VecReal &KES);
    void saveEnergies(Scheduler &scheduler,  VecReal &NPS, VecReal &KES);
    void savePhaseSpaces(Scheduler &scheduler, SpeciesHandler &species);
    void saveMoments(Scheduler &scheduler, SpeciesHandler &species);
    void savePotentialProbes(Scheduler &scheduler, Field &Phi);
    void saveSpectra(Scheduler &scheduler, Field &Rho);

    string getFilename() const {
      return hdf.getFilename();
    }

  protected:

    typedef parser::Parser Parser;

  private:

    enum parser_enum { density=1, potential, energy, _iSpec };

    enum dim_enum { time, x, y, z};

    hdf::Hdf hdf;
    hdf::SDvar V;
    hdf::VecInt start;
    hdf::VecInt edge;

    RangeSpec<int> iSpec;

    bool densityDiag, potentialDiag, energyDiag;

    PhaseSpaceHandler phaseSpaces;
    MomentsHandler    moments;
    ProbesHandler     phiProbes;
    SpectraHandler    spectra;

    void initParsing(int nargs, char *args[]);
    void paramParsing();

    void initHdfVar(const string &name,
                    const rvector &axis1, const string &axisName1);
    void initHdfVar(const string &name,
                    const rvector &axis1, const string &axisName1,
                    const rvector &axis2, const string &axisName2);
    void initHdfVar(const string &name,
                    const rvector &axis1, const string &axisName1,
                    const rvector &axis2, const string &axisName2,
                    const rvector &axis3, const string &axisName3);

    void prepareVector(const string &name, int dim, const Scheduler &scheduler);
    void prepareField(const string &name, Scheduler &scheduler);
    void preparePhaseSpace(const string &name, Scheduler &scheduler,
                           PhaseSpace &space);
    void prepareMoments(const string &name, const string &moms_name,
                        Scheduler &scheduler, Moments &moments);
    void preparePotentialProbes(const string &name, Scheduler &scheduler,
                                Probes &probes);
  };

}

#endif // DIAGNOSTICS_H
