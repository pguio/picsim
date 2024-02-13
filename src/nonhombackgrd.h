/**************************************************************************
 *
 * $Id: nonhombackgrd.h,v 1.6 2011/03/26 15:36:08 patrick Exp $
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
 * o
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef NONHOMBACKGRD_H
#define NONHOMBACKGRD_H

#include <species.h>

namespace picsim {

  class NonHomBackgrd : public Species {
  public:

    enum imode_enum { conserving=0, statistic, noinject };

    friend ostream& operator<<(ostream& os, const NonHomBackgrd &b);

    NonHomBackgrd(int nargs=0, char *args[]=0, const string name="NonHomBackgrd");
    virtual ~NonHomBackgrd();

    virtual void initialise(NumberGenerator &draw, Scheduler &scheduler);
    virtual void injectParticles(NumberGenerator &draw, Scheduler &scheduler);

  protected:

    enum parser_enum { _imodemin=Species::next, _imodemax,
                       _icoefmin, _icoefmax,
                       _nhmfilename, _inhcoefmin, _inhcoefmax,
                       next
                     };

    Boundary<int ,DIMR> injectMode;
    Boundary<real,DIMR> injectCoef;
    Boundary<real,DIMR> injectNHCoef;
    Boundary<int ,DIMR> numInjectStat;
    Boundary<int ,DIMR> numInjectHomPerCellStat;

    LUT imodeMap;

    virtual void printOn(ostream& os) const;

    virtual void fixBoundaryParticles(Scheduler &scheduler);

    // non-homogeneous model primary variables
    std::string nhmFilename;
    RVectori gridSize;
    Field relDens;

    // derived variables
    RVectori numCells;
    Field cellRelDens;
    iField numParticlesPerCell;
    blitz::TinyVector< blitz::Array<int,DIMR-1> ,DIMR> numInjectPerCellAtLB;
    blitz::TinyVector< blitz::Array<int,DIMR-1> ,DIMR> numInjectPerCellAtUB;
    blitz::TinyVector< blitz::Array<int,DIMR-1> ,DIMR> LossAtLB;
    blitz::TinyVector< blitz::Array<int,DIMR-1> ,DIMR> LossAtUB;

  private:

    void loadNonHomModel();

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };

}

#endif // NONHOMBACKGRD_H
