/**************************************************************************
 *
 * $Id: drivennonhombackgrd.h,v 1.2 2011/03/26 15:36:08 patrick Exp $
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

#ifndef DRIVENNONHOMBACKGRD_H
#define DRIVENNONHOMBACKGRD_H

#include <nonhombackgrd.h>

namespace picsim {

  class DrivenNonHomBackgrd : public NonHomBackgrd {
  public:

    friend ostream& operator<<(ostream& os, const DrivenNonHomBackgrd &d);

    DrivenNonHomBackgrd(int nargs, char *args[],
                        const string name="DrivenNonHomBackgrd");
    virtual ~DrivenNonHomBackgrd();

    virtual void injectParticles(NumberGenerator &draw, Scheduler &scheduler);

  protected:

    enum parser_enum { fmin=NonHomBackgrd::next, fmax, amin, amax, posmin, posmax,
                       devmin, devmax, perturbmin, perturbmax,
                       startmin, startmax, stopmin, stopmax, next
                     };

    enum perturb_enum { nope=0, harmonic, pulse };

    Boundary<real,DIMR> Frqcy;
    Boundary<real,DIMR> A;
    Boundary<real,DIMR> Pos;
    Boundary<real,DIMR> Dev;
    Boundary<int,DIMR> Perturbation;
    Boundary<real,DIMR> Start;
    Boundary<real,DIMR> Stop;

    LUT perturbMap;

    virtual void printOn(ostream& os) const;

  private:

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };

}

#endif // DRIVENNONHOMBACKGRD_H
