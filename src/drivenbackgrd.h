/**************************************************************************
 *
 * $Id: drivenbackgrd.h,v 1.43 2017/11/26 18:22:11 patrick Exp $
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

#ifndef DRIVENBACKGRD_H
#define DRIVENBACKGRD_H

#include <backgrd.h>

namespace picsim {

  class DrivenBackgrd : public Backgrd {
  public:

    friend ostream& operator<<(ostream& os, const DrivenBackgrd &d);

    DrivenBackgrd(int nargs, char *args[], const string name="DrivenBackgrd");
    virtual ~DrivenBackgrd();

    virtual void initialise(NumberGenerator &draw, Scheduler &scheduler);
    virtual void injectParticles(NumberGenerator &draw, Scheduler &scheduler);

  protected:

    enum parser_enum { fmin=Backgrd::next, fmax, amin, amax, posmin, posmax,
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

#endif // DRIVENBACKGRD_H
