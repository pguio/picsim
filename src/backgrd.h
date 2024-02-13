/**************************************************************************
 *
 * $Id: backgrd.h,v 1.69 2011/03/26 15:36:08 patrick Exp $
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

#ifndef BACKGRD_H
#define BACKGRD_H

#include <species.h>

namespace picsim {

  class Backgrd : public Species {
  public:

    enum imode_enum { conserving=0, statistic, noinject };

    friend ostream& operator<<(ostream& os, const Backgrd &b);

    Backgrd(int nargs=0, char *args[]=0, const string name="Backgrd");
    virtual ~Backgrd();

    virtual void initialise(NumberGenerator &draw, Scheduler &scheduler);
    virtual void injectParticles(NumberGenerator &draw, Scheduler &scheduler);

  protected:

    enum parser_enum { _imodemin=Species::next, _imodemax,
                       _icoefmin, _icoefmax, next
                     };

    Boundary<int ,DIMR> injectMode;
    Boundary<real,DIMR> injectCoef;
    Boundary<int ,DIMR> numInjectStat;

    LUT imodeMap;

    virtual void printOn(ostream& os) const;

  private:

    void initialiseUniform(NumberGenerator &draw);
    void initialiseExponential(NumberGenerator &draw);

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };

}

#endif // BACKGRD_H
