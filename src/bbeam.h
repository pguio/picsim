/**************************************************************************
 *
 * $Id: bbeam.h,v 1.2 2011/03/26 15:36:08 patrick Exp $
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

#ifndef BBEAM_H
#define BBEAM_H

#include <beam.h>

namespace picsim {

  class BBeam : public Beam {
  public:

    friend ostream& operator<<(ostream& os, const BBeam &b);

    BBeam(int nargs, char *args[], const string name="BBeam");
    virtual ~BBeam();

    virtual void initialise(NumberGenerator &draw, Scheduler &scheduler);
    virtual void injectParticles(NumberGenerator &draw, Scheduler &scheduler);

  protected:

    enum parser_enum { next=Beam::next };

    virtual void initialiseParticles(NumberGenerator &draw);
    virtual void printOn(ostream& os) const;

  private:

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };

}

#endif // BBEAM_H
