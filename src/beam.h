/**************************************************************************
 *
 * $Id: beam.h,v 1.70 2011/03/26 15:36:08 patrick Exp $
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

#ifndef BEAM_H
#define BEAM_H

#include <species.h>

namespace picsim {

  class Beam : public Species {
  public:

    friend ostream& operator<<(ostream& os, const Beam &b);

    Beam(int nargs, char *args[], const string name="Beam");
    virtual ~Beam();

    virtual void initialise(NumberGenerator &draw, Scheduler &scheduler);
    virtual void injectParticles(NumberGenerator &draw, Scheduler &scheduler);

  protected:

    enum parser_enum { _beamDir=Species::next,
                       _beamPos, _beamWidth, _beamInitDepth,
                       _timeSpec, _beamType, next
                     };

    int beamDir;
    RVectorr beamPos, beamWidth, beamInitDepth;
    VecInt timeSpec;

    real         getBeamArea(Scheduler &scheduler);
    virtual void initialiseParticles(NumberGenerator &draw);
    virtual void printOn(ostream& os) const;

  private:

    real beamInitDepthTime;

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };

}

#endif // BEAM_H
