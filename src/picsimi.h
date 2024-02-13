/**************************************************************************
 *
 * $Id: picsimi.h,v 1.26 2011/03/26 15:36:08 patrick Exp $
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


#ifndef PICSIMI_H
#define PICSIMI_H

#include <generic-es-picsim.h>
#include <poisson-boltzmann.h>

namespace picsim {

  class PicSimI : public GenericEsPicSim<mudfas::PoissonBoltzmannSolver> {
  public:

    friend std::ostream& operator<<(std::ostream& os, const PicSimI &s);

    PicSimI(int nargs, char *args[]);
    virtual ~PicSimI();

  protected:

    VecReal ESE, KEE;

    virtual void solveFields();
    virtual void initDiagnostics();
    virtual void processDiagnostics();

  private:

    void initParsing(int nargs, char *args[]);
    void paramParsing();

  };

}

#endif // PICSIMI_H

