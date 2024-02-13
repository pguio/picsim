/**************************************************************************
 *
 * $Id: picsime.h,v 1.25 2011/03/26 15:36:08 patrick Exp $
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


#ifndef PICSIME_H
#define PICSIME_H

#include <generic-es-picsim.h>
#include <poisson.h>

namespace picsim {

  class PicSimE : public GenericEsPicSim<mudfas::PoissonSolver> {
  public:

    friend std::ostream& operator<<(std::ostream& os, const PicSimE &s);

    PicSimE(int nargs, char *args[]);
    virtual ~PicSimE();

  protected:

    virtual void solveFields();
    virtual void initDiagnostics();
    virtual void processDiagnostics();

  private:

    void initParsing(int nargs, char *args[]);
    void paramParsing();

  };

}

#endif // PICSIME_H

