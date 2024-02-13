/**************************************************************************
 *
 * $Id: generic-es-picsim.h,v 1.18 2011/03/26 15:36:08 patrick Exp $
 *
 * Copyright (c) 2003-2011 Patrick Guio <patrick.guio@gmail.com>
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


#ifndef GENERIC_ES_PICSIM_H
#define GENERIC_ES_PICSIM_H

#include <generic-picsim.h>
#include <diagnostics.h>

namespace picsim {

  template <class Solver>
  class GenericEsPicSim : public GenericPicSim {
  public:

    typedef mudfas::Field Field;

    GenericEsPicSim(int nargs, char *args[]);
    virtual ~GenericEsPicSim();

    virtual bool parseHelp()     const;
    virtual bool parseVersion()  const;
    virtual bool parseTemplate() const;

    virtual void initialise();

    string getDiagnosticFilename() const {
      return diagnostics.getFilename();
    }

  protected:

    Solver solver;
    Diagnostics diagnostics;

    // Electrostatic field variables
    Field Rho, Phi;

    virtual void printOn(ostream& os) const;

    virtual void computeSources();
    virtual void solveFields() = 0;
    virtual void initDiagnostics() = 0;
    virtual void processDiagnostics() = 0;

    bool isEforce() const {
      return scheduler.isEforce();
    }

  };

}

#include <generic-es-picsim.cpp>

#endif // GENERIC_ES_PICSIM_H

