/**************************************************************************
 *
 * $Id: surface-handler.h,v 1.13 2011/03/26 15:36:08 patrick Exp $
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


#ifndef SURFACE_HANDLER_H
#define SURFACE_HANDLER_H

#include <surface.h>

namespace picsim {

  class SurfaceHandler : public parser::Parser {
  public:

    typedef std::string string;
    typedef std::ostream ostream;

    typedef mudfas::Field Field;

    typedef std::vector<Surface *> SurfaceList;
    typedef SurfaceList::iterator SurfaceListIter;
    typedef SurfaceList::const_iterator SurfaceListConstIter;

    friend ostream& operator<<(ostream& os, const SurfaceHandler &s);

    SurfaceHandler(int nargs, char *args[]);
    ~SurfaceHandler();

    virtual bool parseHelp()     const;
    virtual bool parseVersion()  const;
    virtual bool parseTemplate() const;

    void initialise(const Scheduler &scheduler);

    void processSpecies(SpeciesHandler &species);
    void processDensity(Field &rho);

    SurfaceListIter      begin() {
      return surfaces.begin();
    }
    SurfaceListConstIter begin() const {
      return surfaces.begin();
    }

    SurfaceListIter      end() {
      return surfaces.end();
    }
    SurfaceListConstIter end()   const {
      return surfaces.end();
    }

    bool            empty()      const {
      return surfaces.empty();
    }

  protected:

    typedef parser::Parser Parser;

  private:

    enum parser_enum { _surfaces };

    SurfaceList surfaces;
    NameList surfaceNames;

    void insert(int nargs, char *args[], const string &name);

    void initParsing(int nargs, char *args[]);
    void paramParsing();
  };

}

#endif
