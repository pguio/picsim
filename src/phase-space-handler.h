/**************************************************************************
 *
 * $Id: phase-space-handler.h,v 1.16 2011/03/26 15:36:08 patrick Exp $
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


#ifndef PHASE_SPACE_HANDLER_H
#define PHASE_SPACE_HANDLER_H

#include <phase-space.h>

namespace picsim {

  class PhaseSpaceHandler : public parser::Parser {
  public:

    typedef std::string string;
    typedef std::ostream ostream;

    typedef std::vector<PhaseSpace *> PhaseSpaceList;
    typedef PhaseSpaceList::iterator PhaseSpaceListIter;
    typedef PhaseSpaceList::const_iterator PhaseSpaceListConstIter;

    friend ostream& operator<<(ostream& os, const PhaseSpaceHandler &h);

    PhaseSpaceHandler(int nargs, char *args[]);
    ~PhaseSpaceHandler();

    virtual bool parseHelp()     const;
    virtual bool parseVersion()  const;
    virtual bool parseTemplate() const;

    void initialise();

    PhaseSpaceListIter      begin() {
      return spaces.begin();
    }
    PhaseSpaceListConstIter begin() const {
      return spaces.begin();
    }

    PhaseSpaceListIter      end() {
      return spaces.end();
    }
    PhaseSpaceListConstIter end()   const {
      return spaces.end();
    }

    bool empty()                    const {
      return spaces.empty();
    }

  protected:

    typedef parser::Parser Parser;

  private:

    enum parser_enum { _spaces };

    PhaseSpaceList spaces;
    NameList spaceNames;

    void insert(int nargs, char *args[], const string &name);

    void initParsing(int nargs, char *args[]);
    void paramParsing();

  };

}

#endif
