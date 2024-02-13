/**************************************************************************
 *
 * $Id: probes-handler.h,v 1.6 2011/03/26 15:36:08 patrick Exp $
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


#ifndef PROBES_HANDLER_H
#define PROBES_HANDLER_H

#include <probes.h>

namespace picsim {

  class ProbesHandler : public parser::Parser {
  public:

    typedef std::string string;
    typedef std::ostream ostream;

    typedef std::vector<Probes *> ProbesList;
    typedef ProbesList::iterator ProbesListIter;
    typedef ProbesList::const_iterator ProbesListConstIter;

    friend ostream& operator<<(ostream& os, const ProbesHandler &p);

    ProbesHandler(int nargs, char *args[]);
    ~ProbesHandler();

    virtual bool parseHelp()     const;
    virtual bool parseVersion()  const;
    virtual bool parseTemplate() const;

    ProbesListIter      begin() {
      return probes.begin();
    }
    ProbesListConstIter begin() const {
      return probes.begin();
    }

    ProbesListIter      end() {
      return probes.end();
    }
    ProbesListConstIter end()   const {
      return probes.end();
    }

    bool                empty() const {
      return probes.empty();
    }

  protected:

    typedef parser::Parser Parser;

  private:

    enum parser_enum { _probes };

    ProbesList probes;
    NameList probesNames;

    void insert(int nargs, char *args[], const string &name);

    void initParsing(int nargs, char *args[]);
    void paramParsing();

  };

#endif

}
