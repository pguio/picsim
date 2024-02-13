/**************************************************************************
 *
 * $Id: spectra-handler.h,v 1.5 2011/03/26 15:36:08 patrick Exp $
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


#ifndef SPECTRA_HANDLER_H
#define SPECTRA_HANDLER_H

#include <spectra.h>

namespace picsim {

  class SpectraHandler : public parser::Parser {
  public:

    typedef std::ostream ostream;

    typedef std::string string;
    typedef std::vector<Spectra *> SpectraList;
    typedef SpectraList::iterator SpectraListIter;
    typedef SpectraList::const_iterator SpectraListConstIter;

    friend ostream& operator<<(ostream& os, const SpectraHandler &s);

    SpectraHandler(int nargs, char *args[]);
    ~SpectraHandler();

    virtual bool parseHelp() const;
    virtual bool parseVersion() const;
    virtual bool parseTemplate() const;

    SpectraListIter      begin() {
      return spectra.begin();
    }
    SpectraListConstIter begin() const {
      return spectra.begin();
    }

    SpectraListIter      end() {
      return spectra.end();
    }
    SpectraListConstIter end()   const {
      return spectra.end();
    }

    bool                 empty() const {
      return spectra.empty();
    }

  protected:

    typedef parser::Parser Parser;

  private:

    enum parser_enum { _spectra };

    SpectraList spectra;
    NameList spectraNames;

    void insert(int nargs, char *args[], const string &name);

    void initParsing(int nargs, char *args[]);
    void paramParsing();

  };

}

#endif
