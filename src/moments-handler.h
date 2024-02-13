/**************************************************************************
 *
 * $Id: moments-handler.h,v 1.6 2011/03/26 15:36:08 patrick Exp $
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


#ifndef MOMENTS_HANDLER_H
#define MOMENTS_HANDLER_H

#include <moments.h>

namespace picsim {

  class MomentsHandler : public parser::Parser {
  public:

    typedef std::ostream ostream;

    typedef std::string string;
    typedef std::vector<Moments *>      MomentsList;
    typedef MomentsList::iterator       MomentsListIter;
    typedef MomentsList::const_iterator MomentsListConstIter;

    friend ostream& operator<<(ostream& os, const MomentsHandler &h);

    MomentsHandler(int nargs, char *args[]);
    ~MomentsHandler();

    virtual bool parseHelp()     const;
    virtual bool parseVersion()  const;
    virtual bool parseTemplate() const;

    void initialise(const Scheduler & scheduler);

    MomentsListIter      begin() {
      return moments.begin();
    }
    MomentsListConstIter begin()    const {
      return moments.begin();
    }

    MomentsListIter      end() {
      return moments.end();
    }
    MomentsListConstIter end()      const {
      return moments.end();
    }

    bool empty()                    const {
      return moments.empty();
    }

  protected:

    typedef parser::Parser Parser;

  private:

    enum parser_enum { _moments };

    MomentsList moments;
    NameList momentsNames;

    void insert(int nargs, char *args[], const string &name);

    void initParsing(int nargs, char *args[]);
    void paramParsing();

  };

}

#endif
