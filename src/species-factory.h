/**************************************************************************
 *
 * $Id: species-factory.h,v 1.16 2011/03/26 15:36:08 patrick Exp $
 *
 * Copyright (c) 2004-2011 Patrick Guio <patrick.guio@gmail.com>
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


#ifndef SPECIES_FACTORY_H
#define SPECIES_FACTORY_H

#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <classexception.h>

namespace picsim {

  typedef std::string defaultIDKeyType;

  namespace factory {
    void dummyBgrd();
    void dummyNonHomBgrd();
    void dummyBeam();
    void dummyDrivenBgrd();
    void dummyDrivenNonHomBgrd();
    void dummyBBeam();
  }

  template <class ManufacturedType, typename ClassIDKey=defaultIDKeyType>
  class SpeciesFactory {
  public:
    typedef ManufacturedType* BasePtr;
  private:

    typedef BasePtr (*BaseCreateFn)(int , char **, const std::string);
    typedef std::map<ClassIDKey, BaseCreateFn> FnRegistry;
    FnRegistry registry;
    SpeciesFactory() {}
    SpeciesFactory(const SpeciesFactory&); // Not implemented
    SpeciesFactory &operator=(const SpeciesFactory&); // Not implemented
  public:
    typedef ClassIDKey IDKeyType;
    typedef typename FnRegistry::const_iterator Species_iterator;

    static SpeciesFactory &instance() {
      static SpeciesFactory<ManufacturedType, ClassIDKey> bf;
      return bf;
    }

    void RegCreateFn(const ClassIDKey & id, BaseCreateFn fn) {
      registry[id] = fn;
    }

    BasePtr create(const ClassIDKey &className,
                   int nargs, char *args[],
                   const std::string name="") const {
      BasePtr theObject(0);
      ClassIDKey objName(name.empty() ? className : name);
      Species_iterator regEntry = registry.find(className);
      if (regEntry != registry.end()) {
        theObject = regEntry->second(nargs,args,objName);
      } else {
        throw ClassException("SpeciesFactory", "Error: unknown ClassIDKey: " + className);
      }
      return theObject;
    }
    Species_iterator begin() const {
      return registry.begin();
    }
    Species_iterator end  () const {
      return registry.end();
    }
    ClassIDKey getKey(Species_iterator &i) const {
      return i->first;
    }

    void content() const {
      using std::cout;
      using std::endl;
      int tab1=30, tab2 =10;

      std::ios::fmtflags f = cout.flags() & std::ios::adjustfield;
      cout << "Factory content:\n" << std::string(tab1+tab2, '=') << '\n'
           << std::setw(tab1) << std::left << "Name"
           << std::setw(tab2) << std::right << "Address" << '\n'
           << std::string(tab1+tab2, '=') << '\n';
      for (Species_iterator i = registry.begin(); i != registry.end(); ++i) {
        cout
            << std::setw(tab1) << std::left << i->first
            << std::setw(tab2) << std::right << std::hex << (unsigned long)(i->second)
            << endl;
      }
      cout << std::string(tab1+tab2, '=') << endl
           << std::resetiosflags(std::ios::left)
           << std::resetiosflags(std::ios::right)
           << std::setiosflags(f);
    }
    void init() const {
      factory::dummyBgrd();
      factory::dummyNonHomBgrd();
      factory::dummyBeam();
      factory::dummyDrivenBgrd();
      factory::dummyDrivenNonHomBgrd();
      factory::dummyBBeam();
    }
  };



  template <class AncestorType,
           class ManufacturedType,
           typename ClassIDKey=defaultIDKeyType>
  class RegisterInFactory {
  private:
    typedef SpeciesFactory<AncestorType, ClassIDKey> Factory;
    typedef typename Factory::BasePtr BasePtr;
  public:
    static BasePtr CreateInstance(int nargs, char *args[], const std::string name) {
      return BasePtr(new ManufacturedType(nargs, args, name));
    }

    RegisterInFactory(const ClassIDKey &id) {
      Factory::instance().RegCreateFn(id, CreateInstance);
    }
  };

}

#endif
