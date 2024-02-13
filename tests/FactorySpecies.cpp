/**************************************************************************
 *
 * $Id: FactorySpecies.cpp,v 1.7 2011/03/26 15:36:08 patrick Exp $
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


#include <species.h>
#include <species-factory.h>


typedef picsim::SpeciesFactory<picsim::Species> Factory;


int main(int nargs, char *args[])
{

  try {
    Factory::instance().init();

    Factory::instance().content();

    Factory::Species_iterator i = Factory::instance().begin();
    Factory::Species_iterator e = Factory::instance().end();

    for (; i!= e ; i++) {
      picsim::defaultIDKeyType key(Factory::instance().getKey(i));
      std::cout << key << std::endl;
      Factory::BasePtr pObj(Factory::instance().create(key,nargs,args));
      std::cout << *pObj << std::endl;
    }

    Factory::BasePtr pObj(Factory::instance().create("foo",nargs,args,"foo"));
    std::cout << *pObj << std::endl;

    return 0;
  } catch(ClassException& c) {
    std::cout << c.what() << std::endl;
    return !0;
  }
}
