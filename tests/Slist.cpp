/**************************************************************************
 *
 * $Id: Slist.cpp,v 1.18 2011/03/26 15:36:08 patrick Exp $
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

#ifdef HAVE_SLIST

#include <time.h>
#include <slist>
#include <iostream>
using namespace std;

typedef slist<double> Container;

int main(int nargs, char *args[])
{
  int n=10;
  if (nargs>1) n = atoi(args[1]);

  cout.precision(10);

  time_t tic;

  tic = time(0);
  // Initialise n random numbers in container
  Container container(n);
  int index = 0;
  Container::iterator i=container.begin();
  Container::iterator last=container.end();
  srand48(1);
  cout << endl;
  for ( ; i != last; i++) {
    *i = drand48();
    cout << "index= " << index++ << " r= " << *i << endl;
  }
  cerr << "Init: " << time(0)-tic << " s" << endl;

  tic = time(0);
  index = 0;
  i = container.begin();
  cout << endl;
  do { // Remove elements with value outside [0.25 0.75[
    if ( *i<=0.25 || *i>0.75 ) {
      cout << "index= " << index << " r= " << *i << " removed." << endl;
      i = container.erase(i);
    } else {
      i++;
    }
    index++;
  } while (i != container.end());
  cerr << "Remove: " << time(0)-tic << " s" << endl;

  index = 0;
  i = container.begin();
  last = container.end();
  cout << endl;
  for ( ; i != last; i++) {
    cout << "index= " << index++ << " r= " << *i << endl;
  }
  return 0;
}
#else

#include <iostream>

using namespace std;


int main(int nargs, char *args[])
{
  cout << "(Cannot test " << __FILE__ <<
       " not supported on this platform)" << endl;

  return 0;
}

#endif
