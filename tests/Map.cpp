/**************************************************************************
 *
 * $Id: Map.cpp,v 1.17 2011/03/26 15:36:08 patrick Exp $
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

#include <time.h>
#include <map>
#include <iostream>
using namespace std;

typedef map<int, double> Container;

int main(int nargs, char *args[])
{
  int n=10;
  if (nargs>1) n = atoi(args[1]);

  cout.precision(10);

  time_t tic;

  tic = time(0);
  // Initialise n random numbers in container
  Container container;
  int index = 0;
  Container::iterator i;
  Container::iterator last;
  srand48(1);
  cout << endl;
  for (index=0; index<n; index++) {
    container.insert(pair<int, double>(index,drand48()));
    cout << "index= " << index << " r= " << container[index] << endl;
  }
  cerr << "Init: " << time(0)-tic << " s" << endl;

  tic = time(0);
  index = 0;
  i = container.begin();
  cout << endl;
  do { // Remove elements with value outside [0.25 0.75[
    if ( (*i).second<=0.25 || (*i).second>0.75 ) {
      cout<< "index= " << index << " r= " << (*i).second << " removed." << endl;
      container.erase(i++);
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
    cout << "index= " << index++ << " r= " << (*i).second <<endl;
  }

  return 0;
}

