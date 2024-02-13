/**************************************************************************
 *
 * $Id: Fftw.cpp,v 1.12 2011/03/26 15:36:08 patrick Exp $
 *
 * Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
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

#include <iostream>
using namespace std;

#if defined(HAVE_CSTDLIB)
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include <fourier.h>

typedef fourier::ODFT1D<fourier::complex,fourier::complex> FT1d;
using fourier::ifftshift;
using fourier::fftshift;

using blitz::Array;
using blitz::tensor::i;

int main(int nargs, char *args[])
{
  try {
    int n=64;
    if (nargs == 2) n = atoi(args[1]);

    FT1d ft1d(n);
    cout << ft1d << endl;

    Array<fourier::complex, 1> A(n), B(n);

    A=1.0;
    cout << "A=" << A << endl;

    ft1d.direct(A, B);
    cout << "B=" << B << endl;

    ft1d.inverse(B, A);
    cout << "A=" << A << endl;

    Array<int,1> signed_index(n);
    signed_index = i-(n>>1);

    cout << "signed_index = " << signed_index << endl;

    signed_index = ifftshift(signed_index);
    cout << "signed_index = " << signed_index << endl;

    signed_index = fftshift(signed_index);
    cout << "signed_index = " << signed_index << endl;

    signed_index.resize(n+1);
    signed_index = i-(n>>1);

    cout << "signed_index = " << signed_index << endl;

    signed_index = ifftshift(signed_index);
    cout << "signed_index = " << signed_index << endl;

    signed_index = fftshift(signed_index);
    cout << "signed_index = " << signed_index << endl;

  } catch(int status) {
    return status;
  } catch(ClassException& c) {
    cerr << c.what() << endl;
    return EXIT_FAILURE;
  }

}
