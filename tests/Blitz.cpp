/**************************************************************************
 *
 * $Id: Blitz.cpp,v 1.29 2011/03/26 15:36:08 patrick Exp $
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

#include <iostream>
#include <fstream>
using namespace std;

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
using namespace blitz;
using namespace blitz::tensor;

#include <bench.h>


inline
TinyVector<int,3> getVector()
{
  return TinyVector<int,3>(1,2,3);
}

inline
Array<float,2> LinearInterpolation2d0(float dx, float dy)
{
  Array<float,1> x(2), y(2);
  x = 1.0-dx, dx;
  y = 1.0-dy, dy;
  Array<float,2> xy(2,2);
  xy = x(i)*y(j);
  return xy;
}

inline
Array<float,2> LinearInterpolation2d1(float dx, float dy)
{
  Array<float,2> xy(2,2);
  float cdx(1.0-dx);
  float cdy(1.0-dy);
  xy =
    cdx*cdy, cdx*dy,
    dx*cdy , dx*dy;
  return xy;
}

inline
TinyVector<float,4> LinearInterpolation2d2(float dx, float dy)
{
  float cdx(1.0-dx);
  float cdy(1.0-dy);
  return TinyVector<float,4>(cdx*cdy, cdx*dy,
                             dx*cdy , dx*dy);
}

inline
Array<float,3> LinearInterpolation3d0(float dx, float dy, float dz)
{
  Array<float,1> x(2), y(2), z(2);
  x = 1.0-dx, dx;
  y = 1.0-dy, dy;
  z = 1.0-dz, dz;
  Array<float,3> xyz(2,2,2);
  xyz = x(i)*(y(j)*z(k));
  return xyz;
}

inline
Array<float,3> LinearInterpolation3d1(float dx, float dy, float dz)
{
  Array<float,3> xyz(2,2,2);
  float cdx(1.0-dx);
  float cdy(1.0-dy);
  float cdz(1.0-dz);
  xyz =
    cdx*cdy*cdz, cdx*cdy*dz,
    cdx*dy*cdz , cdx*dy*dz,
    dx*cdy*cdz , dx*cdy*dz,
    dx*dy*cdz  , dx*dy*dz;
  return xyz;
}

inline
TinyVector<float,8> LinearInterpolation3d2(float dx, float dy, float dz)
{
  float cdx(1.0-dx);
  float cdy(1.0-dy);
  float cdz(1.0-dz);
  return TinyVector<float,8>(cdx*cdy*cdz, cdx*cdy*dz,
                             cdx*dy*cdz , cdx*dy*dz,
                             dx*cdy*cdz , dx*cdy*dz,
                             dx*dy*cdz  , dx*dy*dz);
}

void FixBoundary0(Array<float,3> &Rho)
{
  Range all(Range::all());
  for (int d=0; d<3; d++) {
    switch (d) {
    case 0:
      Rho(0,all,all) += 2.0;
      break;
    case 1:
      Rho(all,0,all) += 2.0;
      break;
    case 2:
      Rho(all,all,0) += 2.0;
      break;
    }
    switch (d) {
    case 0:
      Rho(Rho.ubound(d),all,all) += 2.0;
      break;
    case 1:
      Rho(all,Rho.ubound(d),all) += 2.0;
      break;
    case 2:
      Rho(all,all,Rho.ubound(d)) += 2.0;
      break;
    }
  }
}

void FixBoundary1(Array<float,3> &Rho)
{
  Range all(Range::all());
  for (int d=0; d<3; d++) {
    Rho(0,all,all) += 2.0;
    Rho(Rho.ubound(0),all,all) += 2.0;
    Rho.transposeSelf(secondDim,thirdDim,firstDim);
  }
}

#define blockA0 {                            \
	for (unsigned n=0; n<nli; n++)             \
		A0 += LinearInterpolation2d0(dx, dy);    \
}

#define blockA1 {                             \
	for (unsigned n=0; n<nli; n++)              \
		A1 += LinearInterpolation2d1(dx, dy);     \
}

#define blockB0 {                             \
	for (unsigned n=0; n<nli; n++)              \
		B0 += LinearInterpolation3d0(dx, dy, dz); \
}

#define blockB1 {                             \
	for (unsigned n=0; n<nli; n++)              \
		B1 += LinearInterpolation3d1(dx, dy, dz); \
}

#define blockRho0 {                           \
	for (unsigned n=0; n<nbfix; n++)            \
		FixBoundary0(Rho0);                       \
}

#define blockRho1 {                           \
	for (unsigned n=0; n<nbfix; n++)            \
		FixBoundary0(Rho1);                       \
}

#define blockdot1 {                           \
	for (unsigned n=0; n<ndot; n++)             \
		float res = dot(V1, V2);                  \
}

#define blockdot2 {                           \
	for (unsigned n=0; n<ndot; n++)             \
		float res = sum(V1*V2);                   \
}


int main(int nargs, char *args[])
{
  try {
    const float dx=0.1, dy=0.25, dz=0.6;
    const unsigned nli=1000;
    const unsigned nbfix=10;
    const unsigned ndot=1000;
    double t;
    bench::verbose = true;

    Array<float,2> A0(2,2), A1(2,2);
    BENCH(blockA0, "A0=LinearInterpolation2d0(dx, dy)",5,t)
    BENCH(blockA1, "A1=LinearInterpolation2d1(dx, dy)",5,t)
    A0 = LinearInterpolation2d0(dx, dy);
    A1 = LinearInterpolation2d1(dx, dy);
    if (any(A0-A1 != 0.0)) {
      cerr << "A0 != A1" << endl
           << "A0 = " << A0 << endl << "A1 = " << A1
           << "A0-A1 = " << Array<float,2>(abs(A0-A1)) << endl;
#if 0
      throw EXIT_FAILURE;
#endif
    }

    Array<float,3> B0(2,2,2), B1(2,2,2);
    BENCH(blockB0, "B0=LinearInterpolation3d0(dx, dy, dz)",5,t)
    BENCH(blockB1, "B1=LinearInterpolation3d1(dx, dy, dz)",5,t)
    B0 = LinearInterpolation3d0(dx, dy, dz);
    B1 = LinearInterpolation3d1(dx, dy, dz);
    if (any(B0-B1 != 0.0)) {
      cerr << "B0 != B1" << endl
           << "B0 = " << B0 << endl << "B1 = " << B1
           << "B0-B1 = " << Array<float,3>(abs(B0-B1)) << endl;
#if 0
      throw EXIT_FAILURE;
#endif
    }

    Array<float,3> Rho0(50,200,50), Rho1(50,200,50);
    Rho0 = 1.0;
    Rho1 = 1.0;
    BENCH(blockRho0, "FixBoundary0(Rho0)",5,t)
    BENCH(blockRho1, "FixBoundary1(Rho1)",5,t)
    Rho0 = 1.0;
    Rho1 = 1.0;
    FixBoundary0(Rho0);
    FixBoundary1(Rho1);
    if (any(Rho0-Rho1 != 0.0)) {
      cerr << "Rho0 != Rho1" << endl
           << "Rho0-Rho1 = " << Array<float,3>(abs(Rho0-Rho1)) << endl;
#if 0
      throw EXIT_FAILURE;
#endif
    }

    TinyVector<int,3> V1;
    V1 = getVector();
    TinyVector<int,3> V2(getVector());
    cout << "V1 = " << V1 << endl;
    cout << "V2 = " << V2 << endl;
    cout << "getVector() = " << getVector() << endl;

    BENCH(blockdot1, "dot(V1,V2)",5,t)
    BENCH(blockdot2, "sum(V1*V2)",5,t)
    float dot1 = dot(V1, V2);
    float dot2 = sum(V1* V2);
    cout << "dot(V1,V2)  = " << dot1 << endl;
    cout << "sum(V1*V2)  = " << dot2 << endl;
    if ( dot1-dot2 != 0.0) {
      cerr << "dot1 != dot2" << endl << "dot1-do2 = " << dot1-dot2 << endl;
    }

    Array<float,3> C(60,14,32);
    Array<float,2> D(14,32);
    C = 1.0;
    D = sum(C(k,i,j),k);
    if (any(D != 60.0)) {
      cerr << "D != 60.0" << endl;
#if 0
      throw EXIT_FAILURE;
#endif
    }


  } catch(int status) {
    return status;
  }
  return 0;
}
