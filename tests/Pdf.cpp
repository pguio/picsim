/**************************************************************************
 *
 * $Id: Pdf.cpp,v 1.34 2011/03/26 15:36:08 patrick Exp $
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

#include <new>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <blitz/array.h>
#include <pdf.h>

using blitz::Array;


void my_new_handler()
{
  cerr << "Out of memory" << endl;
  abort();
}

#if defined(FLOAT_COORDINATE)
typedef float REAL;
#else
typedef double REAL;
#endif

template <typename T_num>
T_num stdd(Array<T_num,1> &r)
{
  T_num m(mean(r));
  return sqrt(sum(sqr(r-m))/(r.numElements()-1.0));
}

int main(int nargs, char *args[])
{
  try {

    BZ_STD_SCOPE(set_new_handler(my_new_handler));
    cout.precision(10);
    enum { _nb=1, _a, _b, _mu, _sigma, _lambda, _m, _s, _alpha, _quiet };

    int nb(10);
    REAL a(0.0), b(1.0);
    REAL mu(1.0), sigma(1.0);
    REAL lambda(1.0);
    REAL m(0.5), s(0.1), alpha(0.1);
    bool quiet(false);

    parser::Parser parser(nargs, args);
    parser.registerProgram(args[0]);
    parser.registerPackage(PACKAGE, VERSION, PDF_COPYRIGHT);

    using parser::types::integer;
    using parser::types::real;
    using parser::types::none;

    parser.insertOption(_nb,"-nb",integer,"Number of random numbers",Any(nb));
    parser.insertOption(_a,"a",real,"a / U(a,b), P(a,b,m,s,alpha)",Any(a));
    parser.insertOption(_b,"b",real,"b / U(a,b), P(a,b,m,s,alpha)",Any(b));
    parser.insertOption(_mu, "mu", real,
                        "mu / N(mu,sigma), R(sigma), F(mu,sigma)",Any(mu));
    parser.insertOption(_sigma, "sigma", real,
                        "sigma / N(mu,sigma), R(sigma), F(mu,sigma)",
                        Any(sigma));
    parser.insertOption(_lambda, "lambda", real,
                        "lambda / E(mu,sigma)",Any(lambda));
    parser.insertOption(_m, "m", real, "m / P(a,b,m,s,alpha)",Any(m));
    parser.insertOption(_s, "s", real, "s / P(a,b,m,s,alpha)",Any(s));
    parser.insertOption(_alpha, "alpha", real,
                        "alpha / P(a,b,m,s,alpha)",Any(alpha));
    parser.insertOption(_quiet, "-q", none, "quiet mode", Any());

    parser.parseOption(_nb, nb);
    parser.parseOption(_a, a);
    parser.parseOption(_b, b);
    parser.parseOption(_mu, mu);
    parser.parseOption(_sigma, sigma);
    parser.parseOption(_lambda, lambda);
    parser.parseOption(_m, m);
    parser.parseOption(_s, s);
    parser.parseOption(_alpha, alpha);
    if (parser.parseOption(_quiet)) quiet = !quiet;

    pdf::NumberGenerator<REAL> Draw(nargs, args);

#define PARSE(Fun)      \
if (parser.Fun()) {     \
	Draw.Fun();           \
	return 0;             \
}                       \
 
    PARSE(parseHelp)
    PARSE(parseVersion)
    PARSE(parseTemplate)

#undef PARSE

    Draw.initialise();

    cout << Draw << endl;

    Array<REAL,1> Uniform(nb), Normal(nb), Rayleigh(nb), Exponential(nb),
          Flux(nb), Perturbp(nb), Perturbm(nb);

    for (int i=0; i<nb; i++) {
      //cout << "i = " << i << endl;
      Uniform(i) = Draw.uniform(a, b);
      Normal(i) = Draw.normal(mu, sigma);
      Rayleigh(i) =	Draw.rayleigh(sigma);
      Exponential(i) = Draw.exponential(lambda);
      Flux(i) = Draw.flux(mu, sigma);
      Perturbp(i) = Draw.perturb(a, b, m, s, alpha);
      Perturbm(i) = Draw.perturb(a, b, m, s, -alpha);
    }

    cout.precision(4);
    cout << "U(" << a << ',' << b << ")\t"
         << "N(" << mu << ',' << sigma << ")\t"
         << "R(" << sigma << ")\t"
         << "E(" << lambda << ")\t"
         << "F(" << mu << ',' << sigma << ")\t"
         << "P(" << a << ',' << b << ',' << m << ',' << s << ',' << alpha << ")\t"
         << "P(" << a << ',' << b << ',' << m << ',' << s<< ',' << -alpha<< ")\t"
         << endl;

    if ( !quiet ) {
      cout << string(80, '=') << endl;
      for (int i=0; i<nb; i++) {
        cout << Uniform(i) << '\t' << Normal(i) << '\t'
             << Rayleigh(i) << '\t' << Exponential(i) << '\t'
             << Flux(i) << '\t' << Perturbp(i) << '\t' << Perturbm(i) << endl;
      }
    }

    cout.precision(3);
    cout << string(80,'=') << endl;
    cout << "mean "
         << mean(Uniform) << '\t' << mean(Normal) << '\t'
         << mean(Rayleigh) << '\t' << mean(Exponential) << '\t'
         << mean(Flux) << '\t' << mean(Perturbp) << '\t'
         << mean(Perturbm) << endl;
    cout << "std  "
         << stdd(Uniform) << '\t' << stdd(Normal) << '\t'
         << stdd(Rayleigh) << '\t' << stdd(Exponential) << '\t'
         << stdd(Flux) << '\t' << stdd(Perturbp) << '\t'
         << stdd(Perturbm) << endl;

    return 0;

  } catch(ClassException& x) {
    cerr << x.what() << endl;
    return !0;
  } catch(int status) {
    return status;
  }
}

