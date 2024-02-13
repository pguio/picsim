/**************************************************************************
 *
 * $Id: pdf.cpp,v 1.82 2024/02/13 19:28:18 patrick Exp $
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

#include <pdf.h>
#define DEBUG_LEVEL
#include <picsim-debug.h>
#include <blitz/array.h>

namespace pdf {

#if !defined(DYNAMIC_CREATOR)
  using ranlib::Exponential;
  using ranlib::Normal;
#endif
  using ranlib::Uniform;
  using ranlib::UniformOpen;

  using std::ostream;

#define ID "$Id: pdf.cpp,v 1.82 2024/02/13 19:28:18 patrick Exp $"

#if 0
  template <typename T_num>
  ostream& operator<<(ostream& os, const NumberGenerator<T_num> &d)
  {
    return os << parser::header("Random generator setup")
           << "seed = " << d.seed << '\n';
  }
#endif

  template <typename T_num>
  NumberGenerator<T_num>::NumberGenerator(int nargs, char *args[], rngint s) :
    Parser(nargs, args),  seed(s)
#if defined(DYNAMIC_CREATOR) && defined(HAVE_MPI)
    , rng(rankProc)
#endif
  {
    initParsing (nargs, args);
    paramParsing();
  }

  template <typename T_num>
  NumberGenerator<T_num>::~NumberGenerator()
  {}

  template <typename T_num>
  void NumberGenerator<T_num>::initialise()
  {
#if defined(DYNAMIC_CREATOR)
// The number of independent generators is currently limited to 48
// The dirty trick to remedy is to use different seeding
#if defined(HAVE_MPI)
#if 1
    rngint recv, *send=0;
    if (rankProc == masterProc) {
      rng.seed(seed);
      send = new rngint[nbProc];
      for (int i=0; i<nbProc; ++i) {
        // UINT32_MAX = 4294967295U
        send[i] = rngint(4294967295U*rng.random());
      }
    }

    MPI_Scatter(send, 1, MPI_UNSIGNED, &recv, 1, MPI_UNSIGNED, mpiRoot, MPI_COMM_WORLD);

    BEGIN_DEBUG_OUTPUT(10)
    HEADER_DEBUG_OUTPUT1("NumberGenerator<T_num>::initialise()")
    VAR_DEBUG_OUTPUT1(recv)
    END_DEBUG_OUTPUT

    rng.seed(recv);
#else
    rng.seed(seed + rankProc);
#endif
#else // !defined(HAVE_MPI)
    rng.seed(seed);
#endif // defined(HAVE_MPI)
#else // !defined(DYNAMIC_CREATOR) 
    Uniform<T_num> rng;
#if defined(HAVE_MPI)
    rng.seed(seed + rankProc);
#else // !defined(HAVE_MPI)
    rng.seed(seed);
#endif // defined(HAVE_MPI)
#endif // defined(DYNAMIC_CREATOR)
  }

  template <typename T_num>
  T_num NumberGenerator<T_num>::normal(T_num m, T_num s)
  {
#if defined(DYNAMIC_CREATOR)
    // http://en.wikipedia.org/wiki/Box-Muller_transform
    // polar form of the Box-Muller transformation is both faster
    // and more robust numerically
    // rng should be in [0, 1]
    T_num u, v, w;
    do {
      u = 2.0 * rng.random() - 1.0;
      v = 2.0 * rng.random() - 1.0;
      w = u * u + v * v;
    } while ( w > 1.0 || w == 0.0 );

    return m + s * u * sqrt(-2.0 * log(static_cast<double>(w) ) / w);
#else
    Normal<T_num> rng(m,s);
    return rng.random();
#endif
  }

  template <typename T_num>
  T_num NumberGenerator<T_num>::uniform(T_num a, T_num b)
  {
    if (a > b)
      throw ClassException("NumberGenerator", "In uniform(a, b): a > b");

#if defined(DYNAMIC_CREATOR)
    // rng is in [0, 1)
    T_num x = (b-a)*rng.random()+a;
#else
    UniformOpen<T_num> rng;
    T_num x = (b-a)*rng.random()+a;
#endif
    // check whether the number is equal to the upper limit
    // float truncation possible error
    if (x >= b) {
      std::cout << "NumberGenerator::uniform(" << a << "," << b << ")=" << x
                << " (x-b)=" << (x-b);
      x = b*(1.0-blitz::epsilon(static_cast<T_num>(1.0)));
      std::cout << " --> fixed to "<< x << " (x-b)=" << (x-b) << std::endl;
    }
#if 0 && defined(HAVE_MPI)
    for (int ip=0; ip<nbProc; ++ip) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (ip == rankProc) {
        std::cout <<"rankProc " << rankProc << " rng " << x << std::endl;
      }
    }
#endif
    return x;
  }

  template <typename T_num>
  T_num NumberGenerator<T_num>::linear(T_num x1, T_num x2, T_num y1, T_num y2)
  {
    using blitz::pow2;

    if (x1 > x2)
      throw ClassException("NumberGenerator", "In linear(x1, x2, y1, y2): x1 > x2");

#if defined(DYNAMIC_CREATOR)
    // rng is in [0, 1)
    T_num y = rng.random();
#else
    UniformOpen<T_num> rng;
    T_num y = rng.random();
#endif

    // coefficients for line passing through (x1,y1) and (x2,y2)
    // y = ax + b where
    // a = (y2-y1)/(x2-x1)
    // b = (y1*x2-y2*x1)
    T_num a = (y2-y1)/(x2-x1);
    T_num b = (y1*x2-y2*x1);

    T_num x;
    if (a != 0.0) {
      // inverse cumulated probability distriubion
      T_num A = a/2;
      T_num B = b;
      T_num C = -a/2*(pow2(x1)+y*(pow2(x2)-pow2(x1))) - b*(x1+y*(x2-x1));

      T_num delta = pow2(B)-4*A*C;

      x = (-B+sqrt(delta))/(2*A);

    } else {
      x = (x2-x1)*y+x1;
    }

    return x;
  }

  template <typename T_num>
  T_num NumberGenerator<T_num>::rayleigh(T_num s)
  {
    // rng should be in [0, 1)
#if defined(DYNAMIC_CREATOR)
    return s * sqrt(-2.0*log(static_cast<double>(rng.random())));
#else
    UniformOpen<T_num> rng;
    return s * sqrt(-2.0*log(static_cast<double>(rng.random())));
#endif
  }

  template <typename T_num>
  T_num NumberGenerator<T_num>::exponential(T_num s)
  {
    // rng should be in [0, 1)
#if defined(DYNAMIC_CREATOR)
    return -1.0/s * log(static_cast<double>(rng.random()));
#else
    Exponential<T_num> rng(1.0/s);
    return rng.random();
#endif
  }

  template <typename T_num>
  double NumberGenerator<T_num>::F(double x, double m, double b)
  {
    // -exp(-x^2)+m/2/b*sqrt(2*pi)*erf(x)
    return -exp(-x*x)+m/b/M_2_SQRTPI*M_SQRT2*erf(x);
  }

  template <typename T_num>
  double NumberGenerator<T_num>::dF(double x, double m, double b)
  {
    // 2*x+m/b/sqrt(2)*exp(-x^2)
    return (2.0*x+m/b*M_SQRT2)*exp(-x*x);
  }

  template <typename T_num>
  double NumberGenerator<T_num>::fluxCDF(double x, double m, double b)
  {
    const double z = (x-m)/M_SQRT2/b;
    const double F0 = F(-m/M_SQRT2/b, m, b);
    const double Finf = m/b/M_2_SQRTPI*M_SQRT2;
    return (F(z,m,b) - F0)/(Finf - F0);
  }

  template <typename T_num>
  double NumberGenerator<T_num>::dfluxCDF(double x, double m, double b)
  {
    const double z = (x-m)/M_SQRT2/b;
    const double F0 = F(-m/M_SQRT2/b, m, b);
    const double Finf = m/b/M_2_SQRTPI*M_SQRT2;
    return dF(z,m,b)/(Finf - F0)/M_SQRT2/b;
  }

  template <typename T_num>
  T_num NumberGenerator<T_num>::flux(T_num m, T_num b)
  {
    if (m == 0.0)
      return rayleigh(b);

    double EPS = 1e-2;
    double dy;

#if defined(DYNAMIC_CREATOR)
    double x = static_cast<double>(rng.random());
#else
    UniformOpen<T_num> rng;
    double x = static_cast<double>(rng.random());
#endif
    // starting point is inflection point of CDF
    double y = 0.5*(m+sqrt(static_cast<double>(m*m+4.0*b*b)));
    // Newton iteration rule
    do {
      dy = (x-fluxCDF(y,static_cast<double>(m),static_cast<double>(b)))/
           dfluxCDF(y,static_cast<double>(m),static_cast<double>(b));
      y += dy;
#if defined(isfinite)

      if (!isfinite(y)) // Retry if NaN or Inf
        return flux(m, b);
#elif defined(finite)

      if (!finite(y)) // Retry if NaN or Inf
        return flux(m, b);
#endif

    } while (std::abs(dy) > EPS);

    return static_cast<T_num>(y);
  }

  template <typename T_num>
  double NumberGenerator<T_num>::perturbCDF(double x, double m, double s, double alpha)
  {
    const double erf0 = erf(m/M_SQRT2/s);
    const double erf1 = erf((1.0-m)/M_SQRT2/s);
    const double erfx = erf((x-m)/M_SQRT2/s);

    const double k = 2.0*alpha/(erf1 + erf0);
    const double beta = 1.0 + k*0.5*(erf1 + erf0);
    return 1.0/beta*(x + k*0.5*(erfx + erf0));
  }

  template <typename T_num>
  double NumberGenerator<T_num>::dperturbCDF(double x, double m, double s, double alpha)
  {
    const double erf0 = erf(m/M_SQRT2/s);
    const double erf1 = erf((1.0-m)/M_SQRT2/s);
    const double isqrt2pis = 0.5*M_2_SQRTPI*M_SQRT1_2/s;
    const double expx = exp(-(x-m)*(x-m)/2.0/(s*s));

    const double k = 2.0*alpha/(erf1 + erf0);
    const double beta = 1.0 + k*0.5*(erf1 + erf0);
    return 1.0/beta*(1.0 + k*isqrt2pis*expx);
  }

  template <typename T_num>
  T_num NumberGenerator<T_num>::perturb(T_num a, T_num b, T_num m, T_num s, T_num alpha)
  {
    if (alpha == 0.0)
      return uniform(a, b);

    const double erf0 = erf(m/M_SQRT2/s);
    const double erf1 = erf((1.0-m)/M_SQRT2/s);
    const double M_SQRTPI_2 = 1.0/M_2_SQRTPI*M_SQRT2;
    T_num amax = s*M_SQRTPI_2*(erf1 + erf0);
    if (std::abs(alpha) > amax && alpha<0.0) {
      ostringstream os;
      os << "|alpha| must be smaller than " << amax;
      throw ClassException("NumberGenerator", os.str());
    }

    if (a > b)
      throw ClassException("NumberGenerator", "In uniform(a, b): a > b");

    double EPS = 1e-2;
    double dy;

#if defined(DYNAMIC_CREATOR)
    double x0 = static_cast<double>(rng.random());
#else
    UniformOpen<T_num> rng;
    double x0 = static_cast<double>(rng.random());
#endif
    // starting point is inflection point of CDF
    double y = static_cast<double>(m);
    // Newton iteration rule
    do {
      dy = (x0-perturbCDF(y,m,s,alpha))/dperturbCDF(y,m,s,alpha);
      y += dy;
    } while (std::abs(dy) > EPS);

    T_num x = (b-a)*y+a;
    // check whether the number is equal to the upper limit
    // float truncation possible error
    if (x >= b) {
      x = b*(1.0-blitz::epsilon(static_cast<T_num>(1.0)));
    }
    return x;
  }

  template <typename T_num>
  void NumberGenerator<T_num>::initParsing(int nargs, char *args[])
  {
    registerClass("NumberGenerator");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PDF_COPYRIGHT);

    using parser::types::integer;

    parseLevelDebugOption("NumberGenerator::dl");

    insertOption(_seed, "seed", integer, "Seed for the random generator", Any(seed));
    insertOptionAlias(_seed,"-seed");
  }

  template <typename T_num>
  void NumberGenerator<T_num>::paramParsing()
  {
    parseOption(_seed, seed);
  }


  template class NumberGenerator<float>;
  template class NumberGenerator<double>;

#if 0
  template ostream& operator<<(ostream& os, const NumberGenerator<float> &d);
  template ostream& operator<<(ostream& os, const NumberGenerator<double> &d);
#endif


#undef ID

}
