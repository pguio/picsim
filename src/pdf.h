/**************************************************************************
 *
 * $Id: pdf.h,v 1.51 2016/06/01 15:03:03 patrick Exp $
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

#ifndef PDF_H
#define PDF_H

#include <classexception.h>
#include <parser.h>

#include <blitz/numinquire.h>

#define DYNAMIC_CREATOR 1

#include <random/uniform.h>

#if !defined(DYNAMIC_CREATOR)
#include <random/normal.h>
#include <random/exponential.h>
#endif

namespace pdf {


#define PDF_COPYRIGHT \
	  "Copyright (c) 2000-2005 Patrick Guio <patrick.guio@gmail.com>\n\n"\
	  "This is free software; see the source for copying conditions.  There is NO\n"\
	  "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"


  template <typename T_num>
  class NumberGenerator : public parser::Parser {
  public:

    typedef std::ostream ostream;
    typedef std::ostringstream ostringstream;
    typedef ranlib::IRNG_int rngint;

    friend ostream& operator<<(ostream& os, const NumberGenerator<T_num> &d) {
      return os << parser::header("Random generator setup")
             << "seed = " << d.seed << '\n';
    }


    NumberGenerator(int nargs=0, char *args[]=0, rngint s=1);
    ~NumberGenerator();

    void initialise();

    T_num normal(T_num m=0.0, T_num s=1.0);
    T_num uniform(T_num a=0.0, T_num b=1.0);
    T_num linear(T_num x1=0.0, T_num x2=1.0, T_num y1=0.5, T_num y2=1.5);
    T_num rayleigh(T_num s=1.0);
    T_num exponential(T_num s=1.0);
    T_num flux(T_num m=1.0, T_num b=1.0);
    T_num perturb(T_num a=0.0, T_num b=1.0, T_num m=0.5, T_num s=0.1,
                  T_num alpha=0.1);

    rngint getSeed() {
      return seed;
    }

  protected:

    typedef parser::Parser Parser;

  private:

    enum { _seed=1 };

    rngint seed;

    void initParsing(int nargs, char *args[]);
    void paramParsing();

    double F(double x, double m, double b);
    double dF(double x, double m, double b);
    double fluxCDF(double x, double m, double b);
    double dfluxCDF(double x, double m, double b);

    double perturbCDF(double x, double m, double s, double alpha);
    double dperturbCDF(double x, double m, double s, double alpha);

#if defined(DYNAMIC_CREATOR)
    typedef ranlib::MersenneTwister MersenneTwister;
    typedef ranlib::independentState independentState;
    // uniform random numbers in [0,1) (fastest one)
    ranlib::UniformClosedOpen<T_num, MersenneTwister, independentState> rng;
#endif
  };

}

#endif // PDF_H
