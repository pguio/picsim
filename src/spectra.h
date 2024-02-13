/**************************************************************************
 *
 * $Id: spectra.h,v 1.53 2011/03/26 15:36:08 patrick Exp $
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

#ifndef SPECTRA_H
#define SPECTRA_H

#include <scheduler.h>
#include <fourier.h>

namespace picsim {

  class Spectra : public parser::Parser {
  public:

    typedef std::ostream ostream;
    typedef std::ostringstream ostringstream;
    typedef std::string string;

    typedef blitz::Array<int    , 1> ivector;

    typedef blitz::Array<real  , 1> rvector;
    typedef blitz::Array<real  , 2> rmatrix;
    typedef blitz::Array<real  , 3> rtensor;

    typedef std::complex<double> complex;
    typedef blitz::Array<complex, 1> cvector;
    typedef blitz::Array<complex, 2> cmatrix;
    typedef blitz::Array<complex, 3> ctensor;

    typedef blitz::Range Range;

    typedef mudfas::Field Field;

    friend ostream& operator<<(ostream& os, const Spectra &s);

    Spectra(int nargs, char *args[], const string name="");
    ~Spectra();

    void initialise(Scheduler &scheduler);

    const string &name() const {
      return spectraName;
    }

    void update(Field &phi, Scheduler &scheduler);
    void getSpectra(int &it, rvector &Phik, rmatrix &Phikw);

    bool isSpectraComplete(const Scheduler &scheduler);

    rvector getk()    const {
      return k;
    }
    rvector getkc()   const {
      return kc;
    }
    rvector getSk()   const {
      return Sk;
    }
    rvector getSfrq() const {
      return Sfrq;
    }

    rvector getTimeAxis(const Scheduler &scheduler);

  protected:

    typedef parser::Parser Parser;

  private:

    enum parser_enum { _iterSk=1, _iterSkf, avSkw2nfft,
                       avSkw2kindex, avSkw2kmin, avSkw2kmax, avSkw2nbk,
                       avSrmin, avSrmax
                     };

    const string spectraName;

    int nk;
    double L, tau, fs;
    RangeSpec<int> iterSk, iterSkf;
    int spec_nfft;
    ivector k_index;
    int nb_ks;
    double kmin, kmax;
    Boundary<int,DIMR> av_rlim;
    int nrint;

    rvector     k, kc;
    rvector     av_phik2;
    int         av_number;

    ctensor phikt;
    rtensor phikw;
    rvector Sk;
    rvector Sfrq;
    int spec_number;

    int time_index;

    fourier::ODFT1D<complex,complex> fftTime;
    fourier::ODFT1D<complex,complex> fftSpace;

    ivector signed_index;

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
    void allocateArrays();
  };


  inline
  bool Spectra::isSpectraComplete(const Scheduler &scheduler)
  {
    int citer(scheduler.getCurrentIter());
    return citer >= iterSk.start() && citer <= iterSk.end() &&
           (((citer-iterSkf.start())/iterSkf.stride()) % spec_nfft) == spec_nfft-1;
  }

}

#endif // SPECTRA_H

