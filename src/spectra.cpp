/**************************************************************************
 *
 * $Id: spectra.cpp,v 1.74 2011/03/26 15:36:08 patrick Exp $
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

#include <spectra.h>

namespace picsim {

  using blitz::imag;
  using blitz::pow2;
  using blitz::sum;
  using blitz::zip;

  using std::cout;
  using std::endl;
  using std::ostream;
  using std::setw;

#define ID "$Id: spectra.cpp,v 1.74 2011/03/26 15:36:08 patrick Exp $"

  ostream& operator<<(ostream& os, const Spectra &s)
  {
    os   << "Spectra name                  = " << s.spectraName << '\n'
         << "Spatial period                = " << s.L << '\n'
         << "Grid size configuration space = " << s.nk << '\n'
         << "Time increment                = " << s.tau << '\n'
         << "<|S(k)|^2>   time spec        = " << s.iterSk << '\n'
         << "<|S(k,f)|^2> time spec        = " << s.iterSkf << '\n'
         << "<|S(k,f)|^2> sampling freq    = " << s.fs << '\n'
         << "<|S(k,f)|^2> number of pts    = " << s.spec_nfft << '\n'
         << "Number of spatial integration = " << s.nrint << '\n'
         << "Number of wavevectors k       = " << s.nb_ks << '\n';

    if (s.nb_ks != 0) {
      os << "Index and value of the k's      =\n";
      for (int ik=0; ik<s.nb_ks-1; ++ik) {
        os << "\tk("  << setw(2) << s.k_index(ik) << ") = "
           << s.k(s.k_index(ik)) << '\n';
      }
      os << "\tk(" << setw(2) << s.k_index(s.nb_ks-1) << ") = "
         << s.k(s.k_index(s.nb_ks-1)) << '\n';
    }
    return os;
  }


  Spectra::Spectra(int nargs, char *args[], const string name) :
    Parser(nargs, args), spectraName(name),
    nk(0), L(0.0), tau(0.0), iterSk(-1), iterSkf(-1),
    spec_nfft(0), k_index(0), nb_ks(0), kmin(0.0), kmax(0.0),
    av_rlim(-1,-1), nrint(0),
    k(0), kc(0), av_phik2(0), av_number(0), phikt(0),
    phikw(0), Sk(0), Sfrq(0), spec_number(0), time_index(0)
  {
    initParsing(nargs, args);
    paramParsing();
    checkParam();
  }

  Spectra::~Spectra()
  {}


  void Spectra::initialise(Scheduler &scheduler)
  {
    using blitz::tensor::i;

    RVectori Size(scheduler.getGridSize());
    RVectori Avmin(av_rlim.min);
    RVectori Avmax(av_rlim.max);

    // find out spatial Fourier transform dimension
    for (int d=0; d<DIMR; ++d) {
      if (Avmin(d) == -1 && Avmax(d) == -1) {
        nk = Size(d);
      }
    }

    Boundary<real,DIMR> Domain(scheduler.getDomainBoundary());
    RVectorr rmin(Domain.min);
    RVectorr rmax(Domain.max);
    // find out spatial length
    for (int d=0; d<DIMR; ++d) {
      if (Avmin(d) == -1 && Avmax(d) == -1) {
        L = rmax(d)-rmin(d);
      }
    }

    // sampling period and frequency
    tau = static_cast<double>(scheduler.getTimeInc());
    fs = 1.0/(iterSkf.stride()*tau);

    // find out how many spatial profile to be integrated
    nrint = 1;
    for (int d=0; d<DIMR; ++d) {
      if (Avmin(d)!=-1 && Avmax(d)!=-1) {
        nrint *= (Avmax(d)-Avmin(d)+1);
      }
    }

    allocateArrays();

    // fftw plan for spatial fourier transform (-1)
    fftSpace.resize(nk, -1);

    // symmetric indices [-N/2:1:N/2-1] (even N),  [-N/2:1:N/2] (odd N)
    signed_index = i-(nk>>1);

    // Grid points k in fourier space k=0 centered
    kc = 2.0*M_PI*signed_index/L;
    k = kc;
    fourier::ifftshift(k);

    av_phik2 = 0.0;

    fftSpace.resize(spec_nfft, -1);

    for (int ik=0; ik<nb_ks; ++ik) {
      Sk(ik) = k(k_index(ik));
    }
    cout << "Sk=" << Sk << endl;

    Sfrq = fs*(i-(spec_nfft>>1))/((spec_nfft>>1)<<1);

    phikt = 0.0;
    phikw = 0.0;
  }

  void Spectra::update(Field &phi, Scheduler &scheduler)
  {
    int citer(scheduler.getCurrentIter());
    RVectori Avmin(av_rlim.min);
    RVectori Avmax(av_rlim.max);

    using blitz::real;

    cvector phir(nk);
    cvector phik(nk);
    if (citer >= iterSk.start() && citer <= iterSk.end()) {
      int is(0);
#if (DIMR==2)

      for (int i=Avmin(0); i<=Avmax(0); ++i) {
        for (int j=Avmin(1); j<=Avmax(1); ++j) {
          if (Avmin(0) == -1 && Avmax(0) == -1) {
            phir = zip(phi(Range::all(),j), 0.0, complex());
          }
          if (Avmin(1) == -1 && Avmax(1) == -1) {
            phir = zip(phi(i,Range::all()), 0.0, complex());
          }
          fftSpace.direct(phir, phik);
          av_phik2 += pow2(real(phik))+
                      pow2(imag(phik));
          ++av_number;

          int it(((citer-iterSkf.start())/iterSkf.stride()) % spec_nfft);
          for (int ik=0; ik<nb_ks; ++ik) {
            phikt(is,ik,it) = phik(k_index(ik));
          }
          if (it == spec_nfft-1) { // if time sequence completed
            cvector phi1kt(spec_nfft);
            cvector phi1kw(spec_nfft);
            for (int ik=0; ik<nb_ks; ++ik) {
              // <|phi(k,t)|^2> to <|phi(k,\omega)|^2>
              phi1kt=phikt(is,ik,Range::all());
              fftTime.direct(phi1kt, phi1kw);
              phikw(is,ik,Range::all()) += real(phi1kw)*real(phi1kw)+
                                           imag(phi1kw)*imag(phi1kw);
            }
          }
          ++spec_number;
          ++is;
        }
      }
#elif (DIMR==3)
      for (int i=Avmin(0); i<=Avmax(0); ++i) {
        for (int j=Avmin(1); j<=Avmax(1); ++j) {
          for (int k=Avmin(2); k<=Avmax(2); ++k) {
            if (Avmin(0) == -1 && Avmax(0) == -1) {
              phir = zip(phi(Range::all(), j, k), 0.0, complex());
            }
            if (Avmin(1) == -1 && Avmax(1) == -1) {
              phir = zip(phi(i, Range::all(), k), 0.0, complex());
            }
            if (Avmin(2) == -1 && Avmax(2) == -1) {
              phir = zip(phi(i, j, Range::all()), 0.0, complex());
            }
            fftSpace.direct(phir, phik);
            av_phik2 += real(phik)*real(phik)+imag(phik)*imag(phik);
            ++av_number;

            int it(((citer-iterSkf.start())/iterSkf.stride()) % spec_nfft);
            for (int ik=0; ik<nb_ks; ++ik) {
              phikt(is,ik,it) = phik(k_index(ik));
            }
            if (it == spec_nfft-1) { // if time sequence completed
              cvector phi1kt(spec_nfft);
              cvector phi1kw(spec_nfft);
              for (int ik=0; ik<nb_ks; ++ik) {
                // <|phi(k,t)|^2> to <|n(k,\omega)|^2>
                phi1kt=phikt(is,ik,Range::all());
                fftTime.direct(phi1kt, phi1kw);
                fourier::fftshift(phi1kw);
                phikw(is,ik,Range::all()) += real(phi1kw)*real(phi1kw)+
                                             imag(phi1kw)*imag(phi1kw);
              }
            }
            ++spec_number;
            ++is;
          }
        }
      }
#endif

    }
  }

  void Spectra::getSpectra(int &it, rvector &Phik, rmatrix &Phikw)
  {
    using blitz::firstIndex;
    using blitz::secondIndex;
    using blitz::thirdIndex;

    Phik.resize(nk);
    // divide by av_number AND NFFT^2 due to the definition of the FFT
    Phik = av_phik2 * (1.0/(nk*nk)/av_number);
    fourier::fftshift(Phik);

    Phikw.resize(nb_ks, spec_nfft);

    // do re-ordering manually and take the average over space
    firstIndex i;
    secondIndex j;
    thirdIndex k;
    Phikw = sum(phikw(k,i,j), k) * (1.0/spec_number);

    rvector phi1kw(spec_nfft);
    for (int ik=0; ik<nb_ks; ++ik) {
      phi1kw=Phikw(ik,Range::all());
      fourier::fftshift(phi1kw);
      Phikw(ik,Range::all())=phi1kw;
    }
    it = time_index++;

    av_phik2 = 0.0;
    av_number = 0;

    phikt = 0.0;
    phikw = 0.0;
    spec_number = 0;
  }

  Spectra::rvector Spectra::getTimeAxis(const Scheduler &scheduler)
  {
    rvector Stimes(0);

    if (iterSk.start() < 0 && iterSk.end() < 0)
      return Stimes;

    int ntimes((scheduler.getMaxIter()-iterSkf.start()+1)/
               (iterSkf.stride()*spec_nfft));

    Stimes.resize(ntimes);
    for (int i=0; i<ntimes; ++i)
      Stimes(i) = tau*(iterSkf.start()+i*spec_nfft*iterSkf.stride());

    cout << "Stimes = " << Stimes << endl;

    return Stimes;
  }

  void Spectra::initParsing(int nargs, char *args[])
  {
    registerClass("Spectra");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", PICSIM_COPYRIGHT);

#if 0
    const char ivect[] = "int[3]";
    const char rveci[] = "int["  DIMRSTR "]";
#endif

    using parser::types::integer;
    using parser::types::real;
    using parser::types::intVect;

    insertOption(_iterSk     , "iterSk"      ,intVect, "<|S(k)|^2>   iter spec"        , Any(iterSk.data));

    insertOption(_iterSkf    , "iterSkf"     ,intVect, "<|S(k,f)|^2> iter spec"        , Any(iterSkf.data));

    insertOption(avSkw2nfft  , "avSkw2nfft"  ,integer, "Number of points <|S(k,f)|^2>" , Any(spec_nfft));
    insertOption(avSkw2kindex, "avSkw2kindex",intVect, "Indices of the k"              , Any());

    insertOption(avSkw2kmin  , "avSkw2kmin"  ,real   , "k min for <|S(k,f)|^2>"        , Any(kmin));
    insertOption(avSkw2kmax  , "avSkw2kmax"  ,real   , "k max for <|S(k,f)|^2>"        , Any(kmax));
    insertOption(avSkw2nbk   , "avSkw2nbk"   ,integer, "number of k's for <|S(k,f)|^2>", Any(nb_ks));


    insertOption(avSrmin     , "avSrmin"     ,intVect, "LB space integration"          , Any(av_rlim.min));
    insertOption(avSrmax     , "avSrmax"     ,intVect, "UB space integration"          , Any(av_rlim.max));
  }

  void Spectra::paramParsing()
  {
    parseOption(_iterSk , iterSk.data);
    parseOption(_iterSkf, iterSkf.data);

    parseOption(avSkw2nfft, spec_nfft);
    VecInt kindex(0);
    if (parseOption(avSkw2kindex, kindex)) {
      nb_ks = kindex.size();
      kindex.resize(nb_ks);
      for (int ik=0; ik<nb_ks; ++ik) {
        k_index(ik) = kindex[ik];
      }
    }
    parseOption(avSkw2kmin, kmin);
    parseOption(avSkw2kmax, kmax);
    if (parseOption(avSkw2nbk, nb_ks)) {
      k_index.resize(nb_ks);
    }

    parseOption(avSrmin, av_rlim.min);
    parseOption(avSrmax, av_rlim.max);
  }

  void Spectra::checkParam() const
  {
    if (iterSk.start() < 0 || iterSk.end() < 0 ||
        iterSkf.start() < 0 || iterSkf.end() ) {
      throw ClassException("Spectra","iterSk and iterSk must be initialised");
    }

    if (nb_ks>0 && kmin==0.0 && kmax==0.0) {
      for (int ik=0; ik<nb_ks; ++ik) {
        if (k_index(ik) < 0 || k_index(ik) >= (nk>>1)) {
          ostringstream os;
          os <<"avSkw2kindex=" << k_index(ik) <<
             " must be smaller than nk/2 " << (nk>>1);
          throw ClassException("Spectra", os.str());
        }
      }
      int diff(k_index(1)-k_index(0));
      for (int ik=1; ik<nb_ks-1; ++ik) {
        if (k_index(ik+1)-k_index(ik) != diff ) {
          throw ClassException("Spectra",
                               "avSkw2kindex must be regularly spaced");
        }
      }
    }
    if (nb_ks>0 && kmin!=0.0 && kmax!=0.0) {
      if (nb_ks==1) {
        throw ClassException("Spectra",
                             "avSkw2nbk must be strictly larger than 1");
      }
      if (kmin>=kmax) {
        throw ClassException("Spectra",
                             "avSkw2kmax must be strictly larger than avSkw2kmin");
      }
    }
  }

  void Spectra::allocateArrays()
  {
    k.resize(nk);
    kc.resize(nk);
    signed_index.resize(nk);
    av_phik2.resize(nk);

    phikt.resize(nrint, nb_ks, spec_nfft);
    phikw.resize(nrint, nb_ks, spec_nfft);
    Sk.resize(nb_ks);
    Sfrq.resize(spec_nfft);
  }


#undef ID

}
