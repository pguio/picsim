#!/bin/sh

case `hostname` in

phalanx*)
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib64/hdf --with-fftw3=/usr/include,/usr/lib64 --with-blitz=$HOME/include,$HOME/lib64 --enable-mpi --disable-cxx-flags-preset MPICXX="/software/mpich2/bin/mpicxx" CXX=icpc CXXFLAGS="-xSSE4.1 -O3 -ip -restrict -vec-report1 -no-prec-div -no-ansi-alias" 
  ./configure --enable-cxx-flags-preset --enable-mpi --with-hdf4=/usr/include/hdf,/usr/lib64/hdf --with-fftw3=/usr/include,/usr/lib64 --with-blitz=$HOME/include,$HOME/lib64 --enable-optimize --disable-cxx-flags-preset MPICXX="icpc" CXX="icpc" CXXFLAGS="-xSSE4.1 -O3 -ip -restrict -vec-report1 -no-prec-div -no-ansi-alias" LDFLAGS="-L/usr/lib64" LIBS="-lmpi" AR=xiar
	;;

alun*|dare*)
#	./configure --with-fftw3=/usr/include,/usr/lib64 --with-hdf4=/usr/include/hdf,/usr/lib64/hdf --with-blitz=$HOME --enable-cxx-flags-preset --enable-mpi --enable-optimize MPICXX="/usr/lib64/openmpi/1.3.2-gcc/bin/mpiCC"
	./configure --with-fftw3=/usr/include,/usr/lib64 --with-hdf4=/usr/include/hdf,/usr/lib64/hdf --with-blitz=$HOME --enable-cxx-flags-preset --enable-mpi --enable-optimize MPICXX="/usr/lib64/mpich2/bin/mpic++"
        ;;
       

kryten*|hawke*)
#	./configure --with-fftw3=$HOME --with-hdf4=$HOME --with-blitz=$HOME --disable-cxx-flags-preset CXX="g++" CXXFLAGS="-pedantic -Wall -O -g -DBZ_DEBUG"
#	./configure --with-fftw3=$HOME --with-hdf4=$HOME --with-blitz=$HOME --enable-cxx-flags-preset --enable-optimize CXX="g++"
	./configure --with-fftw3=$HOME --with-hdf4=$HOME --with-blitz=$HOME --enable-cxx-flags-preset --enable-mpi --enable-optimize MPICXX="/usr/share/openmpi/1.2.3-gcc/bin64/mpiCC"
#	./configure --with-fftw3=$HOME --with-hdf4=$HOME --with-blitz=$HOME --disable-cxx-flags-preset --enable-mpi --enable-optimize MPICXX="/usr/share/openmpi/1.2.3-gcc/bin64/mpiCC" CXXFLAGS="-O3 -funroll-loops -fstrict-aliasing -fomit-frame-pointer -ffast-math"
        ;;

login*)
        if test "$LOGNAME" = "ucappgu"; then
	./configure --with-fftw3=/cm/shared/apps/fftw/intel/double/3.2.2 --with-hdf4=/shared/ucl/apps/HDF4.2r4 --with-blitz=$HOME --disable-cxx-flags-preset --enable-mpi --enable-optimize CXX="icpc" MPICXX="mpiCC" CXXFLAGS="-xSSSE3 -O3 -ipo -restrict -vec-report1 -no-prec-div -no-ansi-alias" AR=xiar
        elif test "$LOGNAME" = "patrickg"; then
#	./configure --with-fftw3=$HOME/include,$HOME/lib/gcc --with-hdf4=$HOME/include,$HOME/lib/gcc --with-blitz=$HOME/include,$HOME/lib/gcc --disable-cxx-flags-preset --enable-mpi --enable-optimize CXX="g++" MPICXX="mpic++" CXXFLAGS="-O2"
#	./configure --with-fftw3=$HOME/include,$HOME/lib/gcc --with-hdf4=$HOME/include,$HOME/lib/gcc --with-blitz=$HOME/include,$HOME/lib/gcc --disable-cxx-flags-preset --enable-mpi --enable-optimize CXX="pathCC" MPICXX="mpic++ -ccl pathCC" CXXFLAGS="-Ofast"
	./configure --with-fftw3=$HOME/include,$HOME/lib/gcc --with-hdf4=$HOME/include,$HOME/lib/gcc --with-blitz=$HOME/include,$HOME/lib/gcc --disable-cxx-flags-preset --enable-mpi --enable-optimize --enable-threadsafe-blitz CXX="icc -Kc++" MPICXX="mpic++ -ccl icc -Kc++" CXXFLAGS="-O3 -ip -static -no-prec-div -ansi-alias -xW -axW"
#	./configure --with-fftw3=$HOME/include,$HOME/lib/gcc --with-hdf4=$HOME/include,$HOME/lib/gcc --with-blitz=$HOME/include,$HOME/lib/gcc --disable-cxx-flags-preset --enable-mpi --enable-optimize CXX="pgCC" MPICXX="mpic++ -ccl pgCC" CXXFLAGS="-fast" 
  fi
	;;

eukalyptus*)
#	./configure --with-hdf4=/net/pl/hdf --with-fftw3=/net/pl/fftw --with-blitz=/net/pl/blitz --enable-optimize CXX=g++
	./configure --with-hdf4=/net/pl/hdf --with-fftw3=/net/pl/fftw --with-blitz=/net/pl/blitz --enable-mpi --enable-optimize MPICXX=mpiCC
	;;

fimm*)
        ./configure  --with-blitz=$HOME --with-fftw3=/local/fftw --with-hdf4=/local/HDF --enable-optimize --enable-mpi MPICXX=mpicxx
#        ./configure  --with-blitz=$HOME --with-fftw3=/local/fftw --with-hdf4=/local/HDF --enable-optimize CXX=icpc
        ;;

jazz6*)
./configure CONFIG_SITE=`pwd`/config/sites/f38_x86_64_gcc+openmpi_optim_jazz6
;;

jazz5*)
./configure CONFIG_SITE=`pwd`/config/sites/f32_x86_64_gcc_optim_jazz5
;;

jazz2*)
./configure CONFIG_SITE=`pwd`/config/sites/f31_x86_64_gcc_optim_jazz3
# cross-compiling for alun kryten dare hawke
#./configure --host=x86_64 --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --enable-mpi --disable-cxx-flags-preset MPICXX="mpic++" CXX=icpc CXXFLAGS="-m64 -xSSE2 -vec-report1 -O3 -ip -no-prec-div -ansi-alias -restrict -Wl,-rpath=/usr/lib/openmpi/lib,-rpath=/home/patrick/lib"
# ./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --disable-mpi --disable-cxx-flags-preset CXX=icpc CXXFLAGS="-g -O -DBZ_DEBUG"
#  ./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --disable-mpi --disable-cxx-flags-preset CXX=g++ CXXFLAGS="-pedantic -Wall -g -DBZ_DEBUG"
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local/blitz1 --disable-mpi --disable-cxx-flags-preset CXX=icpc CXXFLAGS="-xSSSE3 -vec-report1 -O3 -no-prec-div -ansi-alias -restrict"
	#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local/blitz2 --disable-mpi --disable-cxx-flags-preset CXX=icpc CXXFLAGS="-xSSSE3 -vec-report1 -O3 -no-prec-div -ansi-alias -restrict"
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --disable-mpi --disable-cxx-flags-preset CXX=icpc CXXFLAGS="-xSSSE3 -vec-report1 -O3 -ip -no-prec-div -ansi-alias -restrict"
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --disable-mpi --disable-cxx-flags-preset CXX=icpc CXXFLAGS="-g -Weffc++ -xSSSE3 -vec-report1 -O3 -ip -no-prec-div -ansi-alias -restrict"
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --enable-mpi --enable-optimize CXX=g++
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --disable-cxx-flags-preset --enable-optimize CXX=g++ CXXFLAGS="-O3 -mtune=core2"
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --disable-cxx-flags-preset  CXX=g++ CXXFLAGS="-Wall -pedantic"
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --enable-mpi --disable-cxx-flags-preset MPICXX="mpic++" CXX=g++ CXXFLAGS="-O3 -mtune=core2 -Wall -std=c++0x -pedantic"
# setenv MPICH2_MPICXX icpc
# setenv MPICH2_MPICXX_FLAGS "-xSSSE3 -O3 -ip -restrict -vec-report1 -no-prec-div -no-ansi-alias"
#setenv OMPI_CXXFLAGS "-xSSE4.2 -ansi -std=c++0x -O3 -ipo -restrict -vec-report1 -no-prec-div -no-ansi-alias"
# setenv MPICH2_MPICXX_FLAGS "-xSSE4.2 -ansi -std=c++0x -O3 -ip -restrict -vec-report1 -no-prec-div -no-ansi-alias"
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local/blitz,/usr/local/blitz/lib --enable-threadsafe-blitz --enable-mpi --disable-cxx-flags-preset MPICXX="mpicxx" CXX=icpc CXXFLAGS="-xSSE4.2 -ansi -std=c++0x -O3 -ip -restrict -vec-report1 -no-prec-div -no-ansi-alias"  AR=xiar
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --enable-threadsafe-blitz --enable-mpi --disable-cxx-flags-preset MPICXX="mpicxx" CXX=icpc CXXFLAGS="-xSSE4.2 -ansi -std=c++0x -O3 -ipo -restrict -vec-report1 -no-prec-div -no-ansi-alias"  AR=xiar
#	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local --enable-mpi --disable-cxx-flags-preset MPICXX="mpicxx" CXX=icpc CXXFLAGS="-xSSSE3 -O3 -ip -restrict -vec-report1 -no-prec-div -no-ansi-alias"
	;;

jazz*|fyspc-rp42*|theta*)
	#./configure --with-blitz=/usr/local --with-hdf4=/usr/local/hdf --with-fftw3=/usr --enable-mpi --enable-dependency-tracking CXX=g++
	#./configure --with-blitz=/usr/local --with-hdf4=/usr/local/hdf --with-fftw=/usr/local/fftw --enable-mpi --enable-dependency-tracking CXX=g++
	#./configure --with-blitz=/usr/local --with-hdf4=/usr/local/hdf --with-fftw=/usr/local/fftw --enable-dependency-tracking CXX=g++
	./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --disable-cxx-flags-preset --enable-debug CXX=g++ CXXFLAGS="-pedantic -Wall -g -DBZ_DEBUG"
	#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz=/usr/local CXX=icpc 
	#./configure --with-hdf4=/usr/include/hdf,/usr/lib/hdf --with-fftw3 --with-blitz --disable-cxx-flags-preset CXX=icpc CXXFLAGS="-xN -O3 -ip -no-prec-div -ansi-alias"
	#./configure --with-hdf4=/usr/local/hdf --with-fftw3 --enable-mpi MPICXX=mpiCC
	#./configure --with-blitz=/usr/local --with-hdf4=/usr/local/hdf --with-fftw=/usr/local/fftw --enable-dependency-tracking CXX=g++
	#./configure --with-blitz=/usr/local --with-hdf4=/usr/local/hdf --with-fftw=/usr/local/fftw --enable-dependency-tracking CXX=icpc
	#./configure --with-blitz=/usr/local --with-hdf4=/usr/local/hdf --with-fftw=/usr/local/fftw --enable-double-field --enable-dependency-tracking CXX=g++
	#./configure --with-blitz=/usr/local --with-hdf4=/usr/local/hdf --with-fftw=/usr/local/fftw --enable-double-coordinate --enable-dependency-tracking CXX=g++
	;;

login*)
  #./configure --with-fftw3=$HOME --with-hdf4=$HOME --with-blitz=$HOME --enable-cxx-flags-preset --enable-optimize CXX="g++"
	 ./configure --with-fftw3=$HOME --with-hdf4=$HOME --with-blitz=$HOME --enable-cxx-flags-preset --enable-mpi --enable-optimize MPICXX="mpiCC"
  ;;

phact*|hyades*|ayil*|barakish*|menkar*|markab*|papsukal*|thuban*)
	if test -z "$1"; then 
#	./configure --with-blitz=$HOME --with-fftw3=$HOME --with-hdf4=$HOME/src/hdf --enable-dependency-tracking --enable-optimize CC=cc CXX=cxx
	./configure --with-blitz=$HOME --with-fftw3=$HOME --with-hdf4=$HOME/src/hdf --enable-dependency-tracking --disable-cxx-flags-preset --enable-optimize CC=cc CXX=cxx CXXFLAGS="-std ansi -D__USE_STD_IOSTREAM -DBZ_ENABLE_XOPEN_SOURCE -D_OSF_SOURCE -ieee -model ansi -accept restrict_keyword -nousing_std -no_implicit_include" CXX_OPTIMIZE_FLAGS="-fast -inline speed -nocleanup"
	#./configure --with-blitz=$HOME --with-fftw3=$HOME --with-hdf4=$HOME/src/hdf --enable-dependency-tracking --enable-optimize --enable-debug CC=cc CXX=cxx 
	elif test "$1"=mpi; then
#	./configure --with-blitz=$HOME --with-fftw3=$HOME --with-hdf4=$HOME/src/hdf --enable-dependency-tracking --enable-optimize --enable-mpi CC=cc CXX=cxx
	./configure --with-blitz=$HOME --with-fftw3=$HOME --with-hdf4=$HOME/src/hdf --enable-dependency-tracking --disable-cxx-flags-preset --enable-optimize --enable-mpi CC=cc CXX=cxx CXXFLAGS="-std ansi -D__USE_STD_IOSTREAM -DBZ_ENABLE_XOPEN_SOURCE -D_OSF_SOURCE -ieee -model ansi -accept restrict_keyword -nousing_std -no_implicit_include" CXX_OPTIMIZE_FLAGS="-fast -inline speed -nocleanup"
	fi
	;;

tre*)
	if test -z "$1"; then
	./configure --with-blitz=$HOME --with-fftw=/usr/local/fftw-2.1.5 --with-hdf4=/usr/local/HDF4 --enable-optimize CXX=xlC
	elif test "$1"=mpi; then
	#./configure --with-blitz=$HOME --with-fftw=/usr/local/fftw-2.1.5 --with-hdf4=/usr/local/HDF4 --enable-mpi --enable-optimize CXX=xlC
	./configure --enable-64bit --with-blitz=$HOME --with-fftw=/usr/local/fftw-2.1.5-64 --with-hdf4=/usr/local/HDF4-64 --enable-mpi --enable-optimize CC=gcc CFLAGS=-maix64 CXX=xlC
	fi
	;;

gridur*|embla*|balder*)
	if test -z "$1"; then
	./configure --enable-64bit --with-blitz=$HOME --with-fftw=$HOME --with-hdf4=$HOME -enable-optimize CC=cc CFLAGS=-64 CXX=CC
	elif test "$1"=mpi; then
	./configure --enable-64bit --with-blitz=$HOME --with-fftw=$HOME --with-hdf4=$HOME -enable-optimize -enable-mpi CC=cc CFLAGS=-64 CXX=CC
	fi
	;;
	
magnum*| pico*)
	if test -z "$1"; then
	./configure --with-blitz=$HOME --with-fftw=/site/fftw --with-hdf4=/site/hdf --enable-dependency-tracking --enable-optimize CC=cc CXX=aCC
	elif test "$1"=mpi; then
	./configure --with-blitz=$HOME --with-fftw=/site/fftw --with-hdf4=/site/hdf --enable-dependency-tracking --enable-optimize --enable-mpi CC=cc CXX=aCC
	fi
	;;

sp0201) # CINECA IBM SP6 sp.sp6.cineca.it
        #module load autoload fftw
        #module load autoload hdf4
	if test -z "$1"; then
	#./configure --enable-64bit --with-blitz=$HOME --with-fftw3=$FFTW_INC,$FFTW_LIB --with-hdf4=$HDF4_HOME --enable-optimize CC=xlc CXX=xlC LDFLAGS="-L$SZLIB_LIB -L$ZLIB_LIB"
	./configure --enable-64bit --with-blitz=$HOME --with-fftw3=$FFTW_INC,$FFTW_LIB --with-hdf4=$HDF4_HOME --disable-cxx-flags-preset CC=xlc CXX=xlC CXXFLAGS="-O0 -g -C" LDFLAGS="-L$SZLIB_LIB -L$ZLIB_LIB"
	#./configure --enable-64bit --with-blitz=$HOME --with-fftw3=$FFTW_INC,$FFTW_LIB --with-hdf4=$HDF4_HOME --disable-cxx-flags-preset CC=xlc CXX=xlC CXXFLAGS="-O3 -qstrict -qstrict_induction -qinline -qansialias -qhot -qunroll=yes -qarch=pwr6 -qtune=pwr6 -qipa" LDFLAGS="-L$SZLIB_LIB -L$ZLIB_LIB"
	elif test "$1"=mpi; then
	#./configure --enable-64bit --with-blitz=$HOME --with-fftw3=$FFTW_INC,$FFTW_LIB --with-hdf4=$HDF4_HOME --enable-mpi --enable-optimize CC=xlc CFLAGS=-maix64 CXX=xlC MPICXX=mpCC
	#./configure --enable-64bit --with-blitz=$HOME --with-fftw3=$FFTW_INC,$FFTW_LIB --with-hdf4=$HDF4_HOME --enable-mpi --disable-cxx-flags-preset CC=xlc CXX=xlC CXXFLAGS="-O0 -g -C" MPICXX=mpCC LDFLAGS="-L$SZLIB_LIB -L$ZLIB_LIB"
	#./configure --enable-64bit --with-blitz=$HOME --with-fftw3=$FFTW_INC,$FFTW_LIB --with-hdf4=$HDF4_HOME --enable-mpi --disable-cxx-flags-preset CC=xlc CXX=xlC CXXFLAGS="-O0 -g" MPICXX="skin mpCC" LDFLAGS="-L$SZLIB_LIB -L$ZLIB_LIB"
	#./configure --enable-64bit --with-blitz=$HOME --with-fftw3=$FFTW_INC,$FFTW_LIB --with-hdf4=$HDF4_HOME --enable-mpi --disable-cxx-flags-preset CC=xlc CXX=xlC CXXFLAGS="-O3 -qstrict -qstrict_induction -qinline -qansialias -qhot -qunroll=yes -qarch=pwr6 -qtune=pwr6 -qipa" MPICXX=mpCC LDFLAGS="-L$SZLIB_LIB -L$ZLIB_LIB"
	./configure --enable-64bit --with-blitz=$HOME --with-fftw3=$FFTW_INC,$FFTW_LIB --with-hdf4=$HDF4_HOME --enable-mpi --disable-cxx-flags-preset CC=xlc CXX=xlC CXXFLAGS="-O3 -qstrict -qstrict_induction -qinline -qansialias -qhot -qunroll=yes -qarch=pwr6 -qtune=pwr6 -qipa" MPICXX="skin mpCC" LDFLAGS="-L$SZLIB_LIB -L$ZLIB_LIB"
	fi
        ;;
*titan*)
  ./configure --with-fftw=$HOME/include,$HOME/lib/gcc --with-blitz=$HOME/include,$HOME/lib/gcc --with-hdf4=$HOME/include,$HOME/lib/gcc --enable-mpi CXX=g++
  ;;

snowstorm*)
	./configure --with-blitz=$HOME --with-fftw=$HOME --with-hdf4=$HOME --enable-dependency-tracking --enable-mpi --enable-optimize CXX=icc
	;;

nana*)
	./configure --with-blitz=$HOME --with-fftw=/usr/local/numerics/fftw --with-hdf4=/usr/local/hdf/hdf --enable-dependency-tracking --enable-optimize --enable-mpi CXX=aCC
	#./configure --with-blitz=$HOME --with-fftw=/usr/local/numerics/fftw --with-hd4f=/usr/local/hdf/hdf --enable-dependency-tracking CXX=aCC
	;;


*) echo No default configuration for machine `hostname`

esac
