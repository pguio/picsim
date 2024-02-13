dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])

AC_INIT([picsim],[2.1.3],[patrick.guio@gmail.com])

AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/species.cpp])
dnl AC_CONFIG_HEADERS([src/picsim-config.h])

AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([1.9 dist-bzip2 subdir-objects])

dnl Avoid generating makefile rules to rebuild maintainer-only files by
dnl default.  Maintainers may override this default and generate these
dnl makefile rules using the `--enable-maintainer-mode' configure option.
AM_MAINTAINER_MODE

AC_MSG_NOTICE
AC_MSG_NOTICE([Configuring $PACKAGE_STRING for $target])
AC_MSG_NOTICE

dnl check for sed required to create Picsim*.cpp
AC_CHECK_PROGS(SED, sed, )

dnl C++ compiler
AC_LANG([C++])
AC_PROG_CXX([CC cxx xlC aCC g++ icc pgCC KCC FCC])
AM_CONDITIONAL(COMPILER_IS_COMPAQ_CXX, [test x"$CXX" = xcxx])
AC_CXX_FLAGS_PRESET
AC_COMPILERS_ENABLE_64BITS
AC_CXX_ENABLE_OPTIMIZE
AC_CXX_ENABLE_DEBUG
AC_CXX_ENABLE_PROFILE

dnl Doxygen documentation
AC_PROG_DOXYGEN

dnl Disable shared library building
AC_DISABLE_SHARED([])
LT_INIT

dnl Checks for libraries

dnl Blitz path and derived variables
AC_CHECK_BLITZ([
CPPFLAGS="$BLITZ_CPPFLAGS $CPPFLAGS"
LDFLAGS="$BLITZ_LDFLAGS $LDFLAGS"
LIBS="$BLITZ_LIBS $LIBS"
],[AC_MSG_ERROR_BLITZ])

dnl HDF library
AC_MSG_CHECKING([whether to enable hdf5 output])
  AC_ARG_ENABLE(hdf5-output, AS_HELP_STRING([--enable-hdf5-output],
               [Enable hdf5 output]),
               [enable_hdf5_output="yes"],
               [enable_hdf5_output="no"])
AC_MSG_RESULT([$enable_hdf5_output])
if test "$enable_hdf5_output" = no ; then
  AC_CHECK_HDF4([
  CPPFLAGS="$HDF4_CPPFLAGS $CPPFLAGS"
  LDFLAGS="$HDF4_LDFLAGS $LDFLAGS"
  LIBS="$HDF4_LIBS $LIBS"
  ],[AC_MSG_ERROR_HDF4])
else
  AC_CHECK_HDF5([
  CPPFLAGS="$HDF5_CPPFLAGS $CPPFLAGS"
  LDFLAGS="$HDF5_LDFLAGS $LDFLAGS"
  LIBS="$HDF5_LIBS $LIBS"
  ],[AC_MSG_ERROR_HDF5])
fi

dnl FFT library
AC_CHECK_FFTW
AC_CHECK_FFTW3
AC_CHECK_MLIB
AC_CHECK_DXML

CPPFLAGS="$FFTW_CPPFLAGS $FFTW3_CPPFLAGS $MLIB_CPPFLAGS $DXML_CPPFLAGS $CPPFLAGS"
LDFLAGS="$FFTW_LDFLAGS $FFTW3_LDFLAGS $MLIB_LFFLAGS $DXML_LDFLAGS $LDFLAGS"
LIBS="$FFTW_LIBS $FFTW3_LIBS $MLIB_LIBS $DXML_LIBS $LIBS"


AC_ARG_ENABLE(double-coordinate,
AS_HELP_STRING([--enable-double-coordinate],[enable the double type for coordinate]),
[case "$enableval" in
yes) AC_MSG_NOTICE([Double enabled for coordinate])
	AC_DEFINE(DOUBLE_COORDINATE,1,[Double enabled for coordinate type])
	;;
no)
	AC_DEFINE(FLOAT_COORDINATE,1,[Float enabled for coordinate type])
	;;
esac],[AC_DEFINE(FLOAT_COORDINATE,1,[Float enabled for coordinate type])])

AC_ARG_ENABLE(double-field,
AS_HELP_STRING([--enable-double-field],[enable the double type for field]),
[case "$enableval" in
yes) AC_MSG_NOTICE([Double enabled for field])
	AC_DEFINE(DOUBLE_FIELD, 1, [Double enabled for field type])
	;;
no)
	AC_DEFINE(FLOAT_FIELD, 1, [Float enabled for field type])
	;;
esac],[AC_DEFINE(FLOAT_FIELD, 1, [Float enabled for field type])])

dnl MPI support
AC_CXX_MPI


dnl Checks for header files.
AC_CXX_HAVE_STD
AC_CHECK_HEADERS([sys/times.h sys/resource.h])
dnl If slist not present do not try to compile Slist.cpp
AC_CHECK_HEADER(slist, AC_DEFINE(HAVE_SLIST, 1, [C++ slist support]))

dnl Config mudfas
AC_CONFIG_SUBDIRS(mudfas)

dnl Print variables
AC_ENV

dnl Generate Makefiles
AC_CONFIG_FILES([
		Makefile
		src/Makefile 
		lib/Makefile 
		tests/Makefile 
		matlab/Makefile
		doc/Makefile 
		doc/doxygen/Doxyfile
		share/Makefile
		share/templates/Makefile]) 
	
AC_CONFIG_FILES([tests/bench_container:tests/bench_container.sh.in], 
		[chmod +x tests/bench_container])
AC_CONFIG_FILES([tests/kh2dgrowth:tests/kh2dgrowth.sh.in], [chmod +x tests/kh2dgrowth])
AC_CONFIG_FILES([makelinks:makelinks.sh.in], [chmod +x makelinks])
AC_CONFIG_FILES([doc/version.sty:doc/version.sty.in])


AC_OUTPUT
