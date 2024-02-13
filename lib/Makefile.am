
picsim_srcdir = $(top_srcdir)/src
picsim_libdir = $(top_srcdir)/lib
mudfas_srcdir = $(top_srcdir)/mudfas/src
cpptools_srcdir = $(top_srcdir)/mudfas/src/c++-tools
mudfas_libdir = $(top_srcdir)/mudfas/lib

AM_CPPFLAGS = -I$(picsim_srcdir) -I$(mudfas_srcdir) -I$(cpptools_srcdir)

picsim_src = \
$(picsim_srcdir)/backgrd.cpp \
$(picsim_srcdir)/beam.cpp  \
$(picsim_srcdir)/bbeam.cpp  \
$(picsim_srcdir)/diagnostics.cpp \
$(picsim_srcdir)/drivenbackgrd.cpp\
$(picsim_srcdir)/drivennonhombackgrd.cpp\
$(picsim_srcdir)/generic-picsim.cpp \
$(picsim_srcdir)/moments.cpp \
$(picsim_srcdir)/moments-handler.cpp \
$(picsim_srcdir)/nonhombackgrd.cpp \
$(picsim_srcdir)/phase-space.cpp \
$(picsim_srcdir)/phase-space-handler.cpp \
$(picsim_srcdir)/picsime.cpp \
$(picsim_srcdir)/picsimei.cpp \
$(picsim_srcdir)/picsimi.cpp \
$(picsim_srcdir)/probes.cpp \
$(picsim_srcdir)/probes-handler.cpp \
$(picsim_srcdir)/scheduler.cpp \
$(picsim_srcdir)/species.cpp \
$(picsim_srcdir)/species-handler.cpp \
$(picsim_srcdir)/spectra.cpp \
$(picsim_srcdir)/spectra-handler.cpp \
$(picsim_srcdir)/surface.cpp \
$(picsim_srcdir)/surface-handler.cpp

pdf_src = \
$(picsim_srcdir)/pdf.cpp

lib_LTLIBRARIES = libpicsim2d.la libpicsim3d.la libpdf.la

if COMPILER_IS_COMPAQ_CXX
templates_picsim2d_dir = -ptr $(picsim_libdir)/picsim2d
templates_picsim3d_dir = -ptr $(picsim_libdir)/picsim3d
templates_pdf_dir      = -ptr $(picsim_libdir)/pdf
endif

## 2-dimensions library
libpicsim2d_la_SOURCES = $(picsim_src)
libpicsim2d_la_CXXFLAGS = -DDIMR=2 -DDIMV=3 $(templates_picsim2d_dir)
libpicsim2d_la_LDFLAGS = -version-info 0:0:0

## 3-dimensions library
libpicsim3d_la_SOURCES = $(picsim_src)
libpicsim3d_la_CXXFLAGS = -DDIMR=3 -DDIMV=3 $(templates_picsim3d_dir)
libpicsim3d_la_LDFLAGS = -version-info 0:0:0

## random generator library
libpdf_la_SOURCES = $(pdf_src)
libpdf_la_CXXFLAGS = $(templates_pdf_dir)
libpdf_la_LDFLAGS = -version-info 0:0:0

## Clean up template repositories
clean-local:
	-rm -rf picsim2d picsim3d pdf ti_files ii_files

