
EXTRA_DIST = Picsim.cpp.in pdftest.m kh2dgrowth.py.in

picsim_srcdir = $(top_srcdir)/src
picsim_libdir = $(top_srcdir)/lib
mudfas_srcdir = $(top_srcdir)/mudfas/src
cpptools_srcdir = $(top_srcdir)/mudfas/src/c++-tools
mudfas_libdir = $(top_srcdir)/mudfas/lib
AM_CPPFLAGS = -I$(picsim_srcdir) -I$(mudfas_srcdir) -I$(cpptools_srcdir)

CLEANFILES = Picsim.cpp Picsimi.cpp Picsime.cpp

# create sources
Picsimi.cpp: Picsim.cpp.in
	$(SED) -e 's/@class_header@/picsimi.h/' -e 's/@class_name@/PicSimI/' $< > $@

Picsime.cpp: Picsim.cpp.in
	$(SED) -e 's/@class_header@/picsime.h/' -e 's/@class_name@/PicSimE/' $< > $@

Picsimei.cpp: Picsim.cpp.in
	$(SED) -e 's/@class_header@/picsimei.h/' -e 's/@class_name@/PicSimEI/' $< > $@

# Executable targets
check_PROGRAMS = List Slist Vector Deque Map Blitz Pdf Fftw FactorySpecies \
testMT1 testMT2 testRNG
bin_PROGRAMS = Picsime2d Picsime3d Picsimei2d Picsimei3d Picsimi2d Picsimi3d 

bin_SCRIPTS = bench_container kh2dgrowth

if COMPILER_IS_COMPAQ_CXX
templates_picsim2d_dir = -ptr $(picsim_libdir)/picsim2d \
-ptr $(mudfas_libdir)/mudfas2d
templates_picsim3d_dir = -ptr $(picsim_libdir)/picsim3d \
-ptr $(mudfas_libdir)/mudfas3d
templates_pdf_dir = -ptr $(picsim_libdir)/pdf
templates_mudfas_utils_dir = -ptr $(mudfas_libdir)/mudfas-utils
templates_local_dir = -ptr cxx_repository
endif

List_SOURCES = List.cpp
Slist_SOURCES = Slist.cpp
Vector_SOURCES = Vector.cpp
Deque_SOURCES = Deque.cpp
Map_SOURCES = Map.cpp
Blitz_SOURCES = Blitz.cpp
Fftw_SOURCES = Fftw.cpp

Pdf_SOURCES = Pdf.cpp 
Pdf_LDFLAGS = $(templates_local_dir) $(templates_pdf_dir) \
$(templates_mudfas_utils_dir)
Pdf_LDADD = $(picsim_libdir)/libpdf.la $(mudfas_libdir)/libmudfas-utils.la

FactorySpecies_SOURCES = FactorySpecies.cpp
FactorySpecies_CPPFLAGS = -DDIMR=3 -DDIMV=3 $(AM_CPPFLAGS)
FactorySpecies_LDFLAGS = $(templates_local_dir) $(templates_picsim2d_dir) 
FactorySpecies_LDADD = $(picsim_libdir)/libpicsim2d.la \
$(picsim_libdir)/libpdf.la $(mudfas_libdir)/libmudfas-utils.la

testMT1_SOURCES = testMT1.cpp
testMT2_SOURCES = testMT2.cpp
testRNG_SOURCES = testRNG.cpp

# Picsime
Picsime2d_SOURCES = Picsime.cpp
Picsime2d_CPPFLAGS = -DDIMR=2 -DDIMV=3 $(AM_CPPFLAGS)
Picsime2d_LDFLAGS = $(templates_local_dir) $(templates_picsim2d_dir) \
$(templates_pdf_dir) $(templates_mudfas_utils_dir)
Picsime2d_LDADD = $(picsim_libdir)/libpicsim2d.la $(picsim_libdir)/libpdf.la \
$(mudfas_libdir)/libmudfas2d.la $(mudfas_libdir)/libmudfas-utils.la

Picsime3d_SOURCES = Picsime.cpp
Picsime3d_CPPFLAGS = -DDIMR=3 -DDIMV=3 $(AM_CPPFLAGS)
Picsime3d_LDFLAGS = $(templates_local_dir) $(templates_picsim3d_dir) \
$(templates_pdf_dir) $(templates_mudfas_utils_dir)
Picsime3d_LDADD = $(picsim_libdir)/libpicsim3d.la $(picsim_libdir)/libpdf.la \
$(mudfas_libdir)/libmudfas3d.la $(mudfas_libdir)/libmudfas-utils.la

# Picsimei
Picsimei2d_SOURCES = Picsimei.cpp
Picsimei2d_CPPFLAGS = -DDIMR=2 -DDIMV=3 $(AM_CPPFLAGS)
Picsimei2d_LDFLAGS = $(templates_local_dir) $(templates_picsim2d_dir) \
$(templates_pdf_dir) $(templates_mudfas_utils_dir)
Picsimei2d_LDADD = $(picsim_libdir)/libpicsim2d.la $(picsim_libdir)/libpdf.la \
$(mudfas_libdir)/libmudfas2d.la $(mudfas_libdir)/libmudfas-utils.la

Picsimei3d_SOURCES = Picsimei.cpp
Picsimei3d_CPPFLAGS = -DDIMR=3 -DDIMV=3 $(AM_CPPFLAGS)
Picsimei3d_LDFLAGS = $(templates_local_dir) $(templates_picsim3d_dir) \
$(templates_pdf_dir) $(templates_mudfas_utils_dir)
Picsimei3d_LDADD = $(picsim_libdir)/libpicsim3d.la $(picsim_libdir)/libpdf.la \
$(mudfas_libdir)/libmudfas3d.la $(mudfas_libdir)/libmudfas-utils.la

# Picsimi
Picsimi2d_SOURCES = Picsimi.cpp
Picsimi2d_CPPFLAGS = -DDIMR=2 -DDIMV=3 $(AM_CPPFLAGS)
Picsimi2d_LDFLAGS = $(templates_local_dir) $(templates_picsim2d_dir) \
$(templates_pdf_dir) $(templates_mudfas_utils_dir)
Picsimi2d_LDADD = $(picsim_libdir)/libpicsim2d.la $(picsim_libdir)/libpdf.la \
$(mudfas_libdir)/libmudfas2d.la $(mudfas_libdir)/libmudfas-utils.la

Picsimi3d_SOURCES = Picsimi.cpp
Picsimi3d_CPPFLAGS = -DDIMR=3 -DDIMV=3 $(AM_CPPFLAGS)
Picsimi3d_LDFLAGS = $(templates_local_dir) $(templates_picsim3d_dir) \
$(templates_pdf_dir) $(templates_mudfas_utils_dir)
Picsimi3d_LDADD = $(picsim_libdir)/libpicsim3d.la $(picsim_libdir)/libpdf.la \
$(mudfas_libdir)/libmudfas3d.la $(mudfas_libdir)/libmudfas-utils.la

## Clean up template repositories
clean-local:
	-rm -rf cxx_repository ti_files ii_files

## Check 
check:
	@echo Running test suite...
	@for prog in $(check_PROGRAMS); do \
		echo $$prog && (./$$prog > /dev/null 2>&1) || exit 1; \
		done
	@echo All tests passed.

