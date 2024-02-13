
ACLOCAL_AMFLAGS  = -I m4

SUBDIRS = mudfas src lib tests matlab doc share

EXTRA_DIST = ChangeLog.0 m4/ChangeLog \
config/sites/f14_ia32_ifort+mpich2_optim_jazz2 \
config/sites/f14_ia32_ifort+mpich2+skin_optim_jazz2 \
config/sites/f14_ia32_ifort_optim_jazz2 \
config/sites/sp6_xlC_cineca \
config/sites/sp6_xlC+mpCC_cineca \
config/sites/sp6_xlC+mpCC+skin_cineca

## Clean up template repositories
clean-local:
	-rm -rf repository ti_files ii_files


update-cl:
	cvs2cl.pl -P --no-times --no-wrap --no-common-dir -U $(HOME)/CVSROOT/users --accum \
	-I m4 -I mudfas -I matlab/@sd -l "-d>2006-08-28"

