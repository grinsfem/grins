BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS    = rt

rt_SOURCES = rayleigh_taylor.C
rt_LDADD = $(top_builddir)/src/libgrins.la

AM_CPPFLAGS = -I$(top_srcdir)/src/solver/include \
              -I$(top_srcdir)/src/physics/include \
              -I$(top_srcdir)/src/boundary_conditions/include \
              -I$(top_srcdir)/src/utilities/include \
              -I$(top_srcdir)/src/visualization/include \
              -I$(top_srcdir)/src/mat_props/include \
              -I$(top_srcdir)/src/bc_handling/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
               $(LIBMESH_CPPFLAGS) \
               $(GRVY_CFLAGS)

MAINTAINERCLEANFILES = Makefile.in

#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(rt_SOURCES)
	echo 'updated source license headers' >$@
