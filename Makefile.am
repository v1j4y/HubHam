AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 --install

VERSION_MAJOR   = @VERSION_MAJOR@
VERSION_MINOR   = @VERSION_MINOR@
VERSION_PATCH   = @VERSION_PATCH@

src_main_c = src/ex1.c src/hubbard.c src/readgraphmllib.c src/utils.c src/get_s2.c
src_main_h = src/hubbard.h src/readgraphmllib.h src/utils.h src/get_s2.h

bin_PROGRAMS = $(top_builddir)/bin/ex1
__top_builddir__bin_ex1_SOURCES = $(src_main_c)
__top_builddir__bin_ex1_CPPFLAGS = $(AM_CPPFLAGS)
__top_builddir__bin_ex1_CFLAGS = $(AM_CFLAGS)
__top_builddir__bin_ex1_CFLAGS += "-I${top_srcdir}/include"  "-I${top_srcdir}/src"
__top_builddir__bin_ex1_LDFLAGS = $(AM_LDFLAGS)
