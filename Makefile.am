AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 --install

VERSION_MAJOR   = @VERSION_MAJOR@
VERSION_MINOR   = @VERSION_MINOR@
VERSION_PATCH   = @VERSION_PATCH@

src_main_c = src/ex1.c src/hubbard.c src/readgraphmllib.c
src_main_h = src/hubbard.h src/readgraphmllib.h

bin_PROGRAMS = ex1
ex1_SOURCES = $(src_main_c)
ex1_CPPFLAGS = $(AM_CPPFLAGS)
ex1_CFLAGS = $(AM_CFLAGS)
ex1_CFLAGS += "-I${top_srcdir}/include"  "-I${top_srcdir}/src"
ex1_LDFLAGS = $(AM_LDFLAGS)
