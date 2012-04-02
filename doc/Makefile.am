MAINTAINERCLEANFILES = Makefile.in html/* latex/*

EXTRA_DIST = html/index.html html latex/refman.tex latex

html/index.html latex/refman.tex: $(top_srcdir)/src/*.cxx $(top_srcdir)/src/*.h
	doxygen
