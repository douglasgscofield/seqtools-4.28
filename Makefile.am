#
# The main SeqTools makefile.am
# Process this file with automake to create Makefile.in.
#

# Define some shortcuts to source-code subdirectories
UTILS = $(top_srcdir)/seqtoolsUtils
PFETCH = $(top_srcdir)/libpfetch
BLIXEM = $(top_srcdir)/blixemApp
DOTTER = $(top_srcdir)/dotterApp
BELVU = $(top_srcdir)/belvuApp
DOC = $(top_srcdir)/doc
USERDOC = $(top_srcdir)/doc/User_doc
TEST = $(top_srcdir)/test
EXAMPLES = $(top_srcdir)/examples


#
# Convenience libraries
#
# These get compiled into the top-level directory to make finding
# make targent dependencies easier, even though the source-code 
# lives in subdirectories
#

if USE_LIBCURL
noinst_LIBRARIES = libSeqtoolsUtils.a libpfetch.a libBlixemApp.a libDotterApp.a libBelvuApp.a
else
noinst_LIBRARIES = libSeqtoolsUtils.a libBlixemApp.a libDotterApp.a libBelvuApp.a
endif

# Define the sources for the convenience libraries
libSeqtoolsUtils_a_SOURCES = $(UTILS)/iupac.h $(UTILS)/version.h $(UTILS)/utilities.h $(UTILS)/utilities.c $(UTILS)/blxmsp.h $(UTILS)/blxmsp.c $(UTILS)/translate.c $(UTILS)/seqtoolsWebBrowser.c $(UTILS)/blxGff3Parser.h $(UTILS)/blxGff3Parser.c $(UTILS)/blxparser.h $(UTILS)/blxparser.c
libpfetch_a_SOURCES = $(PFETCH)/libpfetch-utils.c $(PFETCH)/libpfetch.c $(PFETCH)/libpfetch.h $(PFETCH)/libpfetch-cmarshal.c $(PFETCH)/libpfetch-cmarshal.h $(PFETCH)/libcurlobject.c $(PFETCH)/libcurlobject.h $(PFETCH)/libpfetch_I.h $(PFETCH)/libpfetch_P.h $(PFETCH)/libcurlobject_I.h
libBlixemApp_a_SOURCES = $(BLIXEM)/sequencecellrenderer.c $(BLIXEM)/bigpicture.c $(BLIXEM)/bigpicturegrid.c $(BLIXEM)/detailview.c $(BLIXEM)/detailviewtree.c $(BLIXEM)/blxwindow.c $(BLIXEM)/exonview.c $(BLIXEM)/coverageview.c $(BLIXEM)/blxdotter.c $(BLIXEM)/blxFetchDb.c $(BLIXEM)/blxview.h $(BLIXEM)/blixem_.h $(BLIXEM)/detailview.h $(BLIXEM)/detailviewtree.h $(BLIXEM)/sequencecellrenderer.h $(BLIXEM)/bigpicture.h $(BLIXEM)/bigpicturegrid.h $(BLIXEM)/blxwindow.h $(BLIXEM)/exonview.h $(BLIXEM)/coverageview.h $(BLIXEM)/blxdotter.h
libDotterApp_a_SOURCES = $(DOTTER)/greyramptool.c $(DOTTER)/alignmenttool.c $(DOTTER)/dotplot.c $(DOTTER)/dotter.c $(DOTTER)/dotterKarlin.c $(DOTTER)/seqtoolsExonView.c $(DOTTER)/dotter_.h $(DOTTER)/dotter.h $(DOTTER)/seqtoolsExonView.h
libBelvuApp_a_SOURCES = $(BELVU)/belvu_.h $(BELVU)/belvu.c $(BELVU)/belvuWindow.h $(BELVU)/belvuWindow.c $(BELVU)/belvuAlignment.h $(BELVU)/belvuAlignment.c $(BELVU)/belvuTree.h $(BELVU)/belvuTree.c $(BELVU)/belvuConsPlot.h $(BELVU)/belvuConsPlot.c

# libpfetch marshal stuff
libpfetch_a_marshal = $(PFETCH)/libpfetch-cmarshal.list
libpfetch_a_DEPENDENCIES = $(nodist_libpfetch_a_SOURCES)
libpfetch_a_extra_sources = $(PFETCH)/libpfetch-cmarshal.h $(PFETCH)/libpfetch-cmarshal.c $(NULL)
nodist_libpfetch_a_SOURCES = $(libpfetch_a_extra_sources)
BUILT_SOURCES = $(libpfetch_a_extra_sources)
CLEANFILES = $(BUILT_SOURCES)
PFETCH_EXTRA_DIST = $(libpfetch_a_marshal)
GENMARSHAL_COMMAND = $(GLIB_GENMARSHAL) --prefix=libpfetch_cmarshal
.list.h:
	$(GENMARSHAL_COMMAND) --header $< >$@
.list.c:
	$(GENMARSHAL_COMMAND) --body $< >$@


#
# Documentation and tests
#

# these will be installed alongside the programs
dist_doc_DATA = $(USERDOC)/todo.html $(USERDOC)/release_notes.html \
                $(USERDOC)/belvu_quick_start.html $(USERDOC)/blixem_quick_start.html $(USERDOC)/dotter_quick_start.html \
                $(USERDOC)/belvu_usage.txt $(USERDOC)/blixem_usage.txt $(USERDOC)/dotter_usage.txt \
                $(USERDOC)/Blixem_manual/Blixem_manual.pdf $(USERDOC)/Dotter_manual/Dotter_manual.pdf $(USERDOC)/Belvu_manual/Belvu_manual.pdf

# these are included in the distribution, but are not installed
EXTRA_DIST = $(top_srcdir)/INSTALL.windows $(PFETCH_EXTRA_DIST) $(TEST) $(EXAMPLES) $(DOC)/Design_notes $(DOC)/index.html $(USERDOC)/index.html $(USERDOC)/manuals.html


#
# Executables
#
if USE_LIBCURL
bin_PROGRAMS = dotter blixem blixemh belvu
else
bin_PROGRAMS = dotter blixem belvu
endif


# blixem and blixemh have the same source files, most of which are compiled into libBlixemApp.a. However, there are a 
# few that need to be compiled with different flags/libraries, and those are specified here
blixem_SOURCES = $(top_srcdir)/blixemApp/blxmain.c $(top_srcdir)/blixemApp/blxview.c $(top_srcdir)/blixemApp/blxFetch.c
dotter_SOURCES = $(top_srcdir)/dotterApp/dotterMain.c
belvu_SOURCES = $(top_srcdir)/belvuApp/belvuMain.c

# if sqlite3 is available, compile it into blixem so that sqlite-fetch is enabled
if USE_SQLITE3
AM_CPPFLAGS = $(DEPS_CFLAGS) $(DEPS_SQLITE3) -I$(top_srcdir)
BLX_LIBS = $(DEPS_LIBS) $(DEPS_SQLITE3_LIBS)
blixem_CPPFLAGS = $(AM_CPPFLAGS) -DSQLITE3
else
AM_CPPFLAGS = $(DEPS_CFLAGS) -I$(top_srcdir)
BLX_LIBS = $(DEPS_LIBS)
blixem_CPPFLAGS = $(AM_CPPFLAGS)
endif

AM_LDFLAGS = -L$(top_builddir)

blixem_LDADD = $(top_builddir)/libSeqtoolsUtils.a $(top_builddir)/libBlixemApp.a $(BLX_LIBS) -lSeqtoolsUtils -lBlixemApp
dotter_LDADD = $(top_builddir)/libSeqtoolsUtils.a $(top_builddir)/libDotterApp.a $(DEPS_LIBS) -lSeqtoolsUtils -lDotterApp
belvu_LDADD = $(top_builddir)/libSeqtoolsUtils.a $(top_builddir)/libBelvuApp.a $(DEPS_LIBS) -lSeqtoolsUtils -lBelvuApp

# Only compile the blixemh target if we have the libcurl library
if USE_LIBCURL
blixemh_SOURCES = $(blixem_SOURCES)
blixemh_CPPFLAGS = $(blixem_CPPFLAGS) -DPFETCH_HTML
blixemh_LDADD = $(top_builddir)/libSeqtoolsUtils.a $(top_builddir)/libBlixemApp.a $(top_builddir)/libpfetch.a $(BLX_LIBS) -lcurl -lSeqtoolsUtils -lBlixemApp
endif

