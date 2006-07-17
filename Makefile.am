include $(top_srcdir)/Makefile.common

ACLOCAL_AMFLAGS = -I config

AUTOMAKE_OPTIONS = dist-bzip2

if HTML_DOC
HTMLDOCDIR = doc
endif

if RUN_TESTS
TESTDIR = tests
endif

SUBDIRS = itpp extras win32 $(HTMLDOCDIR) $(TESTDIR)

EXTRA_DIST = autogen.sh itpp.spec Makefile.common TODO

bin_SCRIPTS = @PACKAGE@-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @PACKAGE@.pc

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/@PACKAGE@

disthtml:
	$(MAKE) -C doc disthtml
