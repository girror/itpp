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

EXTRA_DIST = autogen.sh ChangeLog-2005 ChangeLog-2006 itpp.spec \
	Makefile.common NEWS-3.9 TODO

man_MANS = @PACKAGE@-config.1
dist_man_MANS = @PACKAGE@-config.1

bin_SCRIPTS = @PACKAGE@-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @PACKAGE@.pc

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/@PACKAGE@

dist-html:
	$(MAKE) -C doc dist-html

snapshot:
	$(MAKE) dist distdir=@PACKAGE_TARNAME@-`date +"%Y%m%d"`

snapshot-html:
	$(MAKE) -C doc dist-html html_tarname=@PACKAGE_TARNAME@-html-doc-`date +"%Y%m%d"`
