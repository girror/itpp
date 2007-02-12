include $(top_srcdir)/doc/local/sources.mk
include $(top_srcdir)/doc/tutorial/sources.mk
include $(top_srcdir)/doc/tutorial/src/sources.mk
include $(top_srcdir)/itpp/sources.mk
include $(top_srcdir)/itpp/base/sources.mk
include $(top_srcdir)/itpp/base/bessel/sources.mk
include $(top_srcdir)/itpp/comm/sources.mk
include $(top_srcdir)/itpp/fixedpoint/sources.mk
include $(top_srcdir)/itpp/protocol/sources.mk
include $(top_srcdir)/itpp/srccode/sources.mk

SUBDIRS = images tutorial local

docdir = @docdir@

distdir_html = @PACKAGE@-html-doc-@VERSION@

all-local: html-local

html-local: html/index.html

html/index.html: doxygen_html.cfg \
	$(doc_local_sources) $(html_local_sources) \
	$(doc_tutorial_sources) $(cpp_tutorial_sources) \
	$(h_sources) $(cpp_sources) \
	$(h_base_sources) $(cpp_base_sources) \
	$(h_base_bessel_sources) $(cpp_base_bessel_sources) \
	$(h_comm_sources) $(cpp_comm_sources) \
	$(h_fixed_sources) $(cpp_fixed_sources) \
	$(h_protocol_sources) $(cpp_protocol_sources) \
	$(h_srccode_sources) $(cpp_srccode_sources)
	doxygen $<; \
	if test -d html; then \
		cp $(srcdir)/images/itpp_logo.png html; \
		cp $(srcdir)/images/newcom_logo.png html; \
		cp $(srcdir)/images/favicon.ico html; \
	fi

dist-html: html-local
	cp -a $(top_builddir)/doc/html $(top_builddir)/$(distdir_html)
	tar chof - $(top_builddir)/$(distdir_html) | \
		gzip --best -c > $(top_builddir)/$(distdir_html).tar.gz
	tar chof - $(top_builddir)/$(distdir_html) | \
		bzip2 -9 -c > $(top_builddir)/$(distdir_html).tar.bz2
	rm -rf $(top_builddir)/$(distdir_html)

install-data-local: html-local
	$(mkinstalldirs) $(DESTDIR)$(docdir)
	if test -d html; then \
		$(mkinstalldirs) $(DESTDIR)$(docdir)/html; \
		for f in html/*; do \
			$(INSTALL_DATA) $$f $(DESTDIR)$(docdir)/html; \
		done \
	fi

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)

clean-local:
	rm -rf html *.log *.tag
