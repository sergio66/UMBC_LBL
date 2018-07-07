#
# build xsec tools distribution

HDIR = read_xsec

dist:
	cd .. ;                                 \
	tar -zcf $(HDIR)/$(HDIR).tgz            \
	$(HDIR)/README                          \
	$(HDIR)/read_xsec.txt                   \
	$(HDIR)/Makefile                        \
	$(HDIR)/browse_xsec.m			\
	$(HDIR)/calc_xsec.m			\
	$(HDIR)/gid2mol.m			\
	$(HDIR)/list_bands.m			\
	$(HDIR)/read_xsec.m			\
	$(HDIR)/split_xsec.m			\
	$(HDIR)/summary_xsec.m			\
	$(HDIR)/test_extrap.m			\
	$(HDIR)/viewall_xsec.m			\
	$(HDIR)/xsec98.ok

clean:
	rm read_xsec.tgz *.fig 2> /dev/null || true

