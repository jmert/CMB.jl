JULIA ?= $(shell which julia)
usecolor = $(shell test -t 0 && echo "--color=yes" || echo "--color=no")
mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

.PHONY: test

test:
	cd $(realpath $(mkfile_path)..) && \
	$(JULIA) $(usecolor) --startup-file=no --project=test/ \
		-e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()' && \
	$(JULIA) $(usecolor) --check-bounds=yes --startup-file=no --project=test/ \
		test/runtests.jl ${DOCARGS}
