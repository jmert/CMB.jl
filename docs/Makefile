.PHONY: build test

DOCARGS :=
JULIA ?= julia
PORT ?= 8000

mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
package_dir := $(realpath $(mkfile_path)..)

export HASTEST

default: preview

preview: prepare
	cd $(package_dir) && \
		$(JULIA) --startup-file=no --optimize=0 --project=docs/ \
		-e "using LiveServer; servedocs(port=$(PORT))"

build: prepare
	cd $(package_dir) && \
	$(JULIA) --startup-file=no --optimize=0 --project=docs/ \
		docs/make.jl ${DOCARGS}

prepare:
	cd $(package_dir) && \
	$(JULIA) --startup-file=no --optimize=0 --project=docs/ \
		-e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

test: DOCARGS+= --test
test: build
	@true

fix: DOCARGS+= --fix
fix: build
	@true
