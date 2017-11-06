JULIA ?= $(shell which julia)
usecolor = $(shell test -t 0 && echo "--color=yes" || echo "--color=no")

.DUMMY: test

test:
	$(JULIA) $(usecolor) --check-bounds=yes \
		-e 'Pkg.build("CMB"); Pkg.test("CMB");'
