JULIA ?= $(shell which julia)
usecolor = $(shell test -t 0 && echo "--color=yes" || echo "--color=no")

.DUMMY: test

test:
	$(JULIA) $(usecolor) runtests.jl
