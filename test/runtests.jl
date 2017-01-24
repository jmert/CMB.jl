# Make sure the module source will be in the path whether this is run from
# the julia REPL after installation or from a Makefile without being
# installed.
push!(LOAD_PATH, joinpath(dirname(@__FILE__), "..", "src"));
# Then add the testing path so we can initialize the testing module
push!(LOAD_PATH, joinpath(dirname(@__FILE__)));

import CMBTests
CMBTests.loadall!();
CMBTests.runtests();

