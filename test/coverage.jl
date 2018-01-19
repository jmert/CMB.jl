# Only run non-nightlies on Linux
get(ENV, "TRAVIS_OS_NAME", "")       == "linux"   || exit()
get(ENV, "TRAVIS_JULIA_VERSION", "") == "nightly" && exit()

Pkg.add("Coverage")
using Coverage

cd(joinpath(dirname(@__FILE__), "..")) do
    #Coveralls.submit(Coveralls.process_folder())
    Codecov.submit(Codecov.process_folder())
end

