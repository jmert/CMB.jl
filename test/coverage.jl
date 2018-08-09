# Only run non-nightlies on Linux
Sys.ARCH === :x86_64 && Sys.KERNEL === :Linux || exit()
get(ENV, "TRAVIS_JULIA_VERSION", "") == "nightly" && exit()

try
    import Pkg
    Pkg.add("Coverage")
    using Coverage

    cd(joinpath(dirname(@__FILE__), "..")) do
        #Coveralls.submit(Coveralls.process_folder())
        Codecov.submit(Codecov.process_folder())
    end
finally
    Pkg.rm("Coverage")
end
