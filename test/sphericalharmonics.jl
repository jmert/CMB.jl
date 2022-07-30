using CMB.SphericalHarmonics: SHT, HealpixPixelization
include(joinpath(dirname(pathof(SHT)), "..", "test", "testsuite.jl"))

TestSuite.runtests(HealpixPixelization(16), rtol = 1e-6)
