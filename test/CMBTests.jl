module CMBTests
    using Base.Test

    const MODULES = Dict(
        "legendre" => "Legendre",
        "harmonics" => "Harmonics",
        "healpix" => "HEALPix",
        "sphere" => "Sphere"
       )
    TESTLIST = Dict{typeof(""), typeof("")}()

    export load!, loadall!, runtests

    function load!(testlist, id)
        testlist[id] = MODULES[id]
    end
    load!(id::AbstractString) = load!(TESTLIST, id)

    function loadall!(testlist)
        for tt in MODULES
            push!(testlist, tt)
        end
        testlist
    end
    loadall!() = loadall!(TESTLIST)

    function runtests()
        @testset "CMB" begin
            @testset "$desc" for (id,desc) in TESTLIST
                modpath = joinpath(dirname(@__FILE__), "$(id).jl")
                # Include the file and have it processed within this module
                print("running $desc tests... ")
                tic()
                eval(CMBTests, :(include($modpath)) )
                toc()
            end
        end
    end
end

