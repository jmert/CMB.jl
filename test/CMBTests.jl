module CMBTests
    using Compat.Test

    const MODULES = Dict(
        "sphere" => "Sphere",
        "legendre" => "Legendre",
        "healpix" => "HEALPix",
        "harmonics" => "Harmonics",
        "sparse" => "Sparse"
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
                t0 = time_ns()
                eval(CMBTests, :(include($modpath)) )
                t1 = time_ns()
                println( (t1-t0)/1e9, " seconds")
            end
        end
    end
end

