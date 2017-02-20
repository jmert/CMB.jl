module CMBTests
    using Base.Test

    const MODULES = Dict(
        "harmonics"       => :Harmonics,
        "healpix"         => :Healpix,
        "pixelcovariance" => :PixelCovariance)
    const TESTLIST = Dict{typeof(""), Function}()

    export load!, loadall!, runtests

    ###################
    # HELPER FUNCTIONS
    ###################

    function load!(testlist::Dict, id::AbstractString)
        modsym = MODULES[id]
        modpath = joinpath(dirname(@__FILE__), "$(id).jl");
        # Include the file and have it processed within this module
        eval(CMBTests, :( include($modpath) ))
        # Add the list of tests to run to the list
        testlist[id] = eval(CMBTests, modsym).runtests
    end
    load!(id::AbstractString) = load!(TESTLIST, id)

    function loadall!(testlist::Dict; verbose::Bool=true)
        for id in keys(MODULES)
            verbose && print("loading group $(repr(id))..."); tic();
            load!(testlist, id)
            verbose && println("done (took $(toq()) seconds)")
        end
        testlist
    end
    loadall!(; kwargs...) = loadall!(TESTLIST; kwargs...)

    function runtests()
        @testset "$id" for id in keys(TESTLIST)
            TESTLIST[id]()
        end
    end
end

