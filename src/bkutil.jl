module BKUtil

export
    plot_healpix_map, reduc_covmat,
    BKDATA_DIR

using ..PixelCovariance
using ..Mapping
using ..Healpix: UNSEEN
using HDF5
using PyPlot
using PyCall

const BKDATA_DIR = normpath(joinpath(dirname(@__FILE__()), "..", "data"))

const healpy = PyCall.PyNULL()
const m_bicep = BicepMapDefn()

function __init__()
    copy!(healpy, pyimport("healpy"))
end

function plot_healpix_map(m::ECPMapPatchDefn, hmapin::Vector; kwds...)
    # healpy complains about the NaN values, so overwrite NaN with the HEALPix sentinel value
    # UNSEEN
    hmap = copy(hmapin)
    hmap[isnan.(hmap)] .= UNSEEN

    # Working in Julia, we may be given a PyPlot.Figure object instead of the figure number
    # which healpy is actually expecting. Automatically convert these if found.
    _conv_fig(p) = ( (k,v)=p; (v isa PyPlot.Figure) ? (k,v[:number]) : (k,v) )
    map!(_conv_fig, kwds, kwds)

    # Plot the figure
    healpy[:cartview](hmap; lonra=[m.lx,m.hx], latra=[m.ly,m.hy],
                    kwds...)
    ax = gca()
    # Fix the aspect ratio so that it looks much more BK like
    ax[:set_aspect](m.xdos / m.ydos)
end

plot_healpix_map(hmapin::Vector; kwds...) = plot_healpix_map(m_bicep, hmapin; kwds...)

doc"""
    red_spectrum(ℓ)

Returns the reddened $C_ℓ$ spectrum ``2π / ℓ(ℓ+1)``, which is actually flat in $D_ℓ$.
"""
red_spectrum(ℓ) = 2π / (ℓ*(ℓ+1))

module _reduc_covmat
    using ...PixelCovariance
    using HDF5

    const SENTINEL = typemin(Int)

    function wrap_error(f, args...)
        ret = try
            f(args...)
        catch ex
            (ex, catch_stacktrace())
        end
        ret
    end

    struct JobState
        cache::PixelCovarianceCache
        workers::Vector{Int}
        workqueue::RemoteChannel{Channel{Int64}}
        savequeue::RemoteChannel{Channel{Tuple{Int,Int}}}
        savebuffer::SharedArray{Float64,3}
    end

    function feed_pixinds(state::JobState)
        cache = state.cache
        workqueue = state.workqueue
        nprocs = length(state.workers)

        Npix = length(cache.pixels)
        for ii in 1:Npix
            put!(workqueue, ii)
            yield()
        end
        # To signal that all pixels have been consumed, we inject sentinel values into the
        # work queue which the workers will recognize as a signal to end
        for ii in 1:nprocs
            put!(workqueue, SENTINEL)
            yield()
        end
        nothing
    end

    function write_matrix(state::JobState, h5file::HDF5File, nwritten::Ref{Int})
        cache = state.cache
        savequeue = state.savequeue
        savebuffer = state.savebuffer

        Npix = length(cache.pixels)

        # Make the output matrix group exists
        if !exists(h5file, "matrix")
            h5matrix = g_create(h5file, "matrix")
        else
            h5matrix = h5file["matrix"]
        end

        # Use a symbol dict pointing to the matrices in the HDF5 file for the requested
        # fields.
        fielddict = Dict{Symbol,HDF5Dataset}()
        fields = PixelCovariance.FIELDMAP[cache.fields]
        fieldinds = map(f -> find(reshape(f .== PixelCovariance.FIELDMAP, :))[1], fields)
        for fld in fields
            if !exists(h5matrix, string(fld))
                # If the matrix entry doesn't already exist, create it with chunking
                # appropriate to writing out each matrix column independently.
                fielddict[fld] = d_create(h5matrix, string(fld), datatype(Float64),
                        dataspace(Npix,Npix), "chunk", (Npix,1))
            else
                fielddict[fld] = h5matrix[string(fld)]
            end
        end

        loopspec = zip(fieldinds, fields)
        ndone = 0
        # Enter the write loop
        while true
            pixind, off = take!(savequeue)
            if pixind == SENTINEL
                ndone += 1
                if ndone == length(state.workers)
                    break
                else
                    continue
                end
            end

            for (ind,fld) in loopspec
                d = fielddict[fld]
                d[:,pixind] = savebuffer[:,ind,off]
            end
            nwritten[] += 1
            yield()
        end

        nothing
    end

    function calc_covariance(state::JobState)
        # Get our own offset within the workers list so that we know which plane of the
        # shared buffer is ours to write to
        off = find(myid() .== state.workers)[1]
        covmat = view(state.savebuffer, :, :, off)

        cache = state.cache

        while true
            pixind = take!(state.workqueue)
            if pixind == SENTINEL
                put!(state.savequeue, (SENTINEL,off))
                break
            end
            selectpixel!(cache, pixind)
            pixelcovariance!(cache, covmat)
            put!(state.savequeue, (pixind,off))
            yield()
        end
        nothing
    end
end

"""
    reduc_covmat(obspix, covmatfile, [optional arguments...])

Generates a BICEP/Keck-style covariance matrix.

# Optional Arguments

In order:

* `nside = 512`
* `lmax = 700`
* `beamfwhm = 30.0`
* `np = Sys.CPU_CORES`

# Example

```julia
using CMB
obspixfile = joinpath(CMB.BKUtil.BKDATA_DIR, "BK15_obsmat_healpix_pixels.dat");
obspix = reshape(readdlm(obspixfile, '\n', Int), :);
covmatfile = joinpath(expanduser("~"), "BK15_covmat.h5");
CMB.BKUtil.reduc_covmat(obspix[1:10], covmatfile)
```
"""
function reduc_covmat(obspix, covmatfile,
                      nside=512,
                      lmax=700,
                      beamfwhm=30.0,
                      np=Sys.CPU_CORES)

    # Open the output HDF5 file
    h5file = h5open(covmatfile, "w");

    println("Initializing pixel-pixel covariance cache...")
    cache = PixelCovarianceCache(nside, lmax, obspix, [:QQ,:QU,:UQ,:UU])
    makespectra!(cache, red_spectrum)
    applybeam!(cache, beamfwhm)

    # Communication channels
    workqueue = RemoteChannel(() -> Channel{Int}(np))
    savequeue = RemoteChannel(() -> Channel{Tuple{Int,Int}}(np))

    # Startup parallel workers
    println("Initializing $np worker processes...")
    navail = nprocs() > 1 ? nworkers() : 0
    if navail < np
        addprocs(np - navail)
    end
    wrkids = workers()
    @sync for pp in wrkids
        @async remotecall_fetch(Core.eval, pp, Main, :(using CMB, CMB.BKUtil))
    end

    # Allocate the shared buffer
    savebuffer = SharedArray{Float64}(length(obspix), 9, np; pids=wrkids)

    # Group up everything into the state object which all tasks will use.
    state = _reduc_covmat.JobState(cache, wrkids, workqueue, savequeue, savebuffer)
    nwritten = Ref{Int}(0)

    println("Launching jobs...")

    # Start the feed and save tasks
    feedtask = @schedule _reduc_covmat.wrap_error(_reduc_covmat.feed_pixinds, state)
    savetask = @schedule _reduc_covmat.wrap_error(_reduc_covmat.write_matrix, state, h5file, nwritten)

    worktasks = Vector{Task}()
    for pp in wrkids
        push!(worktasks, @async remotecall_fetch(_reduc_covmat.wrap_error, pp,
                                                 _reduc_covmat.calc_covariance, state))
    end

    alltasks = append!([feedtask,savetask], worktasks)
    try
        while !istaskdone(feedtask)
            for task in alltasks
                if istaskdone(task)
                    ret = wait(task)
                    ret isa Tuple{Exception,Array{StackFrame}} && rethrow(ret)
                end
            end
            sleep(1)
            println("$(nwritten[]) of $(length(obspix)) columns computed.")
        end
    catch ex
        interrupt(wrkids)
        if ex isa Tuple{Exception,Array{StackFrame}}
            showerror(STDERR, ex[1], ex[2])
        else
            showerror(STDERR, ex, catch_backtrace())
        end
    finally
        close(h5file)
        close(savequeue)
        close(workqueue)
    end
end

end
