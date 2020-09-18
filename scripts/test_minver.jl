using Logging, Pkg
import Pkg.Types: VersionRange, semver_spec
try
    using TOML
    import TOML: parsefile
catch
    import Pkg.TOML: parsefile
end

function minimize_versions()
    # Read project file and extract the compat entries
    proj = parsefile(joinpath(@__DIR__, "..", "Project.toml"))
    compat = proj["compat"]

    # Transform the semver strings to Pkg's VersionSpec, and then from those extract the
    # minimum bound and set each key to the minimum string.
    map!(semver_spec, values(compat))
    for (proj, vers) in compat
        isempty(vers.ranges) && continue
        v = VersionRange(vers.ranges[1].lower)
        compat[proj] = sprint(print, v)
    end

    # Only continue running this script on the minimum-supported Julia version
    minjulia = VersionNumber(get(compat, "julia", "1"))
    minjulia = VersionNumber("$(minjulia.major).$(minjulia.minor+1)")
    VERSION ≤ minjulia || return

    # Now pin every package version to the minimum bound (excluding Julia itself)
    delete!(compat, "julia")
    for (proj, vers) in compat
        @info "Pinning $proj at v$vers..."
        Pkg.pin(PackageSpec(name = proj, version = vers))
    end
end

minimize_versions()
