"""
dev_pkg_tools.jl

Utility helpers for Julia package development:
1. Create a new package
2. Add dependencies
3. Resolve dependency/version issues
4. Update dependencies cleanly
5. Work in dev mode with Revise
"""

module DevPkgTools

export create_package,
       add_deps,
       resolve_deps,
       update_deps,
       dev_with_revise
using Pkg

# ------------------------------------------------------------
# 1. Create a new package
# ------------------------------------------------------------

"""
    create_package(name; path=pwd())

Create a new Julia package with the given name.
"""
function create_package(name::AbstractString; path=pwd())
    Pkg.activate(path)
    Pkg.generate(name)
    @info "Package `$name` created at $(joinpath(path, name))"
end


# ------------------------------------------------------------
# 2. Add dependencies to an existing package
# ------------------------------------------------------------

"""
    add_deps(pkgs)

Add one or more dependencies to the active package environment.
"""
function add_deps(pkgs)
    pkgs = isa(pkgs, AbstractVector) ? pkgs : [pkgs]
    Pkg.add(pkgs)
    @info "Added dependencies: $(join(pkgs, ", "))"
end


# ------------------------------------------------------------
# 3. Resolve dependency/version conflicts
# ------------------------------------------------------------

"""
    resolve_deps(; force=false)

Re-resolve the dependency graph.
Set `force=true` to rebuild everything from Project.toml.
"""
function resolve_deps(; force::Bool=false)
    if force
        @info "Force-resolving dependencies from Project.toml"
        Pkg.resolve(; force=true)
    else
        @info "Resolving dependencies"
        Pkg.resolve()
    end
end


# ------------------------------------------------------------
# 4. Update all dependencies and regenerate Manifest.toml
# ------------------------------------------------------------

"""
    update_deps(; clean_manifest=true)

Update all dependencies to latest compatible versions.
Optionally delete Manifest.toml first.
"""
function update_deps(; clean_manifest::Bool=true)
    if clean_manifest && isfile("Manifest.toml")
        rm("Manifest.toml")
        @info "Removed Manifest.toml"
    end

    Pkg.update()
    Pkg.resolve()
    Pkg.instantiate()
    Pkg.precompile()

    @info "Dependencies updated and manifest regenerated"
end


# ------------------------------------------------------------
# 5. Use package in dev mode with Revise
# ------------------------------------------------------------

"""
    dev_with_revise()

Activate dev environment, load Revise, and load the package.
"""
function dev_with_revise()
    # using Pkg

    Pkg.activate(".")

    pkgname = Pkg.project().name
    pkgname === nothing && error("Could not determine package name")

    @info "Developing $pkgname with Revise"

    Base.eval(Main, :(using Revise))
    Base.eval(Main, :(using $(Symbol(pkgname))))
end


end # module
