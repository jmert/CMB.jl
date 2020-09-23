# [Observing Matrices](@id man_obsmat)

```@contents
Pages = ["obsmat.md"]
Depth = 2
```

```@raw html
<div style="display:none;">
```
```math
\newcommand\mat[1]{\mathbf{#1}}
\newcommand\expv[1]{\left\langle #1\right\rangle}
```
```@raw html
</div>
```

## Introduction

A fundamental fact about CMB observations is that any given map is filtered in some way.
For instance, ground-based experiments apply timestream filters to mitigate the
contribution of atmospheric noise, and even full-sky experiments such as WMAP and Planck
must mask out the galactic plane.

Generically, any map making procedure can be encoded as some function which transforms
an idealized (pixelized) map ``\mat{\hat m}`` into an observed map ``\mat{\tilde m}``.
If the function is a linear operation, then it can be encoded as a so-called
_observing matrix_ ``\mat R`` with
```math
    \mat{\tilde m} = \mat R \mat{\hat m}
```
This section describes the tools provided by `CMB.jl` for working with obsering matrices.
These tools are specifically motivated by the use of observing matrices in generating
[_reobserved pixel-pixel covariance matrices_](@ref pixelcov_reobs).

A detailed description of composing an observing matrix is given in
[BICEP/Keck Array VII](@ref bib-bicepkeck) (as well as the theses of J. Tolan and
J. Willmert).

## Loading an Observing Matrix

It is assumed that all observing matrices are sparse, and the native format used by this
package is of
[compressed sparse column (CSC)](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS))
matrices in accordance with Julia's `SparseMatrixCSC` type from the `SparseArrays`
standard library.
Observing matrices are read from disk with the [`read_obsmat`](@ref) function, and the
following formats are supported:

1. An HDF5 file written by [`write_obsmat`](@ref). (See the next section.)

2. `numpy` sparse matrices saved to an HDF5 file by
   [`h5sparse`](https://pypi.org/project/h5sparse/) in either CSC or
   [compressed sparse row (CSR)](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format))
   format.

3. MATLAB save files.
   Requires explicitly loading [`MAT.jl`](https://github.com/JuliaIO/MAT.jl) first.

4. Julia sparse matrices saved to an HDF5 file in the `JLD` format.
   Requires explicitly loading [`JLD.jl`](https://github.com/JuliaIO/JLD.jl) first.

5. Julia sparse matrices saved to an HDF5 file in the `JLD2` format.
   Requires explicitly loading [`JLD2.jl`](https://github.com/JuliaIO/JLD.jl) first.

For repeated computation, the first format should be preferred due to its support for
[memory mapping](@ref obsmat_mmap) the matrix from disk.[^1]

[^1]:
    In particular, it is not possible to memory map the on-disk arrays saved by `numpy`
    and MATLAB because the internal data structures use 0-indexed pointers which are
    incompatible with Julia's use of 1-indexing.

The return is a `NamedTuple` with the following fields:
- `R` — the `SparseMatrixCSC` observing matrix ``\mat R``.
- `pixr` — a description of the "right-hand side" map pixelization, i.e. the pixelization of
  a map-vector ``\mat{\hat m}`` for which the product ``\mat R \mat{\hat m}`` is defined.
- `pixl` — a description of the "left-hand side" map pixelization, i.e. the pixelization of
  the resultant map-vector ``\mat{\tilde m} = \mat R \mat{\hat m}``.

The pixelization descriptions are in general data-format and situation specific.
For the HDF5 backend, the pixelization description may be any "native" `HDF5.jl` data
type (such as numbers, strings, or arrays thereof) or an HDF5 group which is loaded
recursively into a `Dict{String,Any}` with leaf-nodes which are a native type.

## Exporting an Observing Matrix

To support [memory mapping](@ref obsmat_mmap) the sparse matrix from disk, the
[`write_obsmat`](@ref) function is provided to ensure proper data alignment and choice of
HDF5 options.

## [Using Memory Mapping](@id obsmat_mmap)

By default with the HDF5 backend, the component arrays of the sparse matrix are memory
mapped.
This choice is made for one primary reason — memory mapping the arrays allows multiple
instances of Julia to share a single copy of the matrix in RAM, with the OS managing this
behavior automatically.
This is especially useful for computations being run on a cluster where the cluster
manager may launch multiple jobs doing similar calculations on the same compute node.

The default memory mapping behavior can be changed by setting the values of
[`CMB.Files.READ_OBSMAT_MMAP`](@ref) and [`CMB.Files.READ_OBSMAT_MMAP_FLAGS`](@ref) —
the former can be set to `false` to disable memory mapping completely:
```julia
julia> CMB.Files.READ_OBSMAT_MMAP[] = false  # disables default memory mapping
```
while the latter dictates the flags which are used during the underlying `mmap` call.
On Linux systems, the default is
```julia
julia> CMB.Files.READ_OBSMAT_MMAP_FLAGS[] = UnixMmap.MMAP_SHARED | UnixMmap.MMAP_POPULATE
```
so that the arrays are pre-faulted into active RAM.
(On other OSs, the arrays can either be manually faulted with a read loop, or the loads
can be left to happen on demand when the matrix is used.)
