# CMB.jl — CMB Analysis

| **Documentation**                                                         | **Build Status**                                             |
|:-------------------------------------------------------------------------:|:------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][ci-img]][ci-url][![][codecov-img]][codecov-url] |

`CMB.jl` is a library of routines for the analysis of cosmic microwave
background (CMB) data. Development of features is being driven by the author's
use cases — at this time, namely the production of “reobserved” pixel-pixel
covariance matrices as used by the BICEP/Keck Array collaboration.

Design goals of this package include:

  * Native Julia implementation of core routines.

  * Numerical stability and efficiency.

  * Parallelism and efficient memory sharing.

### Installation and usage

This library is **not** registered in Julia's [General registry][General.jl],
so the package must be installed either by cloning it directly:

```
(@v1.7) pkg> add https://github.com/jmert/CMB.jl
```

or by making use of my [personal registry][Registry.jl]:

```
(@v1.7) pkg> registry add https://github.com/jmert/Registry.jl
(@v1.7) pkg> add CMB
```

After installing, just load like any other Julia package:

```
julia> using CMB
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://jmert.github.io/CMB.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://jmert.github.io/CMB.jl/dev

[ci-img]: https://github.com/jmert/CMB.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/jmert/CMB.jl/actions

[codecov-img]: https://codecov.io/gh/jmert/CMB.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jmert/CMB.jl

[General.jl]: https://github.com/JuliaRegistries/General
[Registry.jl]: https://github.com/jmert/Registry.jl
