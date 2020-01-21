# CMB.jl Documentation

`CMB.jl` is a library of routines for the analysis of cosmic microwave
background (CMB) data. Development of features is being driven by the author's
use cases — at this time, namely the production of “reobserved” pixel-pixel
covariance matrices as used by the BICEP/Keck Array collaboration.

Design goals of this package include:

  * Native Julia implementation of core routines.

  * Numerical stability and efficiency.

  * Parallelism and efficient memory sharing.

## User Manual and Documentation
```@contents
Pages = [
    "man/healpix.md",
    "man/legendre.md",
    "man/sphere.md",
    "man/references.md"
]
Depth = 1
```

## Library API Reference
```@contents
Pages = [
    "lib/public.md",
    "lib/private.md"
]
Depth = 1
```

### [Index](@id main-index)
```@index
Pages = ["lib/public.md"]
```

