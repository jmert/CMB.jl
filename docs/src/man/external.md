# Re-exported Packages

Some features are simply re-exported functions from other, external packages.
These packages have (typically) been extracted from this package since they can be used
in a self-contained manner, independently of their application to CMB analysis.

For more information, see the documentation provided by each of the following packages:

## [AssociatedLegendrePolynomials.jl](https://github.com/jmert/AssociatedLegendrePolynomials.jl)

The associated Legendre polynomials (ALPs)constructs part of the eigenfunction basis for
functions on the sphere and plays a very important role in analysis of signals on the
sphere.
In addition, the ALPs are required for the computation of
[pixel-pixel covariance matrices](@ref man_pixelcov) that this package uniquely provides.

Note that this package re-exports the `AssociatedLegendrePolynomials` module under the
name `Legendre` for both backwards compatibility with earlier versions of this package
as well as for brevity.

## [SphericalHarmonics.jl](https://github.com/jmert/SphericalHarmonicTransforms.jl)

The spherical harmonic transforms are a key analytical techniques used in analysis of the
CMB.
They are built on top of the associated Legendre polynomials.

## [Healpix.jl](https://github.com/jmert/Healpix.jl)

The HEALPix pixelization is a very common mapping format used in CMB (and other fields of
astronomy more generally).
For instance, this package provides the necessary glue to perform spherical harmonic
transforms on a (scalar, i.e. temperature field) map in the HEALPix format.
