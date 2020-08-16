# Public Documentation

## Contents
```@contents
Pages = ["public.md"]
```

## Sphere Functions
```@autodocs
Modules = [CMB.Sphere]
Private = false
```

## Healpix
```@autodocs
Modules = [CMB.Healpix]
Private = false
```

## Pixel Covariance
```@autodocs
Modules = [CMB.PixelCovariance]
Private = false
```

The following enums/bit flags are not exported globally, but all of the named values can
be imported into a scope by `using` the parent module (e.g. to access all of the
covariance field constants, use `using CMB.PixelCovariance.CovarianceFields`.

```@docs
CMB.PixelCovariance.PolarizationConventions.Convention
CMB.PixelCovariance.CovarianceFields.Field
CMB.PixelCovariance.CovarianceFields.TPol
CMB.PixelCovariance.CovarianceFields.Pol
```
