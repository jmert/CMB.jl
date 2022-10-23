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

## Pixelizations
```@autodocs
Modules = [CMB.Pixelizations]
Private = false
```

## File I/O
```@docs
CMB.Files.read_obsmat
CMB.Files.READ_OBSMAT_MMAP
CMB.Files.READ_OBSMAT_MMAP_FLAGS
CMB.Files.write_obsmat
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
CMB.PolarizationConventions.Convention
CMB.StokesCrossFields.Field
CMB.StokesCrossFields.TPol
CMB.StokesCrossFields.Pol
```
