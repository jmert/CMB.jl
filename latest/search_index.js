var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "CMB.jl Documentation",
    "title": "CMB.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#CMB.jl-Documentation-1",
    "page": "CMB.jl Documentation",
    "title": "CMB.jl Documentation",
    "category": "section",
    "text": "CMB.jl is a library of routines for the analysis of cosmic microwave background (CMB) data. Development of features is being driven by the author's use cases — at this time, namely the production of “reobserved” pixel-pixel covariance matrices as used by the BICEP/Keck Array collaboration.Design goals of this package include:Native Julia implementation of core routines.\nNumerical stability and efficiency.\nParallelism and efficient memory sharing."
},

{
    "location": "index.html#User-Manual-and-Documentation-1",
    "page": "CMB.jl Documentation",
    "title": "User Manual and Documentation",
    "category": "section",
    "text": "Pages = [\n    \"man/legendre.md\",\n    \"man/references.md\"\n]\nDepth = 1"
},

{
    "location": "index.html#Library-API-Reference-1",
    "page": "CMB.jl Documentation",
    "title": "Library API Reference",
    "category": "section",
    "text": "Pages = [\n    \"lib/public.md\",\n    \"lib/private.md\"\n]\nDepth = 1"
},

{
    "location": "index.html#main-index-1",
    "page": "CMB.jl Documentation",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"lib/public.md\"]"
},

{
    "location": "man/healpix.html#",
    "page": "HEALPix Pixelization",
    "title": "HEALPix Pixelization",
    "category": "page",
    "text": ""
},

{
    "location": "man/healpix.html#HEALPix-Pixelization-1",
    "page": "HEALPix Pixelization",
    "title": "HEALPix Pixelization",
    "category": "section",
    "text": "Pages = [\"healpix.md\"]\nDepth = 2The Healpix module implements a selection of functions for interacting with the Hierarchical Equal-Area Iso-Latitude Pixelization (HEALPix) as described in Górski et al. (2005). The emphasis has been on working with pixel index and spherical coordinate conversions for the RING-ordered scheme only. If a complete HEALPix implementation is required, try LibHealpix.jl which provides Julia bindings to the official C++ library or try using the Python healpy interface via PyCall.jl. Additionally, the Healpix.jl package aims for a complete native Julia reimplementation. newcommandNsideN_mathrmside"
},

{
    "location": "man/healpix.html#Basic-pixelization-properties-1",
    "page": "HEALPix Pixelization",
    "title": "Basic pixelization properties",
    "category": "section",
    "text": "A core defining attribute of the HEALPix map is the Nside parameter. The minimum valid value is Nside = 1, representing the coarsest pixelization of the sphere which can be represented by HEALPix. The Nside then increases by factors of two—denoting ever finer resolutions—with the number of pixels covering the sphere being a function of Nside alone. CMB.Healpix.nside2npix returns the number of pixels in a given map:julia> using CMB.Healpix\n\njulia> nside = 4;\n\njulia> nside2npix(nside)\n192Because all pixels are of equal area and the number of pixels is derived from only the Nside, the pixel surface area must be as well. For convenience, CMB.Healpix.nside2pixarea is provided and is equivalent to the simple from-scratch calculation (in steradians).julia> 4π / nside2npix(nside)\n0.06544984694978735\n\njulia> nside2pixarea(nside)\n0.06544984694978735The number of iso-latitude rings is also a function of only Nside and calculated by CMB.Healpix.nside2nring:julia> nside2nring(nside)\n15"
},

{
    "location": "man/healpix.html#Working-with-pixel-indices-1",
    "page": "HEALPix Pixelization",
    "title": "Working with pixel indices",
    "category": "section",
    "text": "The pixels are enumerated as 0-indexed integers from west to east along the iso-latitude rings, from north to south. For example, pixel 0 is the first pixel in the first ring, and pixel 103 is the sixteenth pixel in the eighth ring for an Nside = 4 map:julia> pix = (0, 103);\n\njulia> pix2ring.(nside, pix)    # Ring\n(1, 8)\n\njulia> pix2ringidx.(nside, pix) # Index in ring\n(1, 16)note: Note\nBe careful to note that pixels are 0-indexed, but the rings and indices within a ring are 1-indexed.The HEALPix grid is symmetric about equator, with the equatorial ring considered part of the northern hemisphere by convention. Membership as part of the northern or southern hemisphere can be tested with the CMB.Healpix.isnorth and CMB.Healpix.issouth functions, respectively. Pixel 103 is actually the last pixel in the northern hemisphere, sojulia> isnorth(nside, 103)\ntrue\n\njulia> isnorth(nside, 104)\nfalse\n\njulia> issouth(nside, 104)\ntrueIn fact, each hemisphere is further composed of a so-called polar cap and equatorial belt region of pixels (a property derived from the mathematical details of the HEALPix grid's definition). According to the ring-ordered definition, pixel 0 should be in the polar cap (tested via CMB.Healpix.iscap), while pixel 103 in the equatorial ring is expected to be part of the equitorial belt (tested via CMB.Healpix.isequbelt).julia> iscap.(nside, pix)\n(true, false)\n\njulia> isequbelt.(nside, pix)\n(false, true)Membership in a particular hemisphere's polar cap or equatorial belt is accomplished with variants inserting north and south into the function names, i.e. polar caps are distinguished by CMB.Healpix.isnorthcap and CMB.Healpix.issouthcap, and the halves of the equatorial belt are distinguished by CMB.Healpix.isnorthequbelt and CMB.Healpix.issouthequbelt.julia> pix = (0, 103, 104, 191);\n\njulia> isnorthcap.(nside, pix)\n(true, false, false, false)\n\njulia> isnorthequbelt.(nside, pix)\n(false, true, false, false)\n\njulia> issouthequbelt.(nside, pix)\n(false, false, true, false)\n\njulia> issouthcap.(nside, pix)\n(false, false, false, true)"
},

{
    "location": "man/healpix.html#Working-with-spherical-coordinates-1",
    "page": "HEALPix Pixelization",
    "title": "Working with spherical coordinates",
    "category": "section",
    "text": "Up to now, all the features shown have concerned working with properties of the pixelization scheme, but the utility of the HEALPix grid is its ability to describe the surface of a sphere. Using spherical coordinates is more useful and more natural for more algorithms than the HEALPix-specific indexing scheme.The first method of describing the location of a particular HEALPix pixel is as a colatitude/azimuth pair of angles on the surface of the sphere identifying the pixel center. Colatitude measures the angle (in radians) south of the North Pole, and azimuth measures the angle (in radians) east of the Prime Meridian. To get the colatitude, use CMB.Healpix.pix2theta,julia> pix2theta(nside, 103)\n1.5707963267948966and to get the azimuth, use CMB.Healpix.pix2phijulia> pix2phi(nside, 103)\n6.086835766330224(both named to follow the mathematical convention that colatitude/azimuth pairs in spherical coordinates are the variable pair ( )). When the coordinate pair is required, the method CMB.Healpix.pix2ang returns a 2-tuple with the coordinates:julia> pix2ang(nside, 103)\n(1.5707963267948966, 6.086835766330224)\n\njulia> pix2ang(nside, 103) .|> rad2deg\n(90.0, 348.75)The other common way to represent coordinates on the sphere is via unit vectors. The corresponding vector for a given pixel is retrieved with CMB.Healpix.pix2vec.julia> pix2vec(nside, 103)\n3-element SVector{3,Float64}:\n  0.980785\n -0.19509\n  0.0where the elements correspond to the typical (x y z) right-handed coordinates with the positive z-axis passing through the North Pole and the positive x-axis passing through the Prime Meridian.In reverse, converting an arbitrary spherical coordinate to a pixel index... ...TO BE IMPLEMENTED..."
},

{
    "location": "man/healpix.html#Input-validation-and-error-handling-1",
    "page": "HEALPix Pixelization",
    "title": "Input validation and error handling",
    "category": "section",
    "text": "As stated earlier, the HEALPix Nside parameter takes on values which are powers of two and by convention of the official HEALPix [1] source is limited to the range 1 to 2^29. Validity of any nside parameter can be checked with the CMB.Healpix.ishealpixok function.julia> ishealpixok(4)\ntrue\n\njulia> ishealpixok.((5, 2^30))\n(false, false)Likewise, once given an Nside value, the pixel indices are bounded in 0 to nside2npix(nside) - 1; the two-argument form of ishealpixok returns whether a pixel is valid for the specified nside or not:julia> nside2npix(4)\n192\n\njulia> ishealpixok(4, 191)\ntrue\n\njulia> ishealpixok(4, 192)\nfalseVariants which throw a CMB.Healpix.InvalidNside or CMB.Healpix.InvalidPixel error are provided by CMB.Healpix.checkhealpix:julia> checkhealpix(5)\nERROR: 5 is not a valid Nside parameter (must be power of 2)\n[...]\n\njulia> checkhealpix(4, 192)\nERROR: 192 is not a valid pixel index for Nside = 4 (must be from 0 to 191)\n[...]note: Note\nOnly the functions working with spherical coordinates validate their inputs. The pixel indexing and classification functions are considered low-level routines and assume valid inputs. For example,julia> nside2npix(5)\n300\n\njulia> pix2ring(5, 0)\n1\n\njulia> pix2theta(5, 0)\nERROR: 5 is not a valid Nside parameter (must be power of 2)\n[...]\n\njulia> isnorth(4, -1)\ntrue\n\njulia> pix2ringidx(4, -1)\n0\n\njulia> pix2phi(4, -1)\nERROR: -1 is not a valid pixel index for Nside = 4 (must be from 0 to 191)\n[...]This choise was made for the sake of computational efficiency — the low-level pixel indexing/classification functions are used internally to compute the spherical coordinates."
},

{
    "location": "man/healpix.html#Footnotes:-1",
    "page": "HEALPix Pixelization",
    "title": "Footnotes:",
    "category": "section",
    "text": "[1]: Official HEALPix package: http://healpix.sourceforge.net/"
},

{
    "location": "man/legendre.html#",
    "page": "Legendre Polynomials",
    "title": "Legendre Polynomials",
    "category": "page",
    "text": ""
},

{
    "location": "man/legendre.html#Legendre-Polynomials-1",
    "page": "Legendre Polynomials",
    "title": "Legendre Polynomials",
    "category": "section",
    "text": "Pages = [\"legendre.md\"]\nDepth = 2The Legendre module implementation has been largely based on the approach of Limpanuparb and Milthorpe (2014)."
},

{
    "location": "man/legendre.html#legendre_defn-1",
    "page": "Legendre Polynomials",
    "title": "Definition and Properties",
    "category": "section",
    "text": "The associated Legendre polynomials P_^m(x) are the solution to the differential equationbeginalign\n    (1-x^2) fracd^2dx^2P_^m(x) - 2x fracddxP_^m(x) + left (+1) -\n        fracm^21-x^2 right P_^m(x) = 0\nendalignwhich arises as the colatitude  part of solving Laplace's equation ^2  +  = 0 in spherical coordinates (where x = cos()).There are several different conventions used to define P_^m that provide different properties, but the convention used here is typical of quantum mechanics and obeys the following properties:Solutions only exist for integer  and m, where   0 and m  .\nThe associated Legendre functions are normalized such that P_0^0 is unity and have orthogonality conditions,\nbeginalign\n    int_-1^1 P_^m(x) P_^m(x)mathrmdx\n        = frac22+1 frac(+m)(-m)\n        delta_\nendalign\nfor constant m and\nbeginalign\n    int_-1^1 fracP_^m(x) P_^m(x)1-x^2mathrmdx\n        = frac1m frac(+m)(-m) delta_mm\nendalign\nfor constant , where  is the Kronecker delta.\nThe phase convention for the Legendre functions is chosen such that the negative orders are related to positive orders according to,\nbeginalign\n    P_^-m(x) = (-1)^m frac(-m)(+m) P_^m(x)\nendalign\nThe Legendre functions can be enumerated for non-negative m using the three following recursion relations (given the initial condition P_0^0(x)):\nbeginalign\n    ( - m + 1)P_+1^m(x) = (2+1)xP_^m(x) - (+m)P_-1^m(x)\n    labeleqnstd_rr_2term\n    \n    P_+1^+1(x) = -(2+1)sqrt1-x^2 P_^(x)\n    labeleqnstd_rr_1term_lm\n    \n    P_+1^(x) = x(2+1)P_^(x)\n    labeleqnstd_rr_1term_l\nendalign"
},

{
    "location": "man/legendre.html#legendre_usage-1",
    "page": "Legendre Polynomials",
    "title": "Usage",
    "category": "section",
    "text": ""
},

{
    "location": "man/legendre.html#Calculating-scalar-values-1",
    "page": "Legendre Polynomials",
    "title": "Calculating scalar values",
    "category": "section",
    "text": "At the most basic, the associated Legendre polynomial P_^m(x) is computed by calling CMB.Legendre.Plm. For example, to compute P_2^1(05),julia> using CMB.Legendre\n\njulia> Plm(2, 1, 0.5)\n-1.299038105676658When m = 0 and only the Legendre polynomial P_(x) is needed, CMB.Legendre.Pl can be used instead:julia> Plm(2, 0, 0.5)\n-0.125\n\njulia> Pl(2, 0.5)\n-0.125\n\njulia> Pl(2, 0.5) == Plm(2, 0, 0.5)\ntrueIn the context of CMB analysis, a common use of the associated Legendre polynomials is to compute the spherical harmonics Y_m():beginalign\n    beginaligned\n    Y_m()  (-1)^m N_^m P_^m(cos ) e^im \n    textwhere  N_^m  sqrtfrac2+14 frac(-m)(+m)\n    endaligned\nendalignThe function CMB.Legendre.Nlm calculates the normalization factor N_^m:julia> Nlm(2, 0)\n0.6307831305050401\n\njulia> Nlm(2, 0) * Plm(2, 0, 0.5)\n-0.07884789131313001An important fact about the associated Legendre polynomials is that for m  0, P_^m(x) diverges to  as    [1]. For even moderately large pairs of (m), numerical underflow and overflow make computing the spherical harmonics impossible this way:julia> n = Nlm(157, 150)      # Underflows\n0.0\n\njulia> p = Plm(157, 150, 0.5) # Overflows\nInf\n\njulia> n * p                  # Undefined\nNaNOne way around this would be to just use extended precision arithmeticjulia> n = Nlm(BigFloat, 157, 150)\n4.14800666209481424285411223457923933542541063872695815968861285171699012214351e-314\n\njulia> p = Plm(157, 150, big\"0.5\")\n4.768286486602206390406601862422168575170463348990958242752608686436785229641202e+308\n\njulia> Float64(n * p)\n1.9778884113202627e-5but at the expense of much more computationally expensive calculations.An alternative way forward is to directly calculate the spherical harmonic normalized associated Legendre polynomials _^m(x) so that the spherical harmonics are defined asbeginalign\n    beginaligned\n    Y_m() = (-1)^m _^m(cos ) e^im \n     textwhere  _^m(x)  N_^m P_^m(x)\n    endaligned\nendalignCMB.Legendre.λlm implements this scheme and avoids the under/overflow of computing the normalization separately from the function:julia> λlm(157, 150, 0.5)\n1.977888411320241e-5note: Note\nWe are not just limited to efficient and numerically stable computation of _^m(x); the package supports arbitrary normalizations.  For further information on implementing custom Legendre normalizations, see the Custom normalizations section."
},

{
    "location": "man/legendre.html#Calculating-all-values-up-to-a-given-ℓ_\\mathrm{max}-1",
    "page": "Legendre Polynomials",
    "title": "Calculating all values up to a given _mathrmmax",
    "category": "section",
    "text": "Because calculating a particular Legendre polynomial value is the end result of running a recurrence relation, using Julia's dot broadcasting to compute P_^m(x) for all  is inefficient and redoes a lot of work:julia> λ = zeros(701);\n\njulia> @time λ[3:701] .= λlm.(2:700, 2, 0.5);\n  0.042107 seconds (4.61 k allocations: 257.940 KiB)It's far more efficient to incrementally calculate the +1 term directly from the  term. Both of Plm and λlm have modifying counterparts, Plm! and λlm! respectively, which fill an appropriately sized vector for a specified _mathrmmax.julia> @time λlm!(λ, 700, 2, 0.5);\n  0.000036 seconds (4 allocations: 160 bytes)On my machine, this ends up being roughly 1000 times faster!If all Legendre polynomial values for some x over all   0_mathrmmax and m  0 are required, there are also methods of Plm! and λlm! which fill the entire [lower triangular] matrix of values:julia> Λ = zeros(701, 701);\n\njulia> λlm!(Λ, 700, 0.5);\n\njulia> Λ[701,3] == λlm(700, 2, 0.5)   # N.B. 1-based indexing of the array!\ntrue"
},

{
    "location": "man/legendre.html#Precomputed-recursion-factors-1",
    "page": "Legendre Polynomials",
    "title": "Precomputed recursion factors",
    "category": "section",
    "text": "A final trick to accelerating calculation of any normalization of the associated Legendre polynomials is to pre-compute the appropriate recursion relation coefficients.At a low level, Plm/Plm! and λlm/λlm! are simple wrappers around the general legendre/legendre! functions. The trait type LegendreUnitNorm dispatches internal functions to compute P_^m(x):julia> legendre(LegendreUnitNorm(), 5, 2, 0.5) == Plm(5, 2, 0.5)\ntrueand LegendreSphereNorm does the same for _^m(x):julia> legendre(LegendreSphereNorm(), 5, 2, 0.5) == λlm(5, 2, 0.5)\ntrueThe type LegendreNormCoeff stores the coefficients for a particular normalization (and value type) so that the coefficients must only be calculated once.julia> coeff = LegendreNormCoeff{LegendreSphereNorm, Float64}(700);\n\njulia> legendre(coeff, 5, 2, 0.5)\n-0.15888479843070935On my machine, this results in a further ~50% decrease in computation time compared to λlm!:julia> @time legendre!(coeff, λ, 700, 2, 0.5);\n  0.000020 seconds (4 allocations: 160 bytes)"
},

{
    "location": "man/legendre.html#legendre_customnorm-1",
    "page": "Legendre Polynomials",
    "title": "Custom normalizations",
    "category": "section",
    "text": "CMB.Legendre provides the standard and spherical harmonic normalizations by default, but arbitrary normalizations are also supported. The mile-high overview is that the initial condition and recurrence relation (r.r.) coefficients are all methods which dispatch on a normalization trait type, so a new normalization is added by simply extending appropriate types and methods. The following table lists all of the types to extend and method specialization to implement.Interfaces to extend/implement Brief description\nCMB.Legendre.AbstractLegendreNorm Supertype of normalization trait types\nCMB.Legendre.Plm_00() Value of N_0^0 P_0^0(x) for the given normalization\nCMB.Legendre.Plm_μ() Coefficient for the 1-term r.r. boosting   +1 and m  m+1\nCMB.Legendre.Plm_ν() Coefficient for the 1-term r.r. boosting   +1\nCMB.Legendre.Plm_α() Coefficient for the 2-term r.r. acting on the (m) term\nCMB.Legendre.Plm_β() Coefficient for the 2-term r.r. acting on the (-1m) termAs a concrete example, we'll walk through how _^m(x) is defined to have the spherical harmonic normalization baked in.beginalign\n    _^m(x)  N_^m P_^m(x)\n    \n    N_^m = sqrtfrac2+14 frac(-m)(+m)\nendalignBaking in the normalization happens by changing the coefficients in the recursion relations given in the Definitions and Properties section. For our purposes, they take on the form:beginalign\n    P_ell+1^m(x) = alpha_ell+1^m x P_ell^m(x)\n        - beta_ell+1^m P_ell-1^m(x)\n        labeleqncus_rr_2term\n    \n    P_ell+1^ell+1(x) = mu_ell+1 sqrt1-x^2\n        P_ell^ell(x)\n        labeleqncus_rr_1term_lm\n    \n    P_ell+1^ell(x) = nu_ell+1 x P_ell^ell(x)\n        labeleqncus_rr_1term_l\nendalignThe normalization is encoded in the coefficients _^m, _^m, _, and _. For the standard (unity) normalization, these take on the valuesbeginalign\n    _^m = frac2 - 1 - m \n    _^m = frac + m - 1 - m \n    _ = 2 - 1 \n    _ = 2 - 1\nendalignby simply identifying the coefficients from Eqns. refeqnstd_rr_2term–refeqnstd_rr_1term_l on each of the P_^m(x) terms on the right hand side. For other normalizations, we multiply through by the normalization factor appropriate for the left-hand side of the equations, rearrange terms to correctly normalize the terms on the right, and identify the coefficients left over. For example, _^m and _^m for _^m(x) are determined by starting with Eq. refeqnstd_rr_2term and multiply through by N_+1^m. The left-hand side by definition is _+1^m, leaving us withbeginalign\n    beginsplit\n        _+1^m = frac2 + 1 - m + 1 x\n            sqrtfrac2+34 frac(-m+1)(+m+1) P_^m(x) -\n            \n            quadquad frac+m-m+1 sqrtfrac2+34\n            frac(-m+1)(+m+1) P_-1^m(x)\n    endsplit\nendalignThrough judicious use of algebra, the terms on the right-hand side can be manipulated to gather terms of the form N_^m P_^m(x) and N_-1^m P_-1^m(x), leaving us withbeginalign\n    _+1^m = sqrtfrac2+32-1 frac4^2 - 1(+1)^2 - m^2 x\n        _^m(x) -\n        sqrtfrac2+32-1 frac^2 - m^2(+1)^2 - m^2\n        _-1^m(x)\nendalignWe identify each of the two square root terms as _+1^m and _+1^m since they are the cofficients appropriate for generating _+1^m(x). Doing so with the other two recurrence relation equations, we obtain:beginalign\n    _^m = sqrtfrac2+12-3 frac4(-1)^2 - 1^2 - m^2 \n    _^m = sqrtfrac2+12-3 frac(-1)^2 - m^2^2 - m^2 \n    _ = sqrt1 + frac12 \n    _ = sqrt2 + 1\nendalignThe final math required is to define the initial condition _0^0(x). This is straight forward given the definition:beginalign\n    _0^0(x) = N_0^0 P_0^0(x) = sqrtfrac14  1 \n    _0^0(x) = sqrtfrac14\nendalignWe now have all the information required to define a custom Legendre normalization. Begin by importing the types and methods which will need to be extended:julia> using CMB.Legendre\n\njulia> import CMB.Legendre: AbstractLegendreNorm, Plm_00, Plm_μ, Plm_ν, Plm_α, Plm_βWe'll call our new normalization λNorm, which must be a subclass of AbstractLegendreNorm.julia> struct λNorm <: AbstractLegendreNorm endThe initial condition is specified by providing a method of Plm_00 which takes our normalization trait type as the first argument. (The second argument can be useful if some extra type information is required to set up a type-stable algorithm.)julia> Plm_00(::λNorm, T::Type) = sqrt(1 / 4π)\nPlm_00 (generic function with 4 methods)Finally, we provide methods which encode the cofficients as well:julia> function Plm_α(::λNorm, T::Type, l::Integer, m::Integer)\n           fac1 = (2l + 1) / ((2l - 3) * (l^2 - m^2))\n           fac2 = 4*(l-1)^2 - 1\n           return sqrt(fac1 * fac2)\n       end\nPlm_α (generic function with 4 methods)\n\njulia> function Plm_β(::λNorm, T::Type, l::Integer, m::Integer)\n           fac1 = (2l + 1) / ((2l - 3) * (l^2 - m^2))\n           fac2 = (l-1)^2 - m^2\n           return sqrt(fac1 * fac2)\n       end\nPlm_β (generic function with 4 methods)\n\njulia> Plm_μ(::λNorm, T::Type, l::Integer) = sqrt(1 + 1 / 2l)\nPlm_μ (generic function with 4 methods)\n\njulia> Plm_ν(::λNorm, T::Type, l::Integer) = sqrt(1 + 2l)\nPlm_ν (generic function with 4 methods)With just those 5 methods provided, the full Legendre framework is available, including precomputing the coefficients.julia> legendre(λNorm(), 700, 500, 0.4)\n0.35366224602810997\n\njulia> coeff = LegendreNormCoeff{λNorm,Float64}(700);\n\njulia> legendre(coeff, 700, 500, 0.4)\n0.35366224602810997"
},

{
    "location": "man/legendre.html#Footnotes-1",
    "page": "Legendre Polynomials",
    "title": "Footnotes",
    "category": "section",
    "text": "[1]: Specifically, the envelope of P_^m(x) which bounds the local extrema for all values of x can be shown to be    left P_^m(cos ) right  frac(+m+1)(+frac32)\n        left( frac2 sin  right)^12(see Eq. 8.10.7 (p336) of Abramowitz and Stegun, “Handbook of Mathematical Functions” 10th printing (1972)). For fixed m and any x, we take the asymptotic limit as    and simplify (z) via Stirling's approximation to get the scaling of the associated Legendre polynomial envelope    DeclareMathOperator*envenv\n    env_left( P_^m right)  ^m - 12 text In contrast, the normalization factor N_^m scales as ^12 - m, exactly canceling the scaling of envleft(P_^mright), so overall the spherical harmonic normalized Legendre polynomials _^m(x) asymptote to some constant envelope:    env_ left( _^m right)  ^0 = textconstant "
},

{
    "location": "man/references.html#",
    "page": "References",
    "title": "References",
    "category": "page",
    "text": ""
},

{
    "location": "man/references.html#References-1",
    "page": "References",
    "title": "References",
    "category": "section",
    "text": ""
},

{
    "location": "man/references.html#bib-coordinates-1",
    "page": "References",
    "title": "Coordinate systems and conventions",
    "category": "section",
    "text": "J. P. Hamaker and J. D. Bregman. “Understanding radio polarimetry. III. Interpreting the IAU/IEEE definitions of the Stokes parameters.” In: Astronomy and Astrophysics, Supplement Series 117 (May 1996) pp. 161–165. DOI: 10.1051/aas:1996147"
},

{
    "location": "man/references.html#bib-legendre-1",
    "page": "References",
    "title": "Legendre functions",
    "category": "section",
    "text": "T. Limpanuparb and J. Milthorpe. “Associated Legendre Polynomials and Spherical Harmonics Computation for Chemistry Applications” In: Proceedings of the 40th Congress on Science and Technology of Thailand (Dec 2014) arXiv: 1410.1748\n“Legendre polynomials” on Wikipedia and Wolfram Math World\n“Associated Legendre polynomials” on Wikipedia and Wolfrm Math World"
},

{
    "location": "man/references.html#bib-healpix-1",
    "page": "References",
    "title": "HEALPix",
    "category": "section",
    "text": "K. M. Górski et al. “HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data Distributed on the Sphere” In: The Astrophysical Journal 662, 2 (Apr 2005) p. 759–771 DOI: 10.1086/427976. arXiv: astro-ph/0409513\nErratum: The equatorial-belt longitude phi equation (Eq 9) should have the +1 outside the modulus:\n  s = (i - N_mathrmside) operatornamemod 2 + 1"
},

{
    "location": "man/references.html#bib-pixelcovariance-1",
    "page": "References",
    "title": "Pixel Covariance",
    "category": "section",
    "text": "M. Tegmark and A. de Oliveira-Costa. “How to measure CMB polarization power spectra without losing information” In: Physical Review D 64, 063001 (Sep 2001) DOI: 10.1103/PhysRevD.64.063001 arXiv: astro-ph/0012120\nThe BICEP/Keck Array Collaboration. “BICEP/Keck Array VII: Matrix Based E/B Separation Applied to BICEP2 and the Keck Array” In: The Astrophysical Journal 825, 66 (Jul 2016) DOI: 10.3847/0004-637X/825/1/66 arXiv: 1603.05976"
},

{
    "location": "lib/public.html#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": "DocTestSetup = quote\n    using CMB\nend"
},

{
    "location": "lib/public.html#Public-Documentation-1",
    "page": "Public",
    "title": "Public Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "lib/public.html#Contents-1",
    "page": "Public",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "lib/public.html#CMB.Sphere",
    "page": "Public",
    "title": "CMB.Sphere",
    "category": "Module",
    "text": "Collection of routines for working with coordinates on the sphere.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.bearing",
    "page": "Public",
    "title": "CMB.Sphere.bearing",
    "category": "Function",
    "text": "Calculates the bearing angle (), defined as the angle between the meridian (at the first coordinate) and the great circle connecting the first coordinate to the second. Angles are measured clockwise and will be in the range 0). See also bearing2.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.bearing-NTuple{4,Any}",
    "page": "Public",
    "title": "CMB.Sphere.bearing",
    "category": "Method",
    "text": "α = bearing(θ₁, ϕ₁, θ₂, ϕ₂)\n\nPoints on the sphere are given as coordinate pairs () and () where  is the colatitude angle from the North Pole and  is the azimuthal angle, both in radians.\n\nExamples\n\njulia> bearing(π/2, 0.0, π/4, π/4)\n0.6154797086703871\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.bearing-Tuple{AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Public",
    "title": "CMB.Sphere.bearing",
    "category": "Method",
    "text": "α = bearing(r₁, r₂)\n\nPoints on the sphere are given as unit vectors r and r.\n\nExamples\n\njulia> bearing([1.0, 0.0, 0.0], [0.5, 0.5, sqrt(2)/2])\n0.6154797086703873\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.bearing2",
    "page": "Public",
    "title": "CMB.Sphere.bearing2",
    "category": "Function",
    "text": "Calculates the latitude/longitude vector components of the bearing angle (i.e.  = cos()  = sin()), defined as the angle between the meridian (at the first coordinate) and the great circle connecting the first coordinate to the second. See also bearing.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.bearing2-NTuple{4,Any}",
    "page": "Public",
    "title": "CMB.Sphere.bearing2",
    "category": "Method",
    "text": "(δθ, δϕ) = bearing2(θ₁, ϕ₁, θ₂, ϕ₂)\n\nPoints on the sphere are given as coordinate pairs () and () where  is the colatitude angle from the North Pole and  is the azimuthal angle, both in radians.\n\nExamples\n\njulia> bearing2(π/2, 0.0, π/4, π/4)\n(0.8164965809277261, 0.5773502691896256)\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.bearing2-Tuple{AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Public",
    "title": "CMB.Sphere.bearing2",
    "category": "Method",
    "text": "(δθ, δϕ) = bearing2(r₁, r₂)\n\nPoints on the sphere are given as unit vectors r and r.\n\nExamples\n\njulia> bearing2([1.0, 0.0, 0.0], [0.5, 0.5, sqrt(2)/2])\n(0.816496580927726, 0.5773502691896257)\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.cosdistance",
    "page": "Public",
    "title": "CMB.Sphere.cosdistance",
    "category": "Function",
    "text": "Calculates the cosine of the inner angle (z) between unit vectors pointing from the center of the sphere to two points on its surface.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.cosdistance-NTuple{4,Any}",
    "page": "Public",
    "title": "CMB.Sphere.cosdistance",
    "category": "Method",
    "text": "z = cosdistance(θ₁, ϕ₁, θ₂, ϕ₂)\n\nPoints on the sphere are given as coordinate pairs () and () where  is the colatitude angle from the North Pole and  is the azimuthal angle, both in radians.\n\nExamples\n\njulia> cosdistance(π/2, 0.0, π/4, π/4)\n0.5\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.cosdistance-Tuple{AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Public",
    "title": "CMB.Sphere.cosdistance",
    "category": "Method",
    "text": "z = cosdistance(r₁, r₂)\n\nPoints on the sphere are given as unit vectors r and r.\n\nExamples\n\njulia> cosdistance([1.,0.,0.], [0.5,0.5,sqrt(2)/2])\n0.5\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.distance",
    "page": "Public",
    "title": "CMB.Sphere.distance",
    "category": "Function",
    "text": "Calculates the inner angle () between unit vectors pointing from the center of the sphere to two points on its surface.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.distance-NTuple{4,Any}",
    "page": "Public",
    "title": "CMB.Sphere.distance",
    "category": "Method",
    "text": "σ = distance(θ₁, ϕ₁, θ₂, ϕ₂)\n\nPoints on the sphere are given as coordinate pairs () and () where  is the colatitude angle from the North Pole and  is the azimuthal angle, both in radians.\n\nExamples\n\njulia> distance(π/2, 0.0, π/4, π/4)\n1.0471975511965979\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Sphere.distance-Tuple{AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Public",
    "title": "CMB.Sphere.distance",
    "category": "Method",
    "text": "σ = distance(r₁, r₂)\n\nPoints on the sphere are given as unit vectors r and r.\n\nExamples\n\njulia> distance([1.,0.,0.], [0.5,0.5,sqrt(2)/2])\n1.0471975511965979\n\n\n\n"
},

{
    "location": "lib/public.html#Sphere-Functions-1",
    "page": "Public",
    "title": "Sphere Functions",
    "category": "section",
    "text": "Modules = [CMB.Sphere]\nPrivate = false"
},

{
    "location": "lib/public.html#CMB.Legendre",
    "page": "Public",
    "title": "CMB.Legendre",
    "category": "Module",
    "text": "Collections of functions which compute the associated Legendre functions.\n\nBased on implementation described in Limpanuparb and Milthorpe (2014) “Associated Legendre Polynomials and Spherical Harmonics Computation for Chemistry Applications” arXiv:1410.1748v1\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.AbstractLegendreNorm",
    "page": "Public",
    "title": "CMB.Legendre.AbstractLegendreNorm",
    "category": "Type",
    "text": "abstract type AbstractLegendreNorm end\n\nAbstract supertype for normalization conditions of the Associated Legendre polynomials.\n\nExample\n\njulia> subtypes(AbstractLegendreNorm)\n3-element Array{Union{DataType, UnionAll},1}:\n CMB.Legendre.LegendreNormCoeff\n CMB.Legendre.LegendreSphereNorm\n CMB.Legendre.LegendreUnitNorm\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.LegendreNormCoeff",
    "page": "Public",
    "title": "CMB.Legendre.LegendreNormCoeff",
    "category": "Type",
    "text": "struct LegendreNormCoeff{N<:AbstractLegendreNorm,T<:Real} <: AbstractLegendreNorm\n\nPrecomputed recursion relation coefficients for the normalization N and value type T.\n\nExample\n\njulia> LegendreNormCoeff{LegendreSphereNorm,Float64}(1)\nCMB.Legendre.LegendreNormCoeff{CMB.Legendre.LegendreSphereNorm,Float64} for lmax = 1 with coefficients:\n    μ: [0.0, 1.22474]\n    ν: [0.0, 1.73205]\n    α: [0.0 0.0; 1.73205 0.0]\n    β: [0.0 0.0; -0.0 0.0]\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.LegendreSphereCoeff",
    "page": "Public",
    "title": "CMB.Legendre.LegendreSphereCoeff",
    "category": "Type",
    "text": "LegendreSphereCoeff{T}\n\nTable type of precomputed recursion relation coefficients for the spherical harmonic normalization. Alias for LegendreNormCoeff{LegendreSphereNorm,T}.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.LegendreSphereNorm",
    "page": "Public",
    "title": "CMB.Legendre.LegendreSphereNorm",
    "category": "Type",
    "text": "struct LegendreSphereNorm <: AbstractLegendreNorm end\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.LegendreUnitCoeff",
    "page": "Public",
    "title": "CMB.Legendre.LegendreUnitCoeff",
    "category": "Type",
    "text": "LegendreUnitCoeff{T}\n\nPrecomputed recursion relation coefficients for the standard unit normalization. Alias for LegendreNormCoeff{LegendreUnitNorm,T}.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.LegendreUnitNorm",
    "page": "Public",
    "title": "CMB.Legendre.LegendreUnitNorm",
    "category": "Type",
    "text": "struct LegendreUnitNorm <: AbstractLegendreNorm end\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.Nlm-Union{Tuple{T}, Tuple{Type{T},Integer,Integer}} where T<:Real",
    "page": "Public",
    "title": "CMB.Legendre.Nlm",
    "category": "Method",
    "text": "N = Nlm([T=Float64], l, m)\n\nComputes the normalization constant\n\n    N_^m  sqrtfrac2+14 frac(-m)(+m)\n\nwhich defines the Spherical Harmonic normalized functions _^m(x) in terms of the standard unit normalized P_^m(x)\n\n    _^m(x)  N_^m P_^m(x)\n\nusing numbers of type T.\n\nSee also Plm and λlm.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.Pl!-Union{Tuple{AbstractArray{T,1},Integer,T}, Tuple{T}} where T<:Real",
    "page": "Public",
    "title": "CMB.Legendre.Pl!",
    "category": "Method",
    "text": "Pl!(P::AbstractVector, lmax::Integer, x::Real)\n\nFills the vector P with the Legendre polynomial values P_(x) for all degrees 0 ≤ ℓ ≤ lmax at x.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.Pl-Union{Tuple{Integer,T}, Tuple{T}} where T<:Real",
    "page": "Public",
    "title": "CMB.Legendre.Pl",
    "category": "Method",
    "text": "p = Pl(l::Integer, x::Real)\n\nComputes the scalar value p = P_(x), where P_(x) is the Legendre polynomial of degree l at x.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.Plm!-Union{Tuple{AbstractArray{T,1},Integer,Integer,T}, Tuple{T}} where T<:Real",
    "page": "Public",
    "title": "CMB.Legendre.Plm!",
    "category": "Method",
    "text": "Plm!(P::AbstractVector, lmax::Integer, m::Integer, x::Real)\n\nFills the vector P with the Legendre polynomial values P_^m(x) for all degrees 0 ≤ ℓ ≤ lmax and constant order m at x.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.Plm!-Union{Tuple{AbstractArray{T,2},Integer,T}, Tuple{T}} where T<:Real",
    "page": "Public",
    "title": "CMB.Legendre.Plm!",
    "category": "Method",
    "text": "Plm!(P::AbstractMatrix, lmax::Integer, x::Real)\n\nFills the lower triangle of the matrix P with the associated Legendre polynomial values P_^m(x) for all degrees 0 ≤ ℓ ≤ lmax and all orders 0 ≤ m ≤ ℓ at x.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.Plm-Union{Tuple{Integer,Integer,T}, Tuple{T}} where T<:Real",
    "page": "Public",
    "title": "CMB.Legendre.Plm",
    "category": "Method",
    "text": "p = Plm(l::Integer, m::Integer, x::Real)\n\nComputes the scalar value p = P_^m(x), where P_^m(x) is the associated Legendre polynomial of degree l and order m at x.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.legendre!-Union{Tuple{N,AbstractArray{T,1},Integer,Integer,T}, Tuple{N}, Tuple{T}} where T<:Real where N<:CMB.Legendre.AbstractLegendreNorm",
    "page": "Public",
    "title": "CMB.Legendre.legendre!",
    "category": "Method",
    "text": "legendre!(norm::AbstractLegendreNorm, Λ::AbstractVecotr, lmax::Integer, m::Integer,\n          x::Real)\n\nFills the vector Λ with the pre-normalized Legendre polynomial values N_^m P_^m(x) for all degrees 0 ≤ ℓ ≤ lmax and constant order m at x, where N_^m is the normalization scheme norm.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.legendre!-Union{Tuple{N,AbstractArray{T,2},Integer,T}, Tuple{N}, Tuple{T}} where T<:Real where N<:CMB.Legendre.AbstractLegendreNorm",
    "page": "Public",
    "title": "CMB.Legendre.legendre!",
    "category": "Method",
    "text": "legendre!(norm::AbstractLegendreNorm, P::AbstractMatrix, lmax::Integer, x::Real)\n\nFills the matrix Λ with the pre-normalized Legendre polynomial values N_^m P_^m(x) for all degrees 0 ≤ ℓ ≤ lmax and all orders 0 ≤ m ≤ ℓ at x, where N_^m is the normalization scheme norm.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.legendre-Union{Tuple{N,Integer,Integer,T}, Tuple{N}, Tuple{T}} where T<:Real where N<:CMB.Legendre.AbstractLegendreNorm",
    "page": "Public",
    "title": "CMB.Legendre.legendre",
    "category": "Method",
    "text": "p = legendre(norm::AbstractLegendreNorm, l::Integer, m::Integer, x::Real)\n\nComputes the scalar value p = N_^m P_^m(x), where P_^m(x) is the associated Legendre polynomial of degree l and order m at x and N_^m the normalization scheme norm.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.legendre-Union{Tuple{N,Integer,T}, Tuple{N}, Tuple{T}} where T<:Real where N<:CMB.Legendre.AbstractLegendreNorm",
    "page": "Public",
    "title": "CMB.Legendre.legendre",
    "category": "Method",
    "text": "p = legendre(norm::AbstractLegendreNorm, l::Integer, x::Real)\n\nComputes the scalar value p = N_ P_(x), where P_(x) is the Legendre polynomial of degree l at x and N_ is the normalization scheme norm.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.λlm!-Union{Tuple{AbstractArray{T,1},Integer,Integer,T}, Tuple{T}} where T<:Real",
    "page": "Public",
    "title": "CMB.Legendre.λlm!",
    "category": "Method",
    "text": "λlm!(Λ::AbstractVector, lmax::Integer, m::Integer, x::Real)\n\nFills the vector Λ with the spherical harmonic normalized associated Legendre polynomial values _^m(x) for all degrees 0 ≤ ℓ ≤ lmax and constant order m at x.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.λlm!-Union{Tuple{AbstractArray{T,2},Integer,T}, Tuple{T}} where T<:Real",
    "page": "Public",
    "title": "CMB.Legendre.λlm!",
    "category": "Method",
    "text": "λlm!(Λ::AbstractMatrix, lmax::Integer, x::Real)\n\nFills the lower triangle of the matrix Λ with the spherical harmonic normalized associated Legendre polynomial values _^m(x) for all degrees 0 ≤ ℓ ≤ lmax and all orders 0 ≤ m ≤ ℓ at x.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Legendre.λlm-Union{Tuple{Integer,Integer,T}, Tuple{T}} where T<:Real",
    "page": "Public",
    "title": "CMB.Legendre.λlm",
    "category": "Method",
    "text": "λ = λlm(l::Integer, m::Integer, x::Real)\n\nComputes the scalar value  = _^m(x), where _^m(x) is the spherical-harmonic normalized associated Legendre polynomial of degree l and order m at x.\n\n\n\n"
},

{
    "location": "lib/public.html#Legendre-Functions-1",
    "page": "Public",
    "title": "Legendre Functions",
    "category": "section",
    "text": "Modules = [CMB.Legendre]\nPrivate = false"
},

{
    "location": "lib/public.html#CMB.Healpix",
    "page": "Public",
    "title": "CMB.Healpix",
    "category": "Module",
    "text": "A module of functions implementing function which interact with the HEALPix pixel definitions. In most cases, only the RING ordering functions are being provided.\n\nSee \"HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data Distributed on the Sphere\" Górski, Hivon, & Banday et al (2005) ApJ 622:759–771 arXiv: astro-ph/0409513\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.UNSEEN",
    "page": "Public",
    "title": "CMB.Healpix.UNSEEN",
    "category": "Constant",
    "text": "const UNSEEN = -1.6375e+30\n\nSpecial value recognized by the libhealpix/healpy routines as an unobserved/masked pixel.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.InvalidNside",
    "page": "Public",
    "title": "CMB.Healpix.InvalidNside",
    "category": "Type",
    "text": "InvalidNside(nside)\n\nAn invalid nside value was provided.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.InvalidPixel",
    "page": "Public",
    "title": "CMB.Healpix.InvalidPixel",
    "category": "Type",
    "text": "InvalidPixel(nside, pix)\n\nAn invalid pixel index pix was provided for the given nside.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.checkhealpix-Tuple{Any,Any}",
    "page": "Public",
    "title": "CMB.Healpix.checkhealpix",
    "category": "Method",
    "text": "checkhealpix(nside, pix)\n\nThrows an InvalidNside exception if nside is not a valid value or an InvalidPixel exception if pix is out of range for the given N_mathrmside.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.checkhealpix-Tuple{Any}",
    "page": "Public",
    "title": "CMB.Healpix.checkhealpix",
    "category": "Method",
    "text": "checkhealpix(nside)\n\nThrows an InvalidNside exception if nside is not a valid value.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.iscap",
    "page": "Public",
    "title": "CMB.Healpix.iscap",
    "category": "Function",
    "text": "iscap(nside, pix)\n\nTest for whether the given pixel pix is in either polar cap for an nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.isequbelt",
    "page": "Public",
    "title": "CMB.Healpix.isequbelt",
    "category": "Function",
    "text": "isequbelt(nside, pix)\n\nTest for whether the given pixel pix is in the equatorial belt for an nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.ishealpixok-Tuple{Any,Any}",
    "page": "Public",
    "title": "CMB.Healpix.ishealpixok",
    "category": "Method",
    "text": "isheapixok(nside, pix)\n\nReturns true if nside is valid and pix is in the range 0 to nside2npix(nside) - 1, otherwise false.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.ishealpixok-Tuple{Any}",
    "page": "Public",
    "title": "CMB.Healpix.ishealpixok",
    "category": "Method",
    "text": "ishealpixok(nside)\n\nReturns true if nside is a power of two in the range 1 to 2^29, otherwise false.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.isnorth",
    "page": "Public",
    "title": "CMB.Healpix.isnorth",
    "category": "Function",
    "text": "isnorth(nside, pix)\n\nTest for whether the given pixel pix is in the northern hemisphere (including the equatorial ring) for an nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.isnorthcap",
    "page": "Public",
    "title": "CMB.Healpix.isnorthcap",
    "category": "Function",
    "text": "isnorthcap(nside, pix)\n\nTest for whether the given pixel pix is in the northern polar cap for an nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.isnorthequbelt",
    "page": "Public",
    "title": "CMB.Healpix.isnorthequbelt",
    "category": "Function",
    "text": "isnorthequbelt(nside, pix)\n\nTest for whether the given pixel pix is in the northern equatorial belt (including the equatorial ring) for an nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.issouth",
    "page": "Public",
    "title": "CMB.Healpix.issouth",
    "category": "Function",
    "text": "issouth(nside, pix)\n\nTest for whether the given pixel pix is in the southern hemisphere (excluding the equatorial ring) for an nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.issouthcap",
    "page": "Public",
    "title": "CMB.Healpix.issouthcap",
    "category": "Function",
    "text": "issouthcap(nside, pix)\n\nTest for whether the given pixel pix is in the southern polar cap for an nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.issouthequbelt",
    "page": "Public",
    "title": "CMB.Healpix.issouthequbelt",
    "category": "Function",
    "text": "issouthequbelt(nside, pix)\n\nTest for whether the given pixel pix is in the southern equatorial belt (excluding the equatorial ring) for an nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.npix2nside-Tuple{Any}",
    "page": "Public",
    "title": "CMB.Healpix.npix2nside",
    "category": "Method",
    "text": "Nisde = npix2nside(npix)\n\nReturns the equivalent Nside corresponding to the number of pixels npix.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.nring2nside-Tuple{Any}",
    "page": "Public",
    "title": "CMB.Healpix.nring2nside",
    "category": "Method",
    "text": "Nside = nring2nside(nring)\n\nReturns the equivalent Nside corresponding to the number of iso-latitude rings nring.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.nside2npix",
    "page": "Public",
    "title": "CMB.Healpix.nside2npix",
    "category": "Function",
    "text": "Npix = nside2npix(nside)\n\nReturns the total number of pixels Npix in an nside HEALPix map. Note that HEALPix pixel indexing is 0-based, so valid pixel values are in the range 0 to Npix - 1.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.nside2npixcap",
    "page": "Public",
    "title": "CMB.Healpix.nside2npixcap",
    "category": "Function",
    "text": "Npix = nside2npixcap(nside)\n\nReturns the number of pixels Npix in the polar caps for the given nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.nside2npixequ",
    "page": "Public",
    "title": "CMB.Healpix.nside2npixequ",
    "category": "Function",
    "text": "Npix = nside2npixequ(nside)\n\nReturns the number of pixels Npix in the northern hemisphere, including the equatorial ring, for the given nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.nside2nring",
    "page": "Public",
    "title": "CMB.Healpix.nside2nring",
    "category": "Function",
    "text": "Nring = nside2nring(nside)\n\nReturns the number of iso-latitude rings Nring in the nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.nside2pixarea",
    "page": "Public",
    "title": "CMB.Healpix.nside2pixarea",
    "category": "Function",
    "text": "σ = nside2pixarea(nside)\n\nReturns the surface area σ (in steradians) of each pixel in the given nside HEALPix map.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.pix2ang-Tuple{Any,Any}",
    "page": "Public",
    "title": "CMB.Healpix.pix2ang",
    "category": "Method",
    "text": "(θ,ϕ) = pix2ang(nside, p)\n\nComputes the colatitude and azimuth pair (θ,ϕ) for the given pixel p.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.pix2phi-Union{Tuple{I,I}, Tuple{I}} where I<:Integer",
    "page": "Public",
    "title": "CMB.Healpix.pix2phi",
    "category": "Method",
    "text": "ϕ = pix2phi(nside, p)\n\nComputes the azimuth ϕ for the given pixel p. nside is the Nside resolution factor.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.pix2ring-Union{Tuple{I,I}, Tuple{I}} where I<:Integer",
    "page": "Public",
    "title": "CMB.Healpix.pix2ring",
    "category": "Method",
    "text": "i = pix2ring(nside, p)\n\nComputes the ring index i for the given pixel p. nside is the Nside resolution factor.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.pix2ringidx-Union{Tuple{I,I}, Tuple{I}} where I<:Integer",
    "page": "Public",
    "title": "CMB.Healpix.pix2ringidx",
    "category": "Method",
    "text": "j = pix2ringidx(nside, p)\n\nComputes the index j within the ring for the given pixel p. nside is the Nside resolution factor.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.pix2theta-Tuple{Any,Any}",
    "page": "Public",
    "title": "CMB.Healpix.pix2theta",
    "category": "Method",
    "text": "θ = pix2theta(nside, p)\n\nComputes the colatitude θ for the given pixel p. nside is the Nside resolution factor.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.pix2vec-Union{Tuple{I,I}, Tuple{I}} where I<:Integer",
    "page": "Public",
    "title": "CMB.Healpix.pix2vec",
    "category": "Method",
    "text": "r::SVector{3} = pix2vec(nside, p)\n\nComputes the unit vector r pointing to the pixel center of the given pixel p.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Healpix.pix2z-Union{Tuple{I,I}, Tuple{I}} where I<:Integer",
    "page": "Public",
    "title": "CMB.Healpix.pix2z",
    "category": "Method",
    "text": "z = pix2z(nside, p)\n\nComputes the cosine of the colatitude z for the given pixel p. nside is the Nside resolution factor.\n\n\n\n"
},

{
    "location": "lib/public.html#Healpix-1",
    "page": "Public",
    "title": "Healpix",
    "category": "section",
    "text": "Modules = [CMB.Healpix]\nPrivate = false"
},

{
    "location": "lib/public.html#CMB.PixelCovariance",
    "page": "Public",
    "title": "CMB.PixelCovariance",
    "category": "Module",
    "text": "Collection of functions which compute the pixel-pixel covariance of the CMB sky.\n\nBased on equations given in Tegmark and de Oliveira-Costa (2001) “How to measure CMB polarization power spectra without losing information” arXiv:astro-ph/0012120v3\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.PixelCovariance.PixelCovarianceCache",
    "page": "Public",
    "title": "CMB.PixelCovariance.PixelCovarianceCache",
    "category": "Type",
    "text": "struct PixelCovarianceCache\n\nData structure which contains all the information and buffers required to compute the pixel-pixel covariance terms for a given pixel with respect to all other pixels.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.PixelCovariance.PixelCovarianceCoeff",
    "page": "Public",
    "title": "CMB.PixelCovariance.PixelCovarianceCoeff",
    "category": "Type",
    "text": "struct PixelCovarianceCoeff{T<:Real}\n\nPrecomputed recursion relation coefficients for computing the pixel-pixel covariance.\n\nExample\n\njulia> PixelCovarianceCoeff{Float64}(2)\nCMB.PixelCovariance.PixelCovarianceCoeff{Float64} for lmax = 2 with coefficients:\n    λ: CMB.Legendre.LegendreNormCoeff{CMB.Legendre.LegendreUnitNorm,Float64}\n    η: [0.0795775, 0.238732, 0.397887]\n    α: [0.0, 0.0, 0.324874]\n    β: [0.0, 0.0, 0.0331573]\n\n\n\n"
},

{
    "location": "lib/public.html#Pixel-Covariance-1",
    "page": "Public",
    "title": "Pixel Covariance",
    "category": "section",
    "text": "Modules = [CMB.PixelCovariance]\nPrivate = false"
},

{
    "location": "lib/public.html#CMB.Util",
    "page": "Public",
    "title": "CMB.Util",
    "category": "Module",
    "text": "Miscellaneous utility functions.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Util.outer",
    "page": "Public",
    "title": "CMB.Util.outer",
    "category": "Function",
    "text": "Computes the outer product between a given column of a sparse matrix and a vector.\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Util.outer-Union{Tuple{AbstractArray{Tv,1},SparseMatrixCSC{Tv,Ti},Integer}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "Public",
    "title": "CMB.Util.outer",
    "category": "Method",
    "text": "outer(w::AbstractVector, A::SparseMatrixCSC, n::Integer)\n\nPerforms the equivalent of vec w veca_n^dagger where vec a_n is the column A[:,n].\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Util.outer-Union{Tuple{SparseMatrixCSC{Tv,Ti},Integer,AbstractArray{Tv,1}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "Public",
    "title": "CMB.Util.outer",
    "category": "Method",
    "text": "outer(A::SparseMatrixCSC, n::Integer, w::AbstractVector)\n\nPerforms the equivalent of vec a_n vec w^dagger where vec a_n is the column A[:,n].\n\n\n\n"
},

{
    "location": "lib/public.html#CMB.Util.quadprod",
    "page": "Public",
    "title": "CMB.Util.quadprod",
    "category": "Function",
    "text": "quadprod(A, b, n, dir=:col)\n\nComputes the quadratic product ABA^T efficiently for the case where B is all zero except for the nth column or row vector b, for dir = :col or dir = :row, respectively.\n\n\n\n"
},

{
    "location": "lib/public.html#Miscellaneous-Utilities-1",
    "page": "Public",
    "title": "Miscellaneous Utilities",
    "category": "section",
    "text": "Modules = [CMB.Util]\nPrivate = false"
},

{
    "location": "lib/private.html#",
    "page": "Private",
    "title": "Private",
    "category": "page",
    "text": "DocTestSetup = quote\n    using CMB\nend"
},

{
    "location": "lib/private.html#Private-Documentation-1",
    "page": "Private",
    "title": "Private Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "lib/private.html#Contents-1",
    "page": "Private",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"private.md\"]"
},

{
    "location": "lib/private.html#CMB.Sphere.x̂",
    "page": "Private",
    "title": "CMB.Sphere.x̂",
    "category": "Constant",
    "text": "const x̂ = SVector(1, 0, 0)\nconst ŷ = SVector(0, 1, 0)\nconst ẑ = SVector(0, 0, 1)\n\nConstant unit vectors in the Cartesian directions.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Sphere.ŷ",
    "page": "Private",
    "title": "CMB.Sphere.ŷ",
    "category": "Constant",
    "text": "const x̂ = SVector(1, 0, 0)\nconst ŷ = SVector(0, 1, 0)\nconst ẑ = SVector(0, 0, 1)\n\nConstant unit vectors in the Cartesian directions.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Sphere.ẑ",
    "page": "Private",
    "title": "CMB.Sphere.ẑ",
    "category": "Constant",
    "text": "const x̂ = SVector(1, 0, 0)\nconst ŷ = SVector(0, 1, 0)\nconst ẑ = SVector(0, 0, 1)\n\nConstant unit vectors in the Cartesian directions.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Sphere.:∥-Tuple{Any,Any}",
    "page": "Private",
    "title": "CMB.Sphere.:∥",
    "category": "Method",
    "text": "∥(u, v)\n\nTest whether vector u is parallel to vector v. Assumes that both are unit normalized. See also ⟂.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Sphere.⟂-Tuple{Any,Any}",
    "page": "Private",
    "title": "CMB.Sphere.⟂",
    "category": "Method",
    "text": "⟂(u, v)\n\nTest whether vector u is perpendicular to vector v. Assumes that both are unit normalized. See also ∥.\n\n\n\n"
},

{
    "location": "lib/private.html#Sphere-Functions-1",
    "page": "Private",
    "title": "Sphere Functions",
    "category": "section",
    "text": "Modules = [CMB.Sphere]\nPublic = false"
},

{
    "location": "lib/private.html#CMB.Legendre.Plm_00",
    "page": "Private",
    "title": "CMB.Legendre.Plm_00",
    "category": "Function",
    "text": "Plm_00(::N, ::Type{T}) where {N<:AbstractLegendreNorm, T<:Real}\n\nReturns the initial condition P_0^0(x) for the associated Legendre recursions based on the normalization choice N for numeric type T.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Legendre.Plm_α",
    "page": "Private",
    "title": "CMB.Legendre.Plm_α",
    "category": "Function",
    "text": "Plm_α(norm::N, ::Type{T}, l::Integer, m::Integer) where {N<:AbstractLegendreNorm, T<:Real}\n\nReturns the coefficient _^m for the two-term recursion relation\n\n    P_+1^m(x) = _+1^m x P_^m(x) - _+1^m P_-1^m(x)\n\nwhere _^m is appropriate for the choice of normalization N.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Legendre.Plm_β",
    "page": "Private",
    "title": "CMB.Legendre.Plm_β",
    "category": "Function",
    "text": "Plm_β(norm::N, ::Type{T}, l::Integer, m::Integer) where {N<:AbstractLegendreNorm, T<:Real}\n\nReturns the coefficient _^m for the two-term recursion relation\n\n    P_+1^m(x) = _+1^m x P_^m(x) - _+1^m P_-1^m(x)\n\nwhere _^m is appropriate for the choice of normalization N.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Legendre.Plm_μ",
    "page": "Private",
    "title": "CMB.Legendre.Plm_μ",
    "category": "Function",
    "text": "Plm_μ(norm::N, ::Type{T}, l::Integer) where {N<:AbstractLegendreNorm, T<:Real}\n\nReturns the coefficient _ for the single-term recursion relation\n\n    P_+1^+1(x) = -_+1 sqrt1-x^2 P_^(x)\n\nwhere _ is appropriate for the choice of normalization N.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Legendre.Plm_ν",
    "page": "Private",
    "title": "CMB.Legendre.Plm_ν",
    "category": "Function",
    "text": "Plm_ν(norm::N, ::Type{T}, l::Integer) where {N<:AbstractLegendreNorm, T<:Real}\n\nReturns the coefficient _ for the single-term recursion relation\n\n    P_+1^(x) = _+1 x P_^(x)\n\nwhere _ is appropriate for the choice of normalization N.\n\n\n\n"
},

{
    "location": "lib/private.html#Legendre-Functions-1",
    "page": "Private",
    "title": "Legendre Functions",
    "category": "section",
    "text": "Modules = [CMB.Legendre]\nPublic = false"
},

{
    "location": "lib/private.html#CMB.Healpix.MAX_NSIDE",
    "page": "Private",
    "title": "CMB.Healpix.MAX_NSIDE",
    "category": "Constant",
    "text": "const MAX_NSIDE = 2^29\n\nMaximum valid N_mathrmside parameter value.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Healpix.unsafe_pix2ang-Tuple{Any,Any}",
    "page": "Private",
    "title": "CMB.Healpix.unsafe_pix2ang",
    "category": "Method",
    "text": "(θ,ϕ) = unsafe_pix2ang(nside, p)\n\nLike pix2ang but does not call checkhealpix to check nside and pixel index validity.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Healpix.unsafe_pix2phi-Union{Tuple{I,I}, Tuple{I}} where I<:Integer",
    "page": "Private",
    "title": "CMB.Healpix.unsafe_pix2phi",
    "category": "Method",
    "text": "(θ,ϕ) = unsafe_pix2phi(nside, p)\n\nLike pix2phi but does not call checkhealpix to check nside and pixel index validity.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Healpix.unsafe_pix2theta-Tuple{Any,Any}",
    "page": "Private",
    "title": "CMB.Healpix.unsafe_pix2theta",
    "category": "Method",
    "text": "(θ,ϕ) = unsafe_pix2theta(nside, p)\n\nLike pix2theta but does not call checkhealpix to check nside and pixel index validity.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Healpix.unsafe_pix2vec-Union{Tuple{I,I}, Tuple{I}} where I<:Integer",
    "page": "Private",
    "title": "CMB.Healpix.unsafe_pix2vec",
    "category": "Method",
    "text": "(θ,ϕ) = unsafe_pix2vec(nside, p)\n\nLike pix2vec but does not call checkhealpix to check nside and pixel index validity.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Healpix.unsafe_pix2z-Union{Tuple{I,I}, Tuple{I}} where I<:Integer",
    "page": "Private",
    "title": "CMB.Healpix.unsafe_pix2z",
    "category": "Method",
    "text": "(θ,ϕ) = unsafe_pix2z(nside, p)\n\nLike pix2z but does not call checkhealpix to check nside and pixel index validity.\n\n\n\n"
},

{
    "location": "lib/private.html#Healpix-1",
    "page": "Private",
    "title": "Healpix",
    "category": "section",
    "text": "Modules = [CMB.Healpix]\nPublic = false"
},

{
    "location": "lib/private.html#CMB.PixelCovariance.FIELDMAP",
    "page": "Private",
    "title": "CMB.PixelCovariance.FIELDMAP",
    "category": "Constant",
    "text": "const FIELDMAP\n\nA symbol array which states the canonical ordering of the block matrices within a pixel-pixel covariance matrix.\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.PixelCovariance.SPECTRAMAP",
    "page": "Private",
    "title": "CMB.PixelCovariance.SPECTRAMAP",
    "category": "Constant",
    "text": "const SPECTRAMAP\n\nA symbol array which states the canonical ordering of spectra.\n\n\n\n"
},

{
    "location": "lib/private.html#Pixel-Covariance-1",
    "page": "Private",
    "title": "Pixel Covariance",
    "category": "section",
    "text": "Modules = [CMB.PixelCovariance]\nPublic = false"
},

{
    "location": "lib/private.html#CMB.Util.@abserr-Tuple{Any}",
    "page": "Private",
    "title": "CMB.Util.@abserr",
    "category": "Macro",
    "text": "@abserr fncall(args...)\n\nTakes the function call fncall(args...) and rewrites the expression to calculate the absolute deviation between the regular call and one with all arguments promoted to BigFloat or BigInt.\n\nExample\n\njulia> CMB.Util.@abserr sin(π-1e-5)\n4.1944097844272447e-22\n\nSee Also\n\n@relerr, @absrelerr\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Util.@absrelerr-Tuple{Any}",
    "page": "Private",
    "title": "CMB.Util.@absrelerr",
    "category": "Macro",
    "text": "@absrelerr fncall(args...)\n\nTakes the function call fncall(args...) and rewrites the expression to calculate the absolute and relative deviation between the regular call and one with all arguments promoted to BigFloat or BigInt, returning both as a tuple pair.\n\nExample\n\njulia> CMB.Util.@absrelerr sin(π-1e-5)\n(4.1944097844272447e-22, 0.24759425226749643)\n\nSee Also\n\n@abserr, @relerr\n\n\n\n"
},

{
    "location": "lib/private.html#CMB.Util.@relerr-Tuple{Any}",
    "page": "Private",
    "title": "CMB.Util.@relerr",
    "category": "Macro",
    "text": "@relerr fncall(args...)\n\nTakes the function call fncall(args...) and rewrites the expression to calculate the relative deviation (in ulps) between the regular call and one with all arguments promoted to BigFloat or BigInt.\n\nExample\n\njulia> CMB.Util.@relerr sin(π-1e-5)\n0.24759425226749643\n\nSee Also\n\n@abserr, @absrelerr\n\n\n\n"
},

{
    "location": "lib/private.html#Miscellaneous-Utilities-1",
    "page": "Private",
    "title": "Miscellaneous Utilities",
    "category": "section",
    "text": "Modules = [CMB.Util]\nPublic = false"
},

]}