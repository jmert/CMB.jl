# References

### [BICEP/Keck](@id bib-bicepkeck)
* J. Tolan. “Testing Inflationary Cosmology with BICEP2 and the Keck Array”.
  PhD Thesis. Stanford University, Aug. 2014.
  URL: [http://purl.stanford.edu/sr837cd9569](http://purl.stanford.edu/sr837cd9569)

* The BICEP/Keck Array Collaboration. “BICEP/Keck Array VII: Matrix Based E/B Separation
  Applied to BICEP2 and the Keck Array”
  In: *The Astrophysical Journal* 825, 66 (Jul 2016)
  DOI: [10.3847/0004-637X/825/1/66](http://dx.doi.org/10.3847/0004-637X/825/1/66)
  arXiv: [1603.05976](https://arxiv.org/abs/1603.05976)

* J. Willmert. “Constraining Inflationary B-modes with the BICEP/Keck Array
  Telescopes”.
  PhD Thesis. University of Minnesota, Nov 2019.
  URL: [http://hdl.handle.net/11299/211821](http://hdl.handle.net/11299/211821)

### [Coordinate systems and conventions](@id bib-coordinates)

* J. P. Hamaker and J. D. Bregman. “Understanding radio polarimetry. III.
  Interpreting the IAU/IEEE definitions of the Stokes parameters.”
  In: *Astronomy and Astrophysics, Supplement Series* 117 (May 1996) pp. 161–165.
  DOI: [10.1051/aas:1996147](https://doi.org/10.1051/aas:1996147)

### [HEALPix](@id bib-healpix)

* K. M. Górski et al. “HEALPix: A Framework for High-Resolution Discretization and Fast
  Analysis of Data Distributed on the Sphere”
  In: *The Astrophysical Journal* 662, 2 (Apr 2005) p. 759–771
  DOI: [10.1086/427976](http://dx.doi.org/10.1086/427976).
  arXiv: [astro-ph/0409513](https://arxiv.org/abs/astro-ph/0409513)

  - Erratum: The equatorial-belt longitude ``\phi`` equation (Eq 9) should have the
    ``+1`` outside the modulus:
    ```math
      s = (i - N_\mathrm{side}) \operatorname{mod} 2 + 1
    ```

  - Erratum: Equation 22 differs in signs from equations A2–A3, and neither could be made
    to work. Instead, the implementation here is derived from the system of equations:
    ```math
      \begin{align*}
        z_p(\phi, k_p) &= \frac{2}{3} - \frac{2k_p}{3N_{\mathrm{side}}}
            + \frac{8\phi}{3\pi}
        \\
        z_m(\phi, k_m) &= -\frac{2}{3} + \frac{2k_m}{3N_{\mathrm{side}}}
            - \frac{8\phi}{3\pi}
      \end{align*}
    ```

### [Pixel Covariance](@id bib-pixelcovariance)

* M. Tegmark and A. de Oliveira-Costa. “How to measure CMB polarization power spectra
  without losing information”
  In: *Physical Review D* 64, 063001 (Sep 2001)
  DOI: [10.1103/PhysRevD.64.063001](http://dx.doi.org/10.1103/PhysRevD.64.063001)
  arXiv: [astro-ph/0012120](https://arxiv.org/abs/astro-ph/0012120)

