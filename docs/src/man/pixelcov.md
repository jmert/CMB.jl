# [Pixel-Pixel Covariance](@id man_pixelcov)
```@meta
DocTestFilters = Regex[
        r"Ptr{0x[0-9a-f]+}",
        r"[0-9\.]+ seconds( \(.*\))?",
        ]
```

```@contents
Pages = ["pixelcov.md"]
Depth = 2
```

## [Definition](@id pixelcov_defn)

The following definitions are a lightly-modified reproduction of Appendix A
of [Tegmark & Oliveira-Costa (2001)](@ref bib-pixelcovariance).

```@raw latex
\providecommand{\mat}[1]{\symbf{#1}}
\providecommand{\expv}[1]{\left\langle #1\right\rangle}
\providecommand{\covF}[1]{F_\ell^{#1}(z_{ij})}
\providecommand{\covTT}{\expv{T_i T_j}}
\providecommand{\covTQ}{\expv{T_i Q_j}}
\providecommand{\covTU}{\expv{T_i U_j}}
\providecommand{\covQQ}{\expv{Q_i Q_j}}
\providecommand{\covQU}{\expv{Q_i U_j}}
\providecommand{\covUU}{\expv{U_i U_j}}
\providecommand{\cij}{c_{ij}}
\providecommand{\sij}{s_{ij}}
\providecommand{\cji}{c_{ji}}
\providecommand{\sji}{s_{ji}}
```
```@raw html
<div style="display:none;">
\(
\newcommand{\mat}[1]{\mathbf{#1}}
\newcommand{\expv}[1]{\left\langle #1\right\rangle}
\newcommand{\covF}[1]{F_\ell^{#1}(z_{ij})}
\newcommand{\covTT}{\expv{T_i T_j}}
\newcommand{\covTQ}{\expv{T_i Q_j}}
\newcommand{\covTU}{\expv{T_i U_j}}
\newcommand{\covQQ}{\expv{Q_i Q_j}}
\newcommand{\covQU}{\expv{Q_i U_j}}
\newcommand{\covUU}{\expv{U_i U_j}}
\newcommand{\cij}{c_{ij}}
\newcommand{\sij}{s_{ij}}
\newcommand{\cji}{c_{ji}}
\newcommand{\sji}{s_{ji}}
\)
</div>
```

The pixel-pixel covariance is best constructed in stages.
We start by defining the covariance ``\langle X_i Y_j\rangle`` between any
two Stokes fields ``X`` and ``Y`` at pixels ``i`` and ``j`` for a local coordinate
system where the ``Q`` axis is pointed along the great circle connection
pixel ``i`` to pixel ``j``.
In this simplified case, the 6 unique covariances are
```math
\begin{align}
    \covTT &\equiv \frac{1}{4\pi} \sum_\ell C_\ell^{TT} \covF{00}
        \label{eqn:theorycov:tt}
    \\
    \covTQ &\equiv -\frac{1}{4\pi} \sum_\ell C_\ell^{TE} \covF{10}
        \label{eqn:theorycov:tq}
    \\
    \covTU &\equiv -\frac{1}{4\pi} \sum_\ell C_\ell^{TB} \covF{10}
        \label{eqn:theorycov:tu}
    \\
    \covQQ &\equiv \frac{1}{4\pi} \sum_\ell
        \left[ C_\ell^{EE} \covF{12} - C_\ell^{BB} \covF{22} \right]
        \label{eqn:theorycov:qq}
    \\
    \covUU &\equiv \frac{1}{4\pi} \sum_\ell
        \left[ C_\ell^{BB} \covF{12} - C_\ell^{EE} \covF{22} \right]
        \label{eqn:theorycov:uu}
    \\
    \covQU &\equiv \frac{1}{4\pi} \sum_\ell
        C_\ell^{EB} \left[ \covF{12} + \covF{22} \right]
        \label{eqn:theorycov:qu}
\end{align}
```
for ``z_{ij} = \cos(\sigma_{ij})`` and some fiducial spectrum ``C_\ell``.
The polarization weighting functions are simple functions of the ``P_\ell`` and
``P_\ell^2``
[associated Legendre polynomials](https://jmert.github.io/AssociatedLegendrePolynomials.jl/stable/),
```math
\begin{align}
    \covF{00} &\equiv (2\ell + 1) P_\ell(z_{ij})
        \label{eqn:theorycov:F00}
    \\
    \covF{10} &\equiv \chi_\ell
        \left[
        \frac{z_{ij}}{1-z_{ij}^2} P_{\ell-1}(z_{ij}) - \left(
        \frac{1}{1-z_{ij}^2} + \frac{\ell-1}{2} \right)P_\ell(z_{ij})
        \right]
        \label{eqn:theorycov:F10}
    \\
    \covF{12} &\equiv \gamma_\ell
        \left[
        \frac{\ell+2}{1-z_{ij}^2} z_{ij} P^2_{\ell-1}(z_{ij}) - \left(
        \frac{\ell-4}{1-z_{ij}^2} + \frac{\ell(\ell-1)}{2} \right)
        P^2_\ell(z_{ij})
        \right]
        \label{eqn:theorycov:F12}
    \\
    \covF{22} &\equiv 2\gamma_\ell
        \left[
        \frac{\ell+2}{1-z_{ij}^2} P^2_{\ell-1}(z_{ij}) -
        \frac{\ell-1}{1-z_{ij}^2} z_{ij} P^2_\ell(z_{ij})
        \right]
        \label{eqn:theorycov:F22}
\end{align}
```
where[^1]
```math
\begin{align}
    \chi_\ell &\equiv \frac{2\ell (2\ell + 1)}
        {\sqrt{(\ell-1)\ell(\ell+1)(\ell+2)}}
    \\
    \gamma_\ell &\equiv \frac{2(2\ell + 1)}
        {(\ell-1)\ell(\ell+1)(\ell+2)}
\end{align}
```


For a single pixel-pixel pair, the covariance terms form a symmetric
``3 \times 3`` matrix ``\mat M``
```math
\begin{align}
    \mat M(z_{ij}) &= \begin{bmatrix}
        \covTT & \covTQ & \covTU \\
        \covTQ & \covQQ & \covQU \\
        \covTU & \covQU & \covUU
    \end{bmatrix}
\end{align}
```
which is rotated from the local to a global Stokes coordinate system by
application of the rotation matrix
```math
\begin{align}
    \mat R(\alpha) &= \begin{bmatrix}
        1 & 0 & 0 \\
        0 & \cos(2\alpha) & -\sin(2\alpha) \\
        0 & \sin(2\alpha) &  \cos(2\alpha)
    \end{bmatrix}
\end{align}
```
where ``\alpha`` are bearing angles, such that the pixel-pixel covariance
``{\mat C}_{ij}`` in global coordinates is defined by
```math
\begin{align}
    {\mat C}_{ij} &= \mat R(\alpha_{ij}) \mat M(z_{ij}) \mat R(\alpha_{ji})^\top
        \label{eqn:theorycov:rotate}
\end{align}
```
Pay attention to the fact that the left-hand ``\mat R`` uses the bearing angle
``\alpha_{ij}`` at pixel ``i`` whereas the right-hand ``\mat R`` uses the bearing
angle ``\alpha_{ji}`` at pixel ``j``.
Also note that these rotations move the covariances into the IAU
polarization convention — to move into the `HEALPix` polarization convention,
the rotations both need to be reversed, i.e. ``\mat R \leftrightarrow \mat
R^\top``.

The pixel-pixel covariance matrix ``\mat C`` for some set of pixels is constructed by
calculating the terms in ``{\mat C}_{ij}`` for each pair of pixels and filling in the 9
values as the ``(i,j)``-th entries of the covariance blocks in ``\mat C``.
By blocks, we refer to the fact that ``\mat C`` can be block-decomposed as
```math
\begin{align}
    \mat C &= \begin{bmatrix}
        \mat C^{TT} & \mat C^{TQ} & \mat C^{TU} \\
        \mat C^{QT} & \mat C^{QQ} & \mat C^{QU} \\
        \mat C^{UT} & \mat C^{UQ} & \mat C^{UU}
    \end{bmatrix}
\end{align}
```
For example, if one holds $j$ constant and varies $i \in P$ across all pixels,
then one would form 3 columns of $\mat C$, one in each of the block-columns.

It is also convenient to define the polarization-only block decompositions of the total
covariance matrix as
```math
\begin{align}
    \mat C^{\mathrm{Pol}} &\equiv \begin{bmatrix}
        \mat C^{QQ} & \mat C^{QU} \\
        \mat C^{UQ} & \mat C^{UU} \end{bmatrix}
\end{align}
```

## [Ideal Pixel-Pixel Covariance](@id pixelcov_theory)

## [Reobserved Pixel-Pixel Covariance](@id pixelcov_reobs)

## [Properties and Symmetries](@id pixelcov_properties)

- The entire covariance matrix ``\mat C`` is symmetric.
  Therefore the on-diagonal sub-blocks ``\mat C^{TT}``, ``\mat C^{QQ}``, and
  ``\mat C^{UU}`` are individually symmetric as well, and the off-diagonal blocks must
  be related to one another as
  ```math
  \begin{align*}
      \mat C^{TQ} &= \left(\mat C^{QT}\right)^\top &
      \mat C^{TU} &= \left(\mat C^{UT}\right)^\top &
      \mat C^{QU} &= \left(\mat C^{UQ}\right)^\top
  \end{align*}
  ```
  This directly implies the polarization-only matrix ``\mat C^{\mathrm{Pol}}`` is
  symmetric as well.
- The covariance matrices ``\mat C``, ``\mat C^{TT}``, and ``\mat C^{\mathrm{Pol}}`` are at
  least positive-semidefinite for appropriate non-zero fiducial input spectra.
  Positive-definiteness only occurs when there are at least as many non-zero harmonic modes
  (``C_\ell^{XY} \neq 0`` over all spectra ``XY``) as there are diagonal elements in the
  matrix (the number of pixels in the map(s) that the covariance describes).
- The covariance matrices are linear in the fiducial spectra.
  For example, an ``EE``-only covariance matrix ``[\mat C]^{EE}`` (wherein all ``C_\ell =
  0`` except for ``C_\ell^{EE}`` which is non-zero somewhere) and an ``EB``-only covariance
  matrix ``[\mat C]^{EB}`` can be summed to define
  ```math
  \begin{align*}
      [\mat C]^{EEEB} = [\mat C]^{EE} + [\mat C]^{EB}
  \end{align*}
  ```
  which is equivalent to the covariance matrix which would have been produced if the
  fiducial spectrum had included the both of the ``C_\ell^{EE}`` and ``C_\ell^{EB}``
  spectra from the start.
- Since the IAU and HEALPix polarization conventions differ by the direction of rotation
  in ``\mat R(\alpha)`` which corresponds to a change in the sign of the ``\sin`` terms
  (``\sij`` and ``\sji``), the cosmologically-interesting case where ``C_\ell^{EB} = 0``
  and ``C_\ell^{TB} = 0`` also simplifies such that the following are also true:
  ```math
  \begin{align*}
      {\mat C}_{ij,\mathrm{Healpix}}^{TU} &= -{\mat C}_{ij,\mathrm{IAU}}^{TU} &
      {\mat C}_{ij,\mathrm{Healpix}}^{UT} &= -{\mat C}_{ij,\mathrm{IAU}}^{UT} \\
      {\mat C}_{ij,\mathrm{Healpix}}^{QU} &= -{\mat C}_{ij,\mathrm{IAU}}^{QU} &
      {\mat C}_{ij,\mathrm{Healpix}}^{UQ} &= -{\mat C}_{ij,\mathrm{IAU}}^{UQ}
  \end{align*}
  ```
  and the remaining block components are unchanged.

## [Mathematical Details](@id pixelcov_details)

It is often useful to have the fully-expanded expressions for each of the pixel-pixel
covariance terms after applying the local-to-global coordinate system rotations.
If we define short-cut notation for each of the terms in the rotation matrices
as
```math
\begin{align*}
    \mat R(\alpha_{ij}) &\equiv \begin{bmatrix}
        1 & 0 & 0 \\
        0 & \cij & -\sij \\
        0 & \sij &  \cij
    \end{bmatrix}
    &
    \mat R(\alpha_{ji})^\top &\equiv \begin{bmatrix}
        1 & 0 & 0 \\
        0 &  \cji & \sji \\
        0 & -\sji & \cji
    \end{bmatrix}
\end{align*}
```
then expanding Eqn. ``\ref{eqn:theorycov:rotate}`` explicitly (and grouping the terms
by block-columns):
```math
\begin{align*}
    \begin{bmatrix}
            {\mat C}_{ij}^{TT} \\ {\mat C}_{ij}^{QT} \\ {\mat C}_{ij}^{UT}
        \end{bmatrix} &= \begin{bmatrix}
            \covTT \\
            \covTQ\cij - \covTU\sij \\
            \covTQ\sij + \covTU\cij
        \end{bmatrix}
    \\
    \begin{bmatrix}
            {\mat C}_{ij}^{TQ} \\ {\mat C}_{ij}^{QQ} \\ {\mat C}_{ij}^{UQ}
        \end{bmatrix} &= \begin{bmatrix}
            \covTQ\cji - \covTU\sji \\
            \covQQ\cij\cji - \covQU(\cij\sji + \sij\cji) + \covUU\sij\sji \\
            \covQQ\sij\cji + \covQU(\cij\cji - \sij\sji) - \covUU\cij\sji
        \end{bmatrix}
    \\
    \begin{bmatrix}
            {\mat C}_{ij}^{TU} \\ {\mat C}_{ij}^{QU} \\ {\mat C}_{ij}^{UU}
        \end{bmatrix} &= \begin{bmatrix}
            \covTQ\sji + \covTU\cji \\
            \covQQ\cij\sji + \covQU(\cij\cji - \sij\sji) - \covUU\sij\cji \\
            \covQQ\sij\sji + \covQU(\cij\sji + \sij\cji) + \covUU\cij\cji
        \end{bmatrix}
\end{align*}
```

---

### Footnotes

[^1]:
    In the limit as ``|z| \rightarrow 1`` (i.e. as the two points given by
    pointing vectors are [anti-]parallel), the ``\left(1 - z^2\right)^{-1}``
    terms diverge. The resolution is to special-case the numerical computation
    and make use of the mathematical limits,
    ```math
    \begin{align}
        \covF{10} &= {}\rlap{\,\,0}
            \hphantom{\left\{ (-1)^\ell \frac{2\ell+1}{2} \right.}
            \,\text{as } |z_{ij}| \to 1
        \\
        \covF{12} &= \begin{cases}
                \frac{2\ell+1}{2}            & \text{as } z_{ij} \to +1 \\
                \frac{2\ell+1}{2}  (-1)^\ell & \text{as } z_{ij} \to -1
            \end{cases}
        \\
        \covF{22} &= \begin{cases}
                -\frac{2\ell+1}{2}           & \text{as } z_{ij} \to +1 \\
                 \frac{2\ell+1}{2} (-1)^\ell & \text{as } z_{ij} \to -1
            \end{cases}
    \end{align}
    ```
    Furthermore, the bearing angle is not well defined for [anti-]parallel
    points.
    Thankfully, ``\langle Q_i Q_j\rangle = \langle U_i U_j \rangle`` and
    ``\langle Q_i U_j \rangle = 0``, so the polarization weights are invariant
    under rotation;
    therefore we can take ``\alpha = 0`` without loss of generality.
    (See [Tegmark & Oliveira-Costa (2001)](@ref bib-pixelcovariance).)
