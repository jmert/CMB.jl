# [Pixel-pixel Covariance](@id man_pixelcov)
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

## [Definition and Properties](@id pixelcov_defn)

The following definitions are a lightly-modified reproduction of Appendix A
of [Tegmark & Oliveira-Costa (2001)](@ref bib-pixelcovariance).

```@raw html
<div style="display:none;">
```
```math
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
```
```@raw html
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
``P_\ell^2`` associated Legendre polynomials,
```math
\begin{align}
    \covF{00} &= (2\ell + 1) P_\ell(z_{ij})
        \label{eqn:theorycov:F00}
    \\
    \covF{10} &= \chi_\ell
        \left[
        \frac{z_{ij}}{1-z_{ij}^2} P_{\ell-1}(z_{ij}) - \left(
        \frac{1}{1-z_{ij}^2} + \frac{\ell-1}{2} \right)P_\ell(z_{ij})
        \right]
        \label{eqn:theorycov:F10}
    \\
    \covF{12} &= \gamma_\ell
        \left[
        \frac{\ell+2}{1-z_{ij}^2} z_{ij} P^2_{\ell-1}(z_{ij}) - \left(
        \frac{\ell-4}{1-z_{ij}^2} + \frac{\ell(\ell-1)}{2} \right)
        P^2_\ell(z_{ij})
        \right]
        \label{eqn:theorycov:F12}
    \\
    \covF{22} &= 2\gamma_\ell
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
