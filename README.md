# CMB.jl — CMB Analysis

<!--
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jmert.github.io/CMB.jl/stable)
-->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jmert.github.io/CMB.jl/dev)
[![Build Status](https://travis-ci.com/jmert/CMB.jl.svg?branch=master)](https://travis-ci.com/jmert/CMB.jl)
[![Codecov](https://codecov.io/gh/jmert/CMB.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jmert/CMB.jl)

`CMB.jl` is a library of routines for the analysis of cosmic microwave
background (CMB) data. Development of features is being driven by the author's
use cases — at this time, namely the production of “reobserved” pixel-pixel
covariance matrices as used by the BICEP/Keck Array collaboration.

Design goals of this package include:

  * Native Julia implementation of core routines.

  * Numerical stability and efficiency.

  * Parallelism and efficient memory sharing.


