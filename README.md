[![Build Status](https://travis-ci.com/oscar-system/GroupAtlas.jl.svg?branch=master)](https://travis-ci.com/oscar-system/GroupAtlas.jl)
[![Code Coverage](https://codecov.io/github/oscar-system/GroupAtlas.jl/coverage.svg?branch=master&token=)](https://codecov.io/gh/oscar-system/GroupAtlas.jl)

# GroupAtlas.jl Julia module

## Basic usage

After entering the following in Julia,
```julia
julia> using GroupAtlas
```
one can access the permutation and matrix representations
from the [ATLAS of Group Representations](http://brauer.maths.qmul.ac.uk/Atlas/v3/),
using the permutations and matrices provided by [Nemo](http://www.nemocas.org/) to represent them in Julia.
```julia
julia> info = OneAtlasInfo( Dict( :name => "A5" ) )   # choose a generating set
Dict{Symbol,Any} with 6 entries:
  :id      => "" 
  :degree  => 5
  :name    => "A5"
  :path    => "alt/A5/mtx/A5G1-p5B0.m"
  :tocID   => "core"
  :filenum => 2

julia> AtlasGenerators( info )   # fetch the generators
2-element Array{AbstractAlgebra.Generic.Perm{Int64},1}:
 (1,2)(3,4)
 (1,3,5)  

```

This software is licensed under the LGPL, version 3, or any later version.

