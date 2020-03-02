# GroupAtlas.jl

```@contents
```

## Introduction

GroupAtlas.jl is a Julia interface to a database of generators
of permutation and matrix groups,
including the ones that are available in the
[Atlas of Group Representations](http://brauer.maths.qmul.ac.uk/Atlas/v3/).

More precisely, the interface gives access to the data that are listed for
[the GAP interface](http://www.math.rwth-aachen.de/~Thomas.Breuer/atlasrep/index.html).
Currently these are

- [the original Atlas data](http://www.math.rwth-aachen.de/~Thomas.Breuer/atlasrep/atlasprm.json)

- [additional data available in the GAP interface](http://www.math.rwth-aachen.de/~Thomas.Breuer/atlasrep/datapkg/toc.json)

- [additional multiplicity-free permutation groups from GAP's MFER package](http://www.math.rwth-aachen.de/~mfer/mfertoc.json)

- [additional groups from GAP's CTBlocks package](http://www.math.rwth-aachen.de/~Thomas.Breuer/ctblocks/ctblockstoc.json)

The Atlas of Group Representations provides also shell scripts that can be used
to restrict representations to certain (maximal) subgroups.
The current version of GroupAtlas.jl does not provide this information.

The permutations and matrices that are returned by the interface functions
are currently in formats that are provided by Julia's
[Nemo](http://www.nemocas.org) package.

Here are a few examples how to use GroupAtlas.jl.

```julia
julia> DisplayAtlasInfo()       # show an overview of groups
group                    │  #
─────────────────────────┼───
2.2E6(2)                 │  1
2.2E6(2).2               │  2
2.(2xF4(2)).2            │  1
2.A5                     │ 25
2.A5.2                   │ 11
2.A6                     │ 17
[...]
(S5xS5xS5):S3            │  1
A5                       │ 18
A5.2                     │  9
[...]
U6(3)                    │  2
U7(2)                    │  4
U8(2)                    │  2
W(F4)                    │  0

julia> DisplayAtlasInfo( "A5" )   # show an overview for one group
Representations for A5:
───────────────────────
 1: G ≤ Sym(5)                  3-trans., on cosets of A4 (1st max.)
 2: G ≤ Sym(6)                  2-trans., on cosets of D10 (2nd max.)
 3: G ≤ Sym(10)                 rank 3, on cosets of S3 (3rd max.)
 4: G ≤ GL(4a,2)                φ = 4a
 5: G ≤ GL(4b,2)                φ = 2ab
 6: G ≤ GL(4,3)                 φ = 4a
 7: G ≤ GL(6,3)                 φ = 3ab
 8: G ≤ GL(2a,4)                φ = 2a
 9: G ≤ GL(2b,4)                φ = 2b
10: G ≤ GL(3,5)                 φ = 3a
11: G ≤ GL(5,5)                 φ = 5a
12: G ≤ GL(3a,9)                φ = 3a
13: G ≤ GL(3b,9)                φ = 3b
14: G ≤ GL(4,ℤ)                 χ = 4a
15: G ≤ GL(5,ℤ)                 χ = 5a
16: G ≤ GL(6,ℤ)                 χ = 3ab
17: G ≤ GL(3a,Field([Sqrt(5)])) χ = 3a
18: G ≤ GL(3b,Field([Sqrt(5)])) χ = 3b

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

julia> info = OneAtlasInfo( Dict( :name => "A5",   # choose a generating set
                                  :char => 2,      # in characteristic 2
                                  :dim => 4 ) )    # and with dimension 4
Dict{Symbol,Any} with 7 entries:
  :dim     => 4
  :id      => "a"
  :path    => "alt/A5/mtx/A5G1-f2r4aB0.m"
  :name    => "A5"
  :q       => 2
  :tocID   => "core"
  :filenum => 2

julia> AtlasGenerators( info )   # fetch the generators
2-element Array{Nemo.fq_nmod_mat,1}:
 [1  0  0  0]
[0  0  1  0]
[0  1  0  0]
[1  1  1  1]
 [0  1  0  0]
[0  0  0  1]
[0  0  1  0]
[1  0  0  0]

julia> F = Nemo.GaloisField( UInt(2) )
Galois field with characteristic 2

julia> info = OneAtlasInfo( Dict( :name => "A5",   # choose a generating set
                                  :field => F,     # over a prescribed field
                                  :dim => 4 ) )    # and with dimension 4
Dict{Symbol,Any} with 8 entries:
  :dim     => 4
  :field   => Galois field with characteristic 2
  :id      => "a"
  :path    => "alt/A5/mtx/A5G1-f2r4aB0.m"
  :name    => "A5"
  :q       => 2
  :tocID   => "core"
  :filenum => 2

julia> AtlasGenerators( ans )
2-element Array{gfp_mat,1}:
 [1  0  0  0]
[0  0  1  0]
[0  1  0  0]
[1  1  1  1]
 [0  1  0  0]
[0  0  0  1]
[0  0  1  0]
[1  0  0  0]

```

This package does not depend on 
[the Julia package GAP.jl](https://github.com/oscar-system/GAP.jl),
but if GAP.jl is available then the function [`BrowseAtlasInfo`(@ref) can be used
to choose database entries from an interactive table.
This behaviour is implemented via
[the Julia package Requires.jl](https://github.com/JuliaPackaging/Requires.jl),
and the function is based on the function `BrowseAtlasInfo` from
[the GAP package Browse](http://www.math.rwth-aachen.de/~Browse).

The data files are accessed as follows.

- If [the Julia package GAP.jl](https://github.com/oscar-system/GAP.jl)
  has been loaded before GroupAtlas.jl then it is checked
  whether the requested data file is already available locally,
  in the directory where
  [the GAP package AtlasRep](http://www.math.rwth-aachen.de/~Thomas.Breuer/atlasrep)
  would store the file after download;
  if yes then the contents of this local file is taken.

- If GAP.jl has not been loaded or if the underlying GAP installation does not
  provide a local version of the requested data file then
  [the Julia package HTTP](https://github.com/JuliaWeb/HTTP.jl)
  is used to fetch the file from its web address,
  as stated in one of the tables of contents listed above.
  In this case, each request for the same data file will fetch the file anew,
  no caching is provided.

```@meta
CurrentModule = GroupAtlas
```

## User Functions

```@docs
DisplayAtlasInfo
OneAtlasInfo
AtlasGenerators
BrowseAtlasInfo
```

## Auxiliary Functions

### Read and write MeatAxe text files

```@docs
TuplesIterator
finitefieldlists
finitefieldlist
meataxe_file_header_info
nemo_permutation_from_meataxe_file
nemo_matrix_from_meataxe_file
meataxe_string
```

### Read JSON format files

```@docs
nemo_matrices_from_JSON_file
```

### Get and evaluate the contents of local or remote files

```@docs
contents_textfile
scanfilename
group_atlas_filenames
```

### Administrate of tables of contents

```@docs
GroupAtlasData
AtlasTOC
unite_dictionaries!
isless_groupnames
isless_entries_for_group
notify_data_file
```

## Index

```@index
```

