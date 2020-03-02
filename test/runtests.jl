
using Test, Documenter, Nemo, GroupAtlas

DocMeta.setdocmeta!( GroupAtlas, :DocTestSetup, :(begin using GroupAtlas; using Nemo; end);
                     recursive = true )

include( "testmeataxe.jl" )
include( "testaccess.jl" )
include( "testmanual.jl" )

