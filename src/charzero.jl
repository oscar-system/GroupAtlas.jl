
"""
    nemo_matrices_from_JSON_file( filename::String, z = nothing )

Return the array of Nemo matrices over the ring of integers or a number field
that are described by the JSON format file given by `filename`.

If `filename` starts with `http://` then `filename` is treated as an URL,
and the function attempts to download the file.
Otherwise `filename` is assumed to be the name of a (perhaps gzipped)
local file.

If `z` is given, it must be a root of the rational polynomial given by
`filename`,
and the entries of the returned matrices belong to the same number field as `z`.
If `z` is not given then the number field is constructed anew.

# Example
```jldoctest
julia> testdir = joinpath( dirname( pathof( GroupAtlas ) ), "../test" );

julia> GroupAtlas.nemo_matrices_from_JSON_file(
           joinpath( testdir, "2A5G1-Ar2aB0.json" ) )
2-element Array{AbstractAlgebra.Generic.MatSpaceElem{nf_elem},1}:
 [0  -1]
[1   0]                                      
 [                   -1  z]
[z^3 + 1*z^2 + 1*z + 1  0]

julia> F, z = Nemo.CyclotomicField( 5, "a" );

julia> mats = GroupAtlas.nemo_matrices_from_JSON_file(
                  joinpath( testdir, "2A5G1-Ar2aB0.json" ), z )
2-element Array{AbstractAlgebra.Generic.MatSpaceElem{nf_elem},1}:
 [0  -1]
[1   0]                                      
 [                   -1  a]
[a^3 + 1*a^2 + 1*a + 1  0]

julia> parent( mats[1][1,1] ) == F
true

```
"""
function nemo_matrices_from_JSON_file(filename::String, z = nothing)
    # Try to read the JSON format file.
    str = contents_textfile(filename)

    prs = JSON.parse(str; dicttype = Dict{Symbol,Any})
    ringinfo = prs[:ringinfo]
    polynomial = convert(Array{Int,1}, prs[:polynomial])
    dims = prs[:dimensions]
    generators = prs[:generators]

    if ringinfo[1] == "CyclotomicField"
        if isnothing(z)
            # We need only the primitive root of unity.
            F, z = Nemo.CyclotomicField(ringinfo[2], "z")
        else
            # z must be an appropriate root of unity in F
            F = parent(z)
#T how could one find an appropriate root of unity z in a given field F:
#T R, x = Nemo.PolynomialRing( F, "x" )
#T pol = R( polynomial )
#T but how to take a root of pol in F?
        end
    elseif ringinfo[1] == "IntegerRing"
        if isnothing(z)
            F = Nemo.ZZ
            z = one(F)
        else
            F = parent(z)
        end
    elseif ringinfo[1] == "QuadraticField" || ringinfo[1] == "AbelianNumberField"
        if isnothing(z)
            R, x = Nemo.PolynomialRing(Nemo.QQ, "x")
            F, z = Nemo.NumberField(R(polynomial), "z")
        else
            F = parent(z)
        end
    else
        # We do not have such representations yet.
        error("unknown 'ringinfo' ", ringinfo)
    end

    N = length(polynomial) - 1
    bas = [z^i for i in 0:(N-1)]

    # Decode the non-integral matrix entries.
    for mat in generators
        for i in 1:length(mat)
            if isa(mat[i], Array)
                # coefficient list
                mat[i] = sum(x -> mat[i][x] * bas[x], 1:N)
            end
        end
    end

    # Convert the matrices.
    denom = prs[:denominator]
    return [Nemo.divexact(Nemo.matrix(F, dims[1], dims[2], m), denom) for m in generators]
end
