
export OneAtlasInfo, AllAtlasInfos

# SomeAtlasInfo( conditions::Dict, toc::AtlasTOC, mode::Symbol )
#
# do the work for the user functions `OneAtlasInfo` and `AllAtlasInfos`.
#
function SomeAtlasInfo(conditions::Dict, toc::AtlasTOC, mode::Symbol)

    result = []

    if haskey(conditions, :name)
        # Restrict the search to entries for the group with this name.
        if haskey(GroupAtlasData.gnan, conditions[:name])
            name = conditions[:name]
            list = Dict(:name => toc.toc[name])
        elseif mode === :all
            return []
        else
            return nothing
        end
    else
        # Run over all entries.
        list = toc.toc
    end

    # Run over all group names.
    for k1 in keys(list)
        # Run over all entries for this group.
        for entry in list[k1]
            ok = true
            for k in keys(conditions)
                if k === :char
                    # This is relevant for matrix representations.
                    if haskey(entry, :degree)
                        ok = false
                        break
                    elseif haskey(entry, :q)
                        # finite characteristic
                        if conditions[k] == 0
                            ok = false
                            break
                        elseif mod(entry[:q], conditions[k]) != 0
                            ok = false
                            break
                        end
                    elseif conditions[k] != 0
                        # characteristic zero
                        ok = false
                        break
                    end
                elseif k === :fieldgen
                    # This is relevant for matrices over number fields.
                    if haskey(entry, :ringinfo) && entry[:ringinfo] == "Integers"
                        # This is o.k., remember the desired field element.
                        entry = copy(entry)
                        entry[:gen] = conditions[k]
                    else
                        # Check whether the given generator of the number field
                        # is a zero of the polynomial.
                        if !haskey(entry, :polcoeffs)
                            ok = false
                            break
                        end
                        z = conditions[k]
                        coeffs = reverse(entry[:polcoeffs])
                        val = coeffs[1]
                        for i in 2:length(coeffs)
                            val = val * z + coeffs[i]
                        end
                        if !iszero(val)
                            ok = false
                            break
                        end
                        # Remember the desired field element.
                        entry = copy(entry)
                        entry[:gen] = z
# One would like to cover also the case
# that the matrices live over the parent field of the given generator,
# that is, this field contains a root of the polynomial of `entry`;
# then one can create the matrices in terms of the given generator.
# But currently a factorization of polynomials over number fields
# seems to be not supported.
# The documentation for this case would look as follows.
# `:field`
#     specifies a number field $F$, say;
#     if the polynomial given by the `:polcoeffs` value of the entry
#     has a root in $F$ then regard the entry as a solution,
#     and later use this root in a call of `AtlasGenerators`
#     to fill the matrices.
                    end
                elseif k === :field
                    # This is relevant for matrices over finite fields.
                    if !haskey(entry, :q)
                        ok = false
                        break
                    end
                    f = fieldorder(conditions[k])
                    qq = entry[:q]
                    while f > qq
                        f = f // qq
                    end
                    if f != qq
                        ok = false
                        break
                    end
                    # Remember the desired field.
                    entry = copy(entry)
                    entry[:field] = conditions[k]
                elseif !(
                    (k === :name) || (
                        haskey(entry, k) &&
                        (conditions[k] == Any || entry[k] == conditions[k])
                    )
                )
                    ok = false
                    break
                end
            end
            if ok
                if mode === :all
                    push!(result, entry)
                else
                    return entry
                end
            end
        end
    end

    if mode === :all
        return result
    else
        return nothing
    end
end

"""
    OneAtlasInfo( conditions::Dict )
    AllAtlasInfos( conditions::Dict )

Return one dictionary or all dictionaries describing Group Atlas data,
that match the conditions in `conditions`.

All those keys that occur in the info dictionaries of the database entries
are supported as keys in `conditions`.
If all values in `conditions` match with the corresponding values of
a database entry then this entry is regarded as a solution to the request.

The following keys are most relevant.

`:name`
    the name (a string) of the desired group,

`:degree`
    the degree (an integer) of permutation representations,

`:q`
    the order of the minimal finite field (a prime power)
    over which a matrix representation can be written,

`:dim`
    the dimension (an integer) of the generating matrices,

`:id`
    an identifier (a string) such as `"a"` or `"b"`,
    distinguishing representations of the same degree or dimension.

Besides these keys, the following special keys are supported.

`:char`
    specifies the characteristic (0 or a prime integer)
    of matrix representations;
    if an entry has a `:q` value that is divisible by the `:char` value
    then it matches,
    if an entry has no `:q` value but a `:dim` value then it belongs to
    characteristic 0.

`:field`
    specifies a finite field `F`, say;
    if the `:q` value of an entry is the order of a subfield of `F`
    then regard this entry as a solution, and later use `F`
    in a call of [`AtlasGenerators`](@ref) to fill the matrices.

`:fieldgen`
    specifies an element `z`, say, of a number field;
    if `z` is a root of the polynomial given by the `:polcoeffs` value
    of an entry then regard this entry as a solution, and later use `z`
    in a call of [`AtlasGenerators`](@ref) to fill the matrices.

# Examples
```jldoctest
julia> OneAtlasInfo( Dict( :name => "A5" ) )
Dict{Symbol,Any} with 6 entries:
  :id      => ""
  :degree  => 5
  :name    => "A5"
  :path    => "alt/A5/mtx/A5G1-p5B0.m"
  :tocID   => "core"
  :filenum => 2

julia> OneAtlasInfo( Dict( :name => "A5", :dim => 3 ) )
Dict{Symbol,Any} with 7 entries:
  :dim     => 3
  :id      => ""
  :path    => "alt/A5/mtx/A5G1-f5r3B0.m"
  :name    => "A5"
  :q       => 5
  :tocID   => "core"
  :filenum => 2

julia> OneAtlasInfo( Dict( :name => "A5", :degree => 4 ) )

julia> infos = AllAtlasInfos( Dict( :name => "2.A5", :dim => 2, :q => 9 ) );

julia> length( infos )
2

julia> AllAtlasInfos( Dict( :name => "A5", :degree => 4 ) )
0-element Array{Any,1}

```
"""
function OneAtlasInfo(conditions::Dict)
    return SomeAtlasInfo(conditions, AtlasTOC(), :one)
end

function AllAtlasInfos(conditions::Dict)
    return SomeAtlasInfo(conditions, AtlasTOC(), :all)
end
