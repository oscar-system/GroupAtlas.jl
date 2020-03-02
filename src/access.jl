
export DisplayAtlasInfo, AtlasGenerators

"""
    group_atlas_filenames( entry::Dict )

Return the array of URLs or (if available, and then preferred)
local paths where the data file(s) for `entry` can be found,
which is assumed to be a return value of [`OneAtlasInfo`](@ref)
or an entry in an array returned by `AllAtlasInfos`.
"""
function group_atlas_filenames(entry::Dict)

    tocurl = AtlasTOCs[2][entry[:tocID]]
    toc = AtlasTOCs[1][tocurl]
    path = entry[:path]

    # Search for files that belong to the database.
    if toc.localDirectory != ""
        # Perhaps the data files are stored locally.
        dir = toc.localDirectory
        if isdir(dir)
            # The local directory is available.
            # Check whether the files are already there.
            if entry[:tocID] == "core"
                # In the case of `core` files, remove the directory part from `path`.
                pos = findlast(isequal('/'), path)
                path = SubString(path, pos + 1)
            end
            if (occursin("-Ar", path) || occursin("-Zr", path)) && endswith(path, ".g")
                # In case of matrices in char. zero, use the JSON format file,
                # which is in the same local directory but with a different name.
                path = SubString(path, 1, length(path) - 2) * ".json"
            end
            path = dir * path
            if entry[:filenum] == 1
                if isfile(path)
                    return [path]
                elseif isfile(path * ".gz")
                    return [path * ".gz"]
                end
            else
                list = map(i -> path * string(i), 1:entry[:filenum])
                filesfound = true
                for i in 1:length(list)
                    if isfile(list[i])
                        # o.k.
                    elseif isfile(list[i] * ".gz")
                        list[i] = list[i] * ".gz"
                    else
                        filesfound = false
                        break
                    end
                end
                if filesfound
                    return list
                end
            end
        end
    end

    # Take the remote path(s), always uncompressed.
    filename = toc.dataURL * path
    if (occursin("-Ar", path) || occursin("-Zr", path)) && endswith(path, ".g")
        # In case of matrices in char. zero, use the JSON format file,
        # which is currently stored at another web address.
        if entry[:tocID] == "core"
            # In the case of `core` files, remove the directory part from `path`.
            pos = findlast(isequal('/'), path)
            path = SubString(path, pos + 1)
            filename = source_url_json * path
        end
        filename = SubString(filename, 1, length(filename) - 2) * ".json"
    end
    if entry[:filenum] == 1
        return [filename]
    else
        return map(i -> filename * string(i), 1:entry[:filenum])
    end
end


"""
    AtlasGenerators( entry::Dict )

Return the array of permutations or matrices that generate the group
described by `entry`,
which is assumed to be a return value of [`OneAtlasInfo`](@ref)
or an entry in an array returned by `AllAtlasInfos`.

# Examples

```jldoctest
julia> info = OneAtlasInfo( Dict( :name => "A5", :degree => 5 ) )
Dict{Symbol,Any} with 6 entries:
  :id      => ""
  :degree  => 5
  :name    => "A5"
  :path    => "alt/A5/mtx/A5G1-p5B0.m"
  :tocID   => "core"
  :filenum => 2

julia> gens = AtlasGenerators( info )
2-element Array{Perm{Int64},1}:
 (1,2)(3,4)
 (1,3,5)   
```
"""
function AtlasGenerators(entry::Dict)

    filenam = group_atlas_filenames(entry)

    if length(filenam) != entry[:filenum]
        # something went wrong, perhaps a repres. over a number field,
        # with data in GAP format and with no local JSON file available
        error("cannot find filename(s)")
    elseif haskey(entry, :degree)
        # permutation representation
        return [nemo_permutation_from_meataxe_file(filenam[i]) for i in 1:entry[:filenum]]
    elseif haskey(entry, :q)
        # matrix representation over a finite field
        if haskey(entry, :field)
            return [
                nemo_matrix_from_meataxe_file(filenam[i], entry[:field])
                for i in 1:entry[:filenum]
            ]
        else
            return [nemo_matrix_from_meataxe_file(filenam[i]) for i in 1:entry[:filenum]]
        end
    elseif haskey(entry, :ringinfo)
        # matrix representation in char. zero
        if haskey(entry, :gen)
            # Respect the given field generator.
            return nemo_matrices_from_JSON_file(filenam[1], entry[:gen])
        else
            # Create the needed field anew if needed.
            return nemo_matrices_from_JSON_file(filenam[1])
        end
    else
        # this should never happen
        error("unsupported kind of representation")
    end
end

function ordinal(n::Int)
    str = string(n)
    if n < 0
        n = -n
    end
    if mod(n, 10) == 1 && mod(n, 100) != 11
        return str * "st"
    elseif mod(n, 10) == 2 && mod(n, 100) != 12
        return str * "nd"
    elseif mod(n, 10) == 3 && mod(n, 100) != 13
        return str * "rd"
    else
        return str * "th"
    end
end

function repinfoline(info::Dict, toc::AtlasTOC)
    descr = ""
    spl = split(info[:path], "/")
    spl = split(spl[length(spl)], ".")
    file = spl[1]

    if haskey(GroupAtlasData.char, file)
        char = GroupAtlasData.char[file]
        if length(char) == 5
            if char[3] == 0
                descr = "χ = " * char[5]
            else
                descr = "φ = " * char[5]
            end
        end
    end

    if haskey(info, :degree)
        # permutation repres.
        if haskey(GroupAtlasData.api, file)
            api = GroupAtlasData.api[file]
            if api[1] > 1
                # multiply transitive
                descr = string(api[1]) * "-trans."
            elseif api[1] == 1
                # transitive
                descr = "rank " * string(api[2])
            else
                # intransitive
                descr = "orbit lengths " * join(map(string, api[2]), ", ")
            end
            if length(api) >= 4
                # transitive
                if api[4] != "???"
                    # structure of the point stabilizers is known
                    descr = descr * ", on cosets of " * api[4]
                    if api[3] == "prim" && api[5] != "???"
                        # class of max. subgroups for the point stabilizers is known
                        descr = descr * " (" * ordinal(api[5]) * " max.)"
                    end
                elseif api[3] == "prim"
                    descr = descr * ", primitive"
                end
            end
        end
        return ("G ≤ Sym(" * string(info[:degree]) * info[:id] * ")", descr)
    elseif haskey(info, :dim) && haskey(info, :q)
        # matrix repres. over a finite field
        return (
            "G ≤ GL(" * string(info[:dim]) * info[:id] * "," * string(info[:q]) * ")",
            descr,
        )
    elseif haskey(info, :dim) && haskey(info, :ringinfo)
        # matrix repres. over an infinite ring
        if info[:ringinfo] == "Integers"
            return ("G ≤ GL(" * string(info[:dim]) * info[:id] * ",ℤ)", descr)
        else
            return (
                "G ≤ GL(" * string(info[:dim]) * info[:id] * "," * info[:ringinfo] * ")",
                descr,
            )
        end
    end

    # this should not happen
    return "(unknown info)"
end


"""
    DisplayAtlasInfo( name::String = "all" )

Show an overview table of Group Atlas data.
If no `name` is given or if `name` is the string `"all"`,
one line is shown for each group name for which information is available.
If a group name `name` is given,
one line is shown for each available representation for this group.
"""
function DisplayAtlasInfo(name::String = "all")

    tocinfo = AtlasTOCs[1][AtlasTOCs[2]["all"]]
    toc = tocinfo.toc

    if name == "all"
        # list the overview of groups
        widths = [5, 1]
        for nam in GroupAtlasData.allgroupnames
            widths[1] = max(widths[1], length(nam))
            widths[2] = max(widths[2], length(string(length(toc[nam]))))
        end

        # print the header
        println(rpad("group", widths[1]), " │ ", lpad("#", widths[2]))
        println(repeat("─", widths[1]), "─┼─", repeat("─", widths[2]))

        for nam in GroupAtlasData.allgroupnames
            println(rpad(nam, widths[1]), " │ ", lpad(string(length(toc[nam])), widths[2]))
        end
    elseif haskey(toc, name)
        # list overview for this group
        print("Representations for ", name, ":\n")
        print(repeat("─", 21 + length(name)), "\n")
        n = 0
        nlen = length(string(length(toc[name])))
        lines = map(entry -> repinfoline(entry, tocinfo), toc[name])
        clen = maximum(map(x -> length(x[1]), lines)) + 1
        for l in lines
            n = n + 1
            print(lpad(string(n), nlen), ": ", rpad(l[1], clen), l[2], "\n")
        end
    else
        # do nothing
    end
end
