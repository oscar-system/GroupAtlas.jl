
"""
    unite_dictionaries!( target::Dict, source::Dict )

Set all values from `source` in `target`,
such that the values of keys available in both `target` and `source`
get concatenated with `Base.cat`.

# Examples
```jldoctest
julia> target = Dict( "a" => 1, "b" => [ 1 ] )
Dict{String,Any} with 2 entries:
  "b" => [1]
  "a" => 1

julia> source = Dict( "a" => 0, "c" => 2, "b" => [ 2 ] )
Dict{String,Any} with 3 entries:
  "c" => 2
  "b" => [2]
  "a" => 0

julia> GroupAtlas.unite_dictionaries!( target, source )

julia> target
Dict{String,Any} with 3 entries:
  "c" => 2
  "b" => [1, 2]
  "a" => [1, 0]

```
"""
function unite_dictionaries!( target::Dict, source::Dict )
    for k in keys( source )
      if haskey( target, k )
        target[k] = Base.cat( target[k], source[k]; dims = 1 )
      else
        target[k] = source[k]
      end
    end
end


"""
    scanfilename( path::String )

Return a dictionary that describes the Group Atlas data given by the
data file with the local path `path`.

# Examples
```jldoctest
julia> GroupAtlas.scanfilename( "somewhere/M11G1-p11B0.m1" )
Dict{Symbol,Any} with 4 entries:
  :id     => ""
  :degree => 11
  :name   => "M11"
  :path   => "somewhere/M11G1-p11B0.m1"

julia> GroupAtlas.scanfilename( "somewhere/A5G1-f4r2aB0.m1" )
Dict{Symbol,Any} with 5 entries:
  :dim  => 2
  :id   => "a"
  :path => "somewhere/A5G1-f4r2aB0.m1"
  :name => "A5"
  :q    => 4

julia> GroupAtlas.scanfilename( "somewhere/M11G1-Zr10aB0.g" )
Dict{Symbol,Any} with 6 entries:
  :dim      => 10
  :ringinfo => "Integers"
  :id       => "a"
  :path     => "somewhere/M11G1-Zr10aB0.g"
  :name     => "M11"
  :ring     => Integer Ring

julia> GroupAtlas.scanfilename( "somewhere/M11G1-Ar10bB0.g" )
Dict{Symbol,Any} with 6 entries:
  :dim       => 10
  :ringinfo  => "Field([Sqrt(-2)])"
  :id        => "b"
  :path      => "somewhere/M11G1-Ar10bB0.g"
  :name      => "M11"
  :polcoeffs => Any[2, 0, 1]

```
"""
function scanfilename( path::String )
    # cut off the prefix until the last '/'
    spl = split( path, "/" )

    # cut at the unique '-'
    pair = split( spl[ end ], "-" )

    # extract the group name (as used in filenames)
    name = pair[1]
    i = length( name )
    while isdigit( name[i] )
      i = i-1
    end
    if ! ( name[i] == 'G' || name[i] == 'W' )
      error( "corrupted entry ", path )
    end
    name = name[1:(i-1)]
    if haskey( GroupAtlasData.angn, name )
      # replace the name used in the filename by the group name.
      name = GroupAtlasData.angn[ name ]
    end

    # distinguish representation types
    info = pair[2]
    entry = nothing
    if info[1] == 'p' && isdigit( info[2] )
      # permutation, extract the degree and the id
      i = 2
      while isdigit( info[i] )
        i = i+1
      end
      j = i
      while info[j] != 'B'
        j = j+1
      end
      degree = parse( Int, info[2:(i-1)] )
      entry = Dict( :name => name,
                    :degree => degree,
                    :id => info[i:(j-1)],
                    :path => path )
    elseif info[1] == 'f' && isdigit( info[2] )
      # matrix over a finite field, extract field size, dimension, id
      i = 2
      while isdigit( info[i] )
        i = i+1
      end
      q = parse( Int, info[2:(i-1)] )
      j = i+1
      while isdigit( info[j] )
        j = j+1
      end
      dim = parse( Int, info[(i+1):(j-1)] )
      i = j
      while info[i] != 'B'
        i = i+1
      end
      id = info[j:(i-1)]
      entry = Dict( :name => name,
                    :q => q,
                    :dim => dim,
                    :id => id,
                    :path => path )
    elseif info[1] == 'Z' && info[2] == 'r' && isdigit( info[3] )
      # matrix over the integers, extract dimension and id
      i = 3
      while isdigit( info[i] )
        i = i+1
      end
      dim = parse( Int, info[3:(i-1)] )
      j = i
      while info[j] != 'B'
        j = j+1
      end
      entry = Dict( :name => name,
                    :ringinfo => "Integers",
                    :ring => Nemo.ZZ,
                    :dim => dim,
                    :id => info[i:(j-1)],
                    :path => path )
    elseif info[1] == 'A' && info[2] == 'r' && isdigit( info[3] )
      # matrix over a number field, extract dimension and id
      i = 3
      while isdigit( info[i] )
        i = i+1
      end
      dim = parse( Int, info[3:(i-1)] )
      j = i
      while info[j] != 'B'
        j = j+1
      end
      key = SubString( path, ( findlast( "/", path )[1]+1 ):
                             ( findfirst( ".", path )[1]-1 ) )
      if haskey( GroupAtlasData.rng, key )
        # We need the information about the defining polynomial.
        ringinfo = GroupAtlasData.rng[ key ]
        entry = Dict( :name => name,
                      :ringinfo => ringinfo[1],
                      :dim => dim,
                      :id => info[i:(j-1)],
                      :path => path )
        if length( ringinfo ) == 3
          entry[ :polcoeffs ] = ringinfo[3]
        end
      end
    end

    return entry
end


"""
    isless_groupnames( nam1::String, nam2::String )

Compare the two strings `nam1`, `nam2`,
and return `true` iff `nam1` shall come before `nam2` in the global array
`GroupAtlasData.allgroupnames`.
(The analogon in GAP is `AGR.CompareAsNumbersAndNonnumbers`.)

# Examples
```jldoctest
julia> GroupAtlas.isless_groupnames( "A5", "A6" )
true

julia> GroupAtlas.isless_groupnames( "A10", "A9" )
false

julia> GroupAtlas.isless_groupnames( "A5", "Co1" )
true
```
"""
function isless_groupnames( nam1::String, nam2::String )
    len1 = length( nam1 )
    len2 = length( nam2 )
    len = len1
    if len2 < len
      len = len2
    end
    digit = false
    comparenumber = 0
    for i in 1:len
      if isdigit( nam1[i] )
        if isdigit( nam2[i] )
          digit = true
          if comparenumber == 0
            # first digit of a number, or previous digits were equal
            if nam1[i] < nam2[i]
              comparenumber = 1
            elseif nam1[i] != nam2[i]
              comparenumber = -1
            end
          end
        else
          # if digit then the current number in 'nam2' is shorter,
          # so 'nam2' is smaller;
          # if not digit then a number starts in 'nam1' but not in 'nam2',
          # so 'nam1' is smaller
          return ! digit
        end
      elseif isdigit( nam2[i] )
        # if digit then the current number in 'nam1' is shorter,
        # so 'nam1' is smaller;
        # if not digit then a number starts in 'nam2' but not in 'nam1',
        # so 'nam2' is smaller
        return digit
      else
        # both characters are non-digits
        if digit
          # first evaluate the current numbers (which have the same length)
          if comparenumber == 1
            # nam1 is smaller
            return true
          elseif comparenumber == -1
            # nam2 is smaller
            return false
          end
          digit = false
        end
        # now compare the non-digits
        if nam1[i] != nam2[i]
          return nam1[i] < nam2[i]
        end
      end
    end

    if digit
      # The suffix of the shorter string is a number.
      # If the longer string continues with a digit then it is larger,
      # otherwise the first digits of the number decide.
      if len < len1 && isdigit( nam1[ len+1 ] )
        # nam2 is smaller
        return false
      elseif len < len2 && isdigit( nam2[ len+1 ] )
        # nam1 is smaller
        return true
      elseif comparenumber == 1
        # nam1 is smaller
        return true
      elseif comparenumber == -1
        # nam2 is smaller
        return false
      end
    end

    # Now the longer string is larger.
    return len1 < len2
end


"""
    isless_entries_for_group( entry1::Dict, entry2::Dict )

Return `true` or `false`, depending on whether `entry1` shall appear
before `entry2` in an [`AtlasTOC`](@ref) object that contains both of them;
this defines the ordering of lines shown by [`DisplayAtlasInfo`](@ref).

- Permutation representations precede matrix representations,
  ordered by increasing degree.
- Then follow matrix representations over finite fields,
  ordered by increasing field size and dimension of the matrices.
- Then follow matrix representations over the integers,
  ordered by increasing dimension of the matrices.
- Then follow matrix representations over number fields,
  ordered by increasing dimension of the matrices.

# Examples
```jldoctest
julia> info1 = OneAtlasInfo( Dict( :name => "A5", :degree => 5 ) );

julia> info2 = OneAtlasInfo( Dict( :name => "A5", :q => 2 ) );

julia> GroupAtlas.isless_entries_for_group( info1, info2 )
true

```
"""
function isless_entries_for_group( entry1::Dict, entry2::Dict )
    # perm. repr. come first, ordered by degree
    if haskey( entry1, :degree )
      if haskey( entry2, :degree )
        return entry1[ :degree ] < entry2[ :degree ] ||
               ( entry1[ :degree ] == entry2[ :degree ] &&
                 entry1[ :id ] < entry2[ :id ] )
      else
        return true
      end
    elseif haskey( entry2, :degree )
      return false
    end
    # repr. over finite fields come next, ordered by field size and dimension
    if haskey( entry1, :q )
      if haskey( entry2, :q )
        return entry1[ :q ] < entry2[ :q ] ||
               ( ( entry1[ :q ] == entry2[ :q ] ) && 
                 ( entry1[ :dim ] < entry2[ :dim ] ||
                   ( entry1[ :dim ] == entry2[ :dim ] &&
                     entry1[ :id ] < entry2[ :id ] ) ) )
      else
        return true
      end
    elseif haskey( entry2, :q )
      return false
    end
    # repr. over the integers come next, ordered by dimension
    if haskey( entry1, :ring ) && entry1[ :ring ] == Nemo.ZZ
      if haskey( entry2, :ring ) && entry2[ :ring ] == Nemo.ZZ
        return entry1[ :dim ] < entry2[ :dim ] ||
               ( entry1[ :dim ] == entry2[ :dim ] &&
                 entry1[ :id ] < entry2[ :id ] )
      else
        return true
      end
    elseif haskey( entry2, :ring ) && entry2[ :ring ] == Nemo.ZZ
      return false
    end
    # other repr. come next, ordered by dimension
    return entry1[ :dim ] < entry2[ :dim ] ||
           ( entry1[ :dim ] == entry2[ :dim ] &&
             entry1[ :id ] < entry2[ :id ] )
end


mutable struct GroupAtlasDataObject
    source_urls::Array{String,1}
    allgroupnames::Array{String,1}
    angn::Dict{String,String}
    gnan::Dict{String,String}
    rng::Dict{String,Any}
    api::Dict{String,Any}
    char::Dict{String,Any}
end

"""
    GroupAtlasData

is a global object for storing information read from the table of contents
files of the Group Atlas, with the following fields.

`source_urls`:
    an array of the addresses of known overview files for the Group Atlas,
    each entry is either a web URL or a local directory path.

`allgroupnames`:
    a sorted array of all group names for which the Group Atlas provides
    information, independent of the individual tables of contents.

`angn`:
    a dictionary for mapping group names that appear in filenames
    to the group names used in user functions;
    for example, the value at `"Sz8d3"` is `"Sz(8).3"`,

`gnan`:
    a dictionary for mapping the values in `angn` to their keys,

`rng`:
    a dictionary for mapping the filenames of matrix representations over
    number fields to information about the field of matrix entries,

`api`:
    a dictionary for mapping the filenames of permutation representations
    to information about the representation,

`char`:
    a dictionary for mapping the filenames of matrix representations
    to information about the character of the representation.
"""
const GroupAtlasData = GroupAtlasDataObject(
          String[],
          String[],
          Dict{String,String}(),
          Dict{String,String}(),
          Dict{String,Any}(),
          Dict{String,Any}(),
          Dict{String,Any}() )

"""
    AtlasTOC( url::String )

Return an object that describes a collection of Group Atlas data.
`url` can be either an URL of a JSON format file that describes
the contents of some data directory,
or the string `"all"`, which means the union of information
of the `AtlasTOC` objects that belong to the URLs listed in
the `source_urls` field of [`GroupAtlasData`](@ref).
"""
mutable struct AtlasTOC
    tocURL::String
    dataURL::String
    localDirectory::String
    ID::String
    toc::Dict{String,Any}

    # default: take all sources
    function AtlasTOC( url::String = "all" )

      # Check whether we have stored this already.
      if ! haskey( AtlasTOCs[1], url )

        toc = Dict{String,Any}()

        if url == "all"
          # Join the contents given by `GroupAtlasData.source_urls`.
          notified = Dict()
          tocID = url
          dataURL = "dummy"
          localDirectory = ""
          alltocs = map( AtlasTOC, GroupAtlasData.source_urls )

          for onetoc in alltocs
            # Make sure that each part gets notified only once.
            if ! haskey( notified, onetoc.ID )
              unite_dictionaries!( toc, onetoc.toc )
              notified[ onetoc.ID ] = true
            end
          end

        else
          # Evaluate one JSON format file (local or remote).
          str = contents_textfile( url )
          if isnothing( str )
            error( "could not read file at '", url, "'" )
          end

          prs = JSON.parse( str; dicttype = Dict{Symbol,Any} )
          tocID = prs[ :ID ]
          dataURL = prs[ :DataURL ]
          localDirectory = ""
          if isdefined( Main, :GAP )
            if haskey( prs, :LocalDirectory )
              for path in GAP.gap_to_julia( Array{String,1},
                                            GAP.Globals.GAPInfo.RootPaths )
                dir = path * "/pkg/" * prs[ :LocalDirectory ] * "/"
                if isdir( dir )
                  # *All* data files are (assumed to be) available
                  # in this directory.
                  localDirectory = dir
                  break
                end
              end
            elseif tocID == "core"
              for path in GAP.gap_to_julia( Array{String,1},
                                            GAP.Globals.GAPInfo.RootPaths )
                # *Some* data files may be stored locally.
                # We prefer the local files if they are available.
                dir = path * "/pkg/atlasrep/datagens/"
                if isdir( dir )
                  localDirectory = dir
                  break
                end
              end
            end
          end
          allTOC = []

          for entry in prs[ :Data ]
            if entry[1] == "GNAN"
              # declare a group name
              GroupAtlasData.gnan[ entry[2][1] ] = entry[2][2]
              GroupAtlasData.angn[ entry[2][2] ] = entry[2][1]
              toc[ entry[2][1] ] = []
              push!( GroupAtlasData.allgroupnames, entry[2][1] )
            elseif entry[1] == "RNG"
              # notify the ring of definition for an entry in char. zero;
              # the values of the 2nd entry are
              # - a textual description of the field,
              # - a programmatic description of the field,
              # - the coefficient list of a polynomial defining the field.
              GroupAtlasData.rng[ entry[2][1] ] =
                  [ entry[2][i] for i in 2:length( entry[2] ) ]
            elseif entry[1] == "API"
              # information on permutation representations
              GroupAtlasData.api[ entry[2][1] ] = entry[2][2]
            elseif entry[1] == "CHAR"
              # information on characters of matrix representations
              GroupAtlasData.char[ entry[2][2] ] = entry[2]
            elseif entry[1] == "TOC"
              # notify an entry for a group,
              # but wait until 'gnan' and 'rng' are complete
              push!( allTOC, entry[2] )
            end
          end

          for entry in allTOC
            newentry = scanfilename( entry[2] )
            if ! isnothing( newentry )
              newentry[ :tocID ] = tocID
#T as soon as the crc values can be used in Julia, store them
#             newentry[ :crc ] = entry[3]
              newentry[ :filenum ] = length( entry[3] )
              name = newentry[ :name ]
              if ! haskey( toc, name )
                toc[ name ] = []
              end
              push!( toc[ name ], newentry )
            end
          end

        end

        # Sort the list for each group.
        for nam in keys( toc )
          sort!( toc[ nam ], lt = isless_entries_for_group )
        end

        # Sort the list of group names.
        sort!( GroupAtlasData.allgroupnames, lt = isless_groupnames )

        AtlasTOCs[1][ url ] = new( url, dataURL, localDirectory, tocID, toc )
        AtlasTOCs[2][ tocID ] = url

      end

      return AtlasTOCs[1][ url ]
    end
end

# Cache the tables of contents.
const AtlasTOCs = ( Dict{String,AtlasTOC}(), Dict{String,String}() )

