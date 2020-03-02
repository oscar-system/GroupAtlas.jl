
"""
    contents_textfile( filename::String )

Return the string that is the contents of the text file with name `filename`.
If `filename` starts with `http://` then `filename` is treated
as an URL, and the function attempts to download the file.
Otherwise `filename` is assumed to be the name of a (perhaps gzipped)
local file.

# Examples
```jldoctest
julia> testdir = joinpath( dirname( pathof( GroupAtlas ) ), "../test" );

julia> str = GroupAtlas.contents_textfile(
                  joinpath( testdir, "A5G1-f4r2aB0.m1" ) )
" 1     4     2     2\\n10\\n21\\n"

julia> str = GroupAtlas.contents_textfile(
                  joinpath( testdir, "A5G1-f4r2aB0.m1.gz" ) )
" 1     4     2     2\\n10\\n21\\n"

```
"""
function contents_textfile(filename::String)
    if startswith(filename, "http://")
        # remote file
        req = HTTP.request("GET", filename; verbose = 0)
        if req.status != 200
            throw(ArgumentError("file " * filename * " not found"))
        end
        return String(req.body)
    else
        # local file (note that HTTP.request does not support file urls)
        if !isfile(filename)
            throw(ArgumentError("no local file " * filename * "(.gz) found"))
        elseif endswith(filename, ".gz")
            # note that read does not automatically read gzipped files
            stream = CodecZlib.GzipDecompressorStream(open(filename))
            res = read(stream, String)
            close(stream)
            return res
        else
            return read(filename, String)
        end
    end
end


"""
    TuplesIterator( arr, n )

implements the iteration over the `n`-fold Cartesian product
of the array `arr`, where the ordering is lexicographic,
induced by the ordering in `arr`.

This is used in the module `GroupAtlas` for creating element lists
of finite fields that are ordered as defined in the C-MeatAxe.

# Examples
```jldoctest
julia> iter = GroupAtlas.TuplesIterator( [ 1, 2, 3 ], 2 )
GroupAtlas.TuplesIterator([1, 2, 3], 3, 2)

julia> for x in iter print( x, " " ); end
[1, 1] [1, 2] [1, 3] [2, 1] [2, 2] [2, 3] [3, 1] [3, 2] [3, 3] 
```
"""
mutable struct TuplesIterator
    arr::Array{T,1} where {T}
    len::Int
    n::Int

    function TuplesIterator(arr::Array{T,1} where {T}, n::Int)
        return new(arr, length(arr), n)
    end
end

# Initialize the iterator.
function Base.iterate(ci::TuplesIterator)
    if ci.len == 0
        # The iterator is empty.
        return
    end
    return ([ci.arr[1] for i in 1:ci.n], [1 for i in 1:ci.n])
end

# Define the iteration step.
# Note that `state` is changed in place.
function Base.iterate(ci::TuplesIterator, state::Array{Int,1})
    # Find the largest position with a value less than `ci.len`.
    for i in ci.n:-1:1
        if state[i] < ci.len
            state[i] = state[i] + 1
            return ([ci.arr[j] for j in state], state)
        else
            state[i] = 1
        end
    end

    # There is no such position, we are done.
    return
end

function Base.length(ci::TuplesIterator)
    return ci.len^ci.n
end


# order_mod( i::Int, p::Int )
#
# Return the multiplicative order of `i` modulo `p`
# if `i` and `p` are coprime, and `0` otherwise.
#
function order_mod(i::Int, p::Int)
    if gcd(i, p) != 1
        o = 0
    else
        o = 1
        x = i
        while x > 1
            x = mod(x * i, p)
            o = o + 1
        end
    end

    return o
end


# primitive_root_mod_prime( p::Int )
#
# Return the smallest integer whose multiplicative order modulo the
# prime integer `p` is `p-1`.
#
# (The primes that occur in the Group Atlas are small.)
#
function primitive_root_mod_prime(p::Int)
    if p == 2
        return 1
    end

    # Run through all candidates
    for root in 2:(p-1)
        if order_mod(root, p) == p - 1
            return root
        end
    end

    error(p, " is not a prime")
end


"""
    const finitefieldlists = Dict()

We cache the lookup lists for small finite fields.
Extensions over different subfields must have separate entries,
and there are also various possible representations of a finite prime field.
"""
const finitefieldlists = Dict()


@doc Markdown.doc"""
    finitefieldlist( F::Union{Nemo.NmodRing,Nemo.GaloisField,Nemo.FqNmodFiniteField},
                     p::Int, d::Int )

Return the array of elements of the subfield of order `p^d`
in the finite field `F` in characteristic `p`.

The elements are ordered as follows.
The bijection between the elements in the field
with $q = p^d$ elements and the set $\{ 0, 1, \ldots, q-1 \}$
of integers by assigning the field element
$\sum_{{i=0}}^{{d-1}} c_i z^i$ to the integer
$\sum_{{i=0}}^{{d-1}} c_i p^i$,
where the $c_i$ are in the set $\{ 0, 1, \ldots, p-1 \}$
and $z$ is the primitive root of the field with $q$ elements
that corresponds to the residue class of the indeterminate,
modulo the ideal spanned by the Conway polynomial of degree $d$
over the field with $p$ elements.
(This is the ordering defined in the C-MeatAxe.)

# Examples
```jldoctest
julia> F = Nemo.NmodRing( UInt64( 3 ) );

julia> GroupAtlas.finitefieldlist( F, 3, 1 )
3-element Array{Nemo.nmod,1}:
 0
 1
 2

julia> p = 2;  F, z = Nemo.FiniteField( p, 2, "z" );

julia> GroupAtlas.finitefieldlist( F, p, 1 )
2-element Array{fq_nmod,1}:
 0
 1

julia> GroupAtlas.finitefieldlist( F, p, 2 )
4-element Array{fq_nmod,1}:
 0  
 1  
 z  
 z+1

```
""" function finitefieldlist(
    F::Union{Nemo.NmodRing,Nemo.GaloisField,Nemo.FqNmodFiniteField},
    p::Int,
    d::Int,
)
    size = p^d
    Fsize = fieldorder(F)
    if !haskey(finitefieldlists, (F, size))
        if isa(F, Nemo.NmodRing)
            # residue class ring, no prim. root and no `order` method
            if Fsize != size
                error(
                    "cannot regard a matrix over the field with ",
                    size,
                    " elements as a matrix over the field with ",
                    Int(Fsize),
                    " elements",
                )
            end
            root = primitive_root_mod_prime(Fsize) * one(F)
        elseif isa(F, Nemo.GaloisField)
            # prime field, no prim. root
            if mod(Fsize, size) != 0
                error(
                    "given field F of size ",
                    Fsize,
                    " does not contain the field of size ",
                    size,
                )
            end
            root = primitive_root_mod_prime(Fsize) * one(F)
        elseif isa(F, Nemo.FqNmodFiniteField)
            # field defined by a Conway polynomial
            if mod(Fsize, size) != 0
                error(
                    "given field F of size ",
                    Fsize,
                    " does not contain a subfield of size ",
                    size,
                )
            end

            # Here we may have to deal with a proper subfield.
            root = Nemo.gen(F)
            if Fsize != size
                # Handle the embedding of the subfield.
                root = root^(((Fsize - 1) // (size - 1)).num)
            end
        else
            error(
                "F must be in Nemo.NmodRing or Nemo.GaloisField or ",
                "Nemo.FqNmodFiniteField",
            )
        end

        pow = one(root)
        powers = [pow]
        for i in 2:d
            pow = pow * root
            push!(powers, pow)
        end
        powers = reverse(powers)
        itr = TuplesIterator(collect(0:(p-1)), d)

        finitefieldlists[(F, size)] =
            [sum(c * z for (c, z) in zip(elm, powers)) for elm in itr]
    end

    return finitefieldlists[(F, size)]
end


"""
    meataxe_file_header_info( headerline::String )

Return a dictionary with the keys `mode`, `size`, `dimr`, `dimc`
describing the information in `headerline`,
which is assumed to be the header line of a MeatAxe text format file.

# Examples
```jldoctest
julia> GroupAtlas.meataxe_file_header_info( "12 1 100 1" )
Dict{Symbol,Int64} with 4 entries:
  :mode => 12
  :dimc => 1
  :size => 1
  :dimr => 100

julia> GroupAtlas.meataxe_file_header_info( "matrix field=7 rows=5 cols=6" )
Dict{Symbol,Int64} with 4 entries:
  :mode => 1
  :dimc => 6
  :size => 7
  :dimr => 5

```
"""
function meataxe_file_header_info(headerline::String)
    header = split(headerline, " ", keepempty = false)
    if !('=' in headerline)
        # assume old header format, normally four integers
        if length(header) == 3 && header[1] == "12"
            # If the header is valid then we are in the case of permutations
            # of degree requiring at least 6 digits,
            # an example is '12     1100000     1'.
            if header[2][1] == '1'
                push!(header, header[3])
                header[3] = SubString(header[2], 2)
                header[2] = "1"
            end
        end
        if length(header) != 4
            error("corrupted header ", headerline, " of MeatAxe file")
        end

        mode = parse(Int, header[1])
        size = parse(Int, header[2])
        dimr = parse(Int, header[3])
        dimc = parse(Int, header[4])
    else
        # assume new header format, textual information
        if length(header) == 2 &&
           header[1] == "permutation" && startswith(header[2], "degree=")
            mode = 12
            size = 1
            dimr = parse(Int, SubString(header[2], 8))
            dimc = 1
        elseif length(header) == 4 &&
               header[1] == "matrix" && startswith(header[2], "field=") &&
               startswith(header[3], "rows=") && startswith(header[4], "cols=")
            mode = 6
            size = parse(Int, SubString(header[2], 7))
            if size < 10
                mode = 1
            end
            dimr = parse(Int, SubString(header[3], 6))
            dimc = parse(Int, SubString(header[4], 6))
        elseif length(header) == 4 && header[1] == "integer" && header[2] == "matrix" &&
               startswith(header[3], "rows=") && startswith(header[4], "cols=")
            mode = 8
            size = nothing
            dimr = parse(Int, SubString(header[3], 6))
            dimc = parse(Int, SubString(header[4], 6))
        elseif length(header) == 3 && header[1] == "integer-matrix" &&
               startswith(header[2], "rows=") && startswith(header[3], "cols=")
            mode = 8
            size = nothing
            dimr = parse(Int, SubString(header[2], 6))
            dimc = parse(Int, SubString(header[3], 6))
        else
            error("not supported file header")
        end
    end

    return Dict(:mode => mode, :size => size, :dimr => dimr, :dimc => dimc)
end


"""
    nemo_permutation_from_meataxe_file( filename::String, degree::Int = 0 )

Return the Nemo permutation that is encoded by the MeatAxe text format file
with name `filename`,
on the points 1:`degree`, where the default is the degree stored in the file.

# Examples
```jldoctest
julia> testdir = joinpath( dirname( pathof( GroupAtlas ) ), "../test" );

julia> GroupAtlas.nemo_permutation_from_meataxe_file(
           joinpath( testdir, "M11G1-p11B0.m1" ) )
(2,10)(4,11)(5,7)(8,9)

```
"""
function nemo_permutation_from_meataxe_file(filename::String, degree::Int = 0)
    # Try to read the file.
    str = contents_textfile(filename)
    lines = split(str, "\n", keepempty = false)

    # Inspect the file header.
    header = meataxe_file_header_info(string(popfirst!(lines)))

    mode = header[:mode]
    num = header[:size]
    deg = header[:dimr]

    if mode != 12
        error("MeatAxe file ", filename, " does not contain permutations")
    elseif num != 1
        error("MeatAxe file ", filename, " contains more than one permutation?")
    end

    # Turn lines into numbers.
    imgs = map(x -> parse(Int, x), lines)

    # Extend the permutation if necessary.
    if degree != 0
        if deg > degree
            error(
                "the largest moved point of the permutation is ",
                deg,
                ", degree ",
                degree,
                " is not valid",
            )
        elseif deg < degree
            append!(imgs, collect((deg+1):degree))
        end
    end

    return Nemo.perm(imgs)
end


"""
    nemo_matrix_from_meataxe_file( filename::String, F = nothing )

Return the Nemo matrix that is encoded by the MeatAxe text format file
with name `filename`, over the finite field `F`,
which may have been constructed with `Nemo.FiniteField` or `Nemo.NmodRing`
or `Nemo.GaloisField`.
If `F` is not given then it is constructed with `Nemo.FiniteField`.

# Examples
```jldoctest
julia> testdir = joinpath( dirname( pathof( GroupAtlas ) ), "../test" );

julia> GroupAtlas.nemo_matrix_from_meataxe_file(
           joinpath( testdir, "A5G1-f4r2aB0.m1" ) )
[1  0]
[z  1]

julia> F, z = Nemo.FiniteField( 2, 4, "z" );

julia> GroupAtlas.nemo_matrix_from_meataxe_file(
           joinpath( testdir, "A5G1-f4r2aB0.m1" ), F )
[    1  0]
[z^2+z  1]

```
"""
function nemo_matrix_from_meataxe_file(filename::String, F = nothing)

    # Try to read the file.
    str = contents_textfile(filename)
    lines = split(str, "\n", keepempty = false)

    # Inspect the file header.
    header = meataxe_file_header_info(string(popfirst!(lines)))

    mode = header[:mode]
    size = header[:size]
    dimr = header[:dimr]
    dimc = header[:dimc]

    if isnothing(size)
        # This happens only for integer matrices in MeatAxe format (mode 8)
        if isnothing(F)
            F = Nemo.ZZ
        end
    else
        factors = collect(Primes.factor(size))
        if length(factors) != 1
            error("size of the field is not a prime power?")
        end
        power = factors[1]
        char = power[1]
        degree = power[2]

        if isnothing(F)
            # No field was given, construct one of the given size.
            # Assume that the relevant Conway polynomial is available.
            F, z = Nemo.FiniteField(char, degree, "z")
        end
    end

    T = typeof(one(F))
    l = T[]

    if mode == 1
        # fixed format without whitespace (field of size at most 9)
        ffl = finitefieldlist(F, char, degree)
        lookup = Dict{Char,Any}()
        for i in 0:(size-1)
            lookup[string(i)[1]] = ffl[i+1]
        end
        for line in lines
            append!(l, [lookup[char] for char in line])
        end
    elseif mode == 5
        # matrix of integers to be reduced modulo 'char'
        Fone = one(F)
        for line in lines
            append!(
                l,
                [parse(Int, obj) * Fone for obj in split(line, " ", keepempty = false)],
            )
        end
    elseif mode in 3:6
        # We assume that the field is not large,
        # compared to, e. g., the number of matrix entries
        ffl = finitefieldlist(F, char, degree)
        lookup = Dict{String,Any}()
        for i in 0:(size-1)
            lookup[string(i)] = ffl[i+1]
        end
        for line in lines
            append!(l, [lookup[obj] for obj in split(line, " ", keepempty = false)])
        end
    elseif mode == 8
        for line in lines
            append!(
                l,
                [parse(Nemo.fmpz, obj) for obj in split(line, " ", keepempty = false)],
            )
        end
    else
        error("mode must be one of 1, 3, 4, 5, 6, 8")
#T support also mode = 2 (a permutation, to be converted to a matrix)
    end

    return Nemo.matrix(F, dimr, dimc, l)
end


#############################################################################
##
##  Create MeatAxe text format strings from Nemo permutations and matrices.
##

"""
    meataxe_string( perm::Nemo.Perm, degree::Int = 0,
                    headerformat::String = "numeric" )

Return a string in MeatAxe text format that encodes the permutation `perm`,
viewed as a permutation of degree `degree`.
If `degree` is zero then the degree of `perm` is taken.

The admissible `headerformat` values are "numeric", "numeric (fixed)",
and "textual".

# Examples
```jldoctest
julia> GroupAtlas.meataxe_string( 
                     Nemo.perm( [1,10,3,11,7,6,5,9,8,2,4] ) )
"12 1 11 1\\n1\\n10\\n3\\n11\\n7\\n6\\n5\\n9\\n8\\n2\\n4\\n"

julia> GroupAtlas.meataxe_string( 
                     Nemo.perm( [1,10,3,11,7,6,5,9,8,2,4] ), 13 )
"12 1 13 1\\n1\\n10\\n3\\n11\\n7\\n6\\n5\\n9\\n8\\n2\\n4\\n12\\n13\\n"

```
"""
function meataxe_string(perm::Nemo.Perm, degree::Int = 0, headerformat::String = "numeric")

    imgs = perm.d
    if degree == 0
        degree = length(imgs)
    elseif degree < length(imgs)
        error("degree must be at least ", length(imgs))
    end

    if headerformat == "numeric"
        str = "12 1 " * string(degree) * " 1\n"
    elseif headerformat == "numeric (fixed)"
        str = "    12     1" * lpad(degree, 6) * "     1 \n"
    elseif headerformat == "textual"
        str = "permutation degree=" * string(degree) * "\n"
    else
        error(
            "`headerformat` must be one of \"numeric\", \"numeric (fixed)\"",
            ", or \"textual\"",
        )
    end

    for img in imgs
        str = str * string(img) * "\n"
    end

    if degree > length(imgs)
        for i in (length(imgs)+1):degree
            str = str * string(i) * "\n"
        end
    end

    return str
end


"""
    meataxe_string( mat::Union{Nemo.fq_nmod_mat,Nemo.nmod_mat,Nemo.gfp_mat},
                    headerformat::String = "numeric" )

Return a string in MeatAxe text format that encodes the matrix `mat`.
The admissible `headerformat` values are "numeric", "numeric (fixed)",
and "textual".

# Examples
```jldoctest
julia> arr = [ 1 2 3 ; 4 5 6 ; 7 8 9 ];

julia> F, z = Nemo.FiniteField( 11, 1, "z" );

julia> GroupAtlas.meataxe_string( Nemo.matrix( F, arr ) )
"6 11 3 3\\n1\\n2\\n3\\n4\\n5\\n6\\n7\\n8\\n9\\n"

julia> F = Nemo.NmodRing( UInt64( 5 ) );

julia> GroupAtlas.meataxe_string( Nemo.matrix( F, arr ) )
"1 5 3 3\\n123\\n401\\n234\\n"

julia> F = Nemo.GaloisField( UInt64( 5 ) );

julia> GroupAtlas.meataxe_string( Nemo.matrix( F, arr ) )
"1 5 3 3\\n123\\n401\\n234\\n"

```
"""
function meataxe_string(
    mat::Union{Nemo.fq_nmod_mat,Nemo.nmod_mat,Nemo.gfp_mat},
    headerformat::String = "numeric",
)
    if Nemo.nrows(mat) == 0 || Nemo.ncols(mat) == 0
        error("MeatAxe matrices must have at least one row and column")
    end
    F = mat.base_ring
    q = fieldorder(F)
    if q < 10
        mode = "1"
    else
        mode = "6"
    end
    factors = collect(Primes.factor(q))
    if length(factors) != 1
        error("size of the field is not a prime power?")
    end

    power = factors[1]
    char = power[1]
    degree = power[2]

    if headerformat == "numeric"
        str =
            string(mode) *
            " " *
            string(q) *
            " " *
            string(Nemo.nrows(mat)) *
            " " *
            string(Nemo.ncols(mat)) *
            "\n"
    elseif headerformat == "numeric (fixed)"
        str =
            lpad(mode, 6) *
            lpad(string(q), 6) *
            lpad(string(Nemo.nrows(mat))) *
            lpad(string(Nemo.ncols(mat))) *
            "\n"
    elseif headerformat == "textual"
        "matrix field=" * string(q) * " rows=" * string(Nemo.nrows(mat)) * " cols=",
        string(Nemo.ncols(mat)) * "\n"
    else
        error(
            "`headerformat` must be one of \"numeric\", \"numeric (fixed)\"",
            ", or \"textual\"",
        )
    end

    ffl = finitefieldlist(F, char, degree)
    lookup = Dict{Any,String}()
    for i in 0:(q-1)
        lookup[ffl[i+1]] = string(i)
    end

    if mode == "1"
        for i in 1:Nemo.nrows(mat)
            len = 0
            for j in 1:Nemo.ncols(mat)
                str = str * lookup[mat[i, j]]
                len = len + 1
                if len == 80
                    str = str * "\n"
                    len = 0
                end
            end
            if len != 0
                str = str * "\n"
            end
        end
    else  #  mode == "6" (free format for all q > 9)
        for i in 1:Nemo.nrows(mat)
            for j in 1:Nemo.ncols(mat)
                str = str * lookup[mat[i, j]]
                str = str * "\n"
            end
        end
    end

    return str
end
