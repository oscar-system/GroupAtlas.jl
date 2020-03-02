using Nemo

@testset "element lists of finite fields" begin
    # iterator of tuples
    itr = GroupAtlas.TuplesIterator( [ 1, 2, 3 ], 0 )
    @test Base.length( itr ) == 1
    @test collect( itr ) == [ [] ]

    itr = GroupAtlas.TuplesIterator( [ 1, 2, 3 ], 1 )
    @test Base.length( itr ) == 3
    @test collect( itr ) == [ [1], [2], [3] ]

    itr = GroupAtlas.TuplesIterator( [ 1, 2, 3 ], 2 )
    @test Base.length( itr ) == 9
    @test collect( itr ) ==
          [ [1,1], [1,2], [1,3], [2,1], [2,2], [2,3], [3,1], [3,2], [3,3] ]

    # order mod prime, prom. root mod prime
    @test [ GroupAtlas.order_mod( i, 7 ) for i in 0:6 ] ==
          [ 0, 1, 3, 6, 3, 6, 2 ]
    @test map( GroupAtlas.primitive_root_mod_prime, [ 2, 3, 5, 7, 11 ] ) ==
          [ 1, 2, 2, 3, 2 ]

    # element lists ...
    # ... for prime fields ...
    p = 7
    F = Nemo.NmodRing( UInt64( p ) )
    @test GroupAtlas.finitefieldlist( F, p, 1 ) ==
          [ one( F ) * i for i in 0:(p-1) ]

    F = Nemo.GaloisField( UInt64( p ) )
    @test GroupAtlas.finitefieldlist( F, p, 1 ) ==
          [ one( F ) * i for i in 0:(p-1) ]
    F, z = Nemo.FiniteField( p, 1, "z" )
    @test GroupAtlas.finitefieldlist( F, p, 1 ) ==
          [ one( F ) * i for i in 0:(p-1) ]

    # ... and for extension fields
    F, z = Nemo.FiniteField( p, 2, "z" )
    @test GroupAtlas.finitefieldlist( F, p, 1 ) ==
          [ one( F ) * i for i in 0:(p-1) ]
    root = z
    @test GroupAtlas.finitefieldlist( F, p, 2 ) ==
          [ i*z + j for i in 0:(p-1) for j in 0:(p-1) ]
end

@testset "interpret MeatAxe text format file headers" begin
    # permutation
    @test GroupAtlas.meataxe_file_header_info( "12 1 100 1" ) ==
          Dict( :mode => 12, :size => 1, :dimr => 100, :dimc => 1 )
    @test GroupAtlas.meataxe_file_header_info(
              "    12     1   100     1" ) ==
          Dict( :mode => 12, :size => 1, :dimr => 100, :dimc => 1 )
    @test GroupAtlas.meataxe_file_header_info( "12 1100 1" ) ==
          Dict( :mode => 12, :size => 1, :dimr => 100, :dimc => 1 )
    @test GroupAtlas.meataxe_file_header_info(
              "permutation degree=100" ) ==
          Dict( :mode => 12, :size => 1, :dimr => 100, :dimc => 1 )

    # matrix over a finite field
    @test GroupAtlas.meataxe_file_header_info( "1 7 5 6" ) ==
          Dict( :mode => 1, :size => 7, :dimr => 5, :dimc => 6 )
    @test GroupAtlas.meataxe_file_header_info(
              "matrix field=7 rows=5 cols=6" ) ==
          Dict( :mode => 1, :size => 7, :dimr => 5, :dimc => 6 )
    @test GroupAtlas.meataxe_file_header_info( "6 11 5 6" ) ==
          Dict( :mode => 6, :size => 11, :dimr => 5, :dimc => 6 )

    # integer matrix
    @test GroupAtlas.meataxe_file_header_info(
              "integer matrix rows=5 cols=6" ) ==
          Dict( :mode => 8, :size => nothing, :dimr => 5, :dimc => 6 )
    @test GroupAtlas.meataxe_file_header_info(
              "integer-matrix rows=5 cols=6" ) ==
          Dict( :mode => 8, :size => nothing, :dimr => 5, :dimc => 6 )

    # error
    @test_throws ErrorException GroupAtlas.meataxe_file_header_info(
                                    "matrix field=7" )
end

@testset "read MeatAxe text format files" begin
    # permutations
    @test GroupAtlas.nemo_permutation_from_meataxe_file(
              abspath( @__DIR__, "M11G1-p11B0.m1" ) ) ==
          Nemo.perm( [1,10,3,11,7,6,5,9,8,2,4] )
    @test GroupAtlas.nemo_permutation_from_meataxe_file(
              abspath( @__DIR__, "M11G1-p11B0.m1" ), 11 ) ==
          Nemo.perm( [1,10,3,11,7,6,5,9,8,2,4] )
    @test GroupAtlas.nemo_permutation_from_meataxe_file(
              abspath( @__DIR__, "M11G1-p11B0.m1" ), 12 ) ==
          Nemo.perm( [1,10,3,11,7,6,5,9,8,2,4,12] )

    # errors for permutations
    @test_throws ErrorException GroupAtlas.nemo_permutation_from_meataxe_file(
              abspath( @__DIR__, "M11G1-f11r11B0.m1" ) )
    @test_throws ErrorException GroupAtlas.nemo_permutation_from_meataxe_file(
              abspath( @__DIR__, "M11G1-p11B0.m1" ), 10 )

    # matrix over a finite field (q < 10)


    # matrix over a finite field (reduce mod p)


    # matrix over a finite field (q > 10)
    arr = [ 10 0 0 ; 0 10 0 ; 9 9 1 ]
    F, z = Nemo.FiniteField( 11, 1, "z" )
    @test GroupAtlas.nemo_matrix_from_meataxe_file(
              abspath( @__DIR__, "L211d2G1-f11r3B0.m1" ) ) ==
    Nemo.matrix( F, arr )
    @test GroupAtlas.nemo_matrix_from_meataxe_file(
              abspath( @__DIR__, "L211d2G1-f11r3B0.m1" ), F ) ==
    Nemo.matrix( F, arr )

    F = Nemo.NmodRing( UInt64( 11 ) )
    @test GroupAtlas.nemo_matrix_from_meataxe_file(
              abspath( @__DIR__, "L211d2G1-f11r3B0.m1" ), F ) ==
    Nemo.matrix( F, arr )

    F = Nemo.GaloisField( UInt64( 11 ) )
    @test GroupAtlas.nemo_matrix_from_meataxe_file(
              abspath( @__DIR__, "L211d2G1-f11r3B0.m1" ), F ) ==
    Nemo.matrix( F, arr )

# matrix over the integers (mode 8)

end

@testset "create MeatAxe text format strings" begin
    # permutation
    @test GroupAtlas.meataxe_string(
              Nemo.perm( [1,10,3,11,7,6,5,9,8,2,4] ) ) ==
          "12 1 11 1\n1\n10\n3\n11\n7\n6\n5\n9\n8\n2\n4\n"
    @test GroupAtlas.meataxe_string(
              Nemo.perm( [1,10,3,11,7,6,5,9,8,2,4] ), 13 ) ==
          "12 1 13 1\n1\n10\n3\n11\n7\n6\n5\n9\n8\n2\n4\n12\n13\n"

    # matrix over a finite field (q > 10)
    F, z = Nemo.FiniteField( 11, 1, "z" )
    arr = [ 10 0 0 ; 0 10 0 ; 9 9 1 ]
    mat = Nemo.matrix( F, arr )
    @test GroupAtlas.meataxe_string( mat ) ==
              "6 11 3 3\n10\n0\n0\n0\n10\n0\n9\n9\n1\n"

    F = Nemo.NmodRing( UInt64( 11 ) )
    mat = Nemo.matrix( F, arr )
    @test GroupAtlas.meataxe_string( mat ) ==
              "6 11 3 3\n10\n0\n0\n0\n10\n0\n9\n9\n1\n"

    F = Nemo.GaloisField( UInt64( 11 ) )
    mat = Nemo.matrix( F, arr )
    @test GroupAtlas.meataxe_string( mat ) ==
              "6 11 3 3\n10\n0\n0\n0\n10\n0\n9\n9\n1\n"

    # matrix over a finite field (q < 10)
    arr = [ 1 2 3 ; 4 5 6 ; 7 8 9 ]
    F, z = Nemo.FiniteField( 5, 1, "z" )
    mat = Nemo.matrix( F, arr )
    @test GroupAtlas.meataxe_string( mat ) ==
              "1 5 3 3\n123\n401\n234\n"

    F = Nemo.NmodRing( UInt64( 5 ) )
    mat = Nemo.matrix( F, arr )
    @test GroupAtlas.meataxe_string( mat ) ==
              "1 5 3 3\n123\n401\n234\n"

    F = Nemo.GaloisField( UInt64( 5 ) )
    mat = Nemo.matrix( F, arr )
    @test GroupAtlas.meataxe_string( mat ) ==
              "1 5 3 3\n123\n401\n234\n"

    # error
    @test_throws ErrorException GroupAtlas.meataxe_string(
              Nemo.perm( [1,10,3,11,7,6,5,9,8,2,4] ), 10 )
# error for matrix!
end

