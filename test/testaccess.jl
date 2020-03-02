
@testset "permutation representations" begin

    info = OneAtlasInfo( Dict( :name => "A5", :degree => 5 ) )
    gens = AtlasGenerators( info )
    @test order( gens[1] ) == 2
    @test order( gens[2] ) == 3
    @test order( gens[1]*gens[2] ) == 5

    info = OneAtlasInfo( Dict( :name => "A5", :degree => 4 ) )
    @test isnothing( info )

end

@testset "matrix representations over finite fields" begin

    info = OneAtlasInfo( Dict( :name => "A5", :char => 2, :dim => 4 ) )
    gens = AtlasGenerators( info )
    @test isone( gens[1]^2 )
    @test isone( gens[2]^3 )
    @test isone( (gens[1]*gens[2])^5 )
    @test mod( Nemo.order( gens[1].base_ring ), 2 ) == 0

    F = Nemo.GaloisField( UInt(2) )
    info = OneAtlasInfo( Dict( :name => "A5", :field => F, :dim => 4 ) )
    gens = AtlasGenerators( info )
    @test isone( gens[1]^2 )
    @test isone( gens[2]^3 )
    @test isone( (gens[1]*gens[2])^5 )
    @test parent( gens[1][1,1] ) == F

    info = OneAtlasInfo( Dict( :name => "A5", :char => 101 ) )
    @test isnothing( info )

end

@testset "repres. over the integers" begin

    info = OneAtlasInfo( Dict( :name => "A5", :char => 0, :dim => 4 ) )
    gens = AtlasGenerators( info )
    @test isone( gens[1]^2 )
    @test isone( gens[2]^3 )
    @test isone( (gens[1]*gens[2])^5 )
    @test gens[1].base_ring === Nemo.ZZ

end

@testset "repres. over number fields" begin

    # field not prescribed
    info = OneAtlasInfo( Dict( :name => "A5",
                               :ringinfo => "Field([Sqrt(5)])" ) )
    gens = AtlasGenerators( info )
    @test isone( gens[1]^2 )
    @test isone( gens[2]^3 )
    @test isone( (gens[1]*gens[2])^5 )

    # field prescribed, but in fact integral matrices
    R,x = PolynomialRing( Nemo.QQ, "x" )
    pol = R([-1,1,1] )
    F,z = NumberField( pol, "z" )
    info = OneAtlasInfo( Dict( :name => "A5", :dim => 4, :fieldgen => z ) )
    gens = AtlasGenerators( info )
    @test gens[1].base_ring === F

    # proper extension field prescribed
    info = OneAtlasInfo( Dict( :name => "A5", :dim => 3, :fieldgen => z,
                               :ringinfo => "Field([Sqrt(5)])" ) )
    gens = AtlasGenerators( info )
    @test gens[1].base_ring === F

end

