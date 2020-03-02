
# uniform access to the number of elements in a finite field

# residue class ring, no prim. root and no `order` method
fieldorder( F::Nemo.NmodRing ) = Int( F.n )

# prime field, no prim. root
fieldorder( F::Nemo.GaloisField ) = Int( Nemo.order( F ) )

# field defined by a Conway polynomial
fieldorder( F::Nemo.FqNmodFiniteField ) = Int( Nemo.order( F ) )

