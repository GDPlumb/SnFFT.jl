Symmetric Group Functionality
=============================
*SnFFT* is designed to take the Fourier transform of functions over |Sn|.  
Although not strictly necessary to do this, having the basic functionality of the group can make testing and development much easier.  
These functions are in the file `Element.jl <https://github.com/GDPlumb/SnFFT.jl/blob/master/src/Element.jl>`_.  

Group Operations
----------------

.. function:: sn_multiply(P1,P2)

::

# Parameters:
#	P1::Array{Int, 1}
#	- the first permutation
#	P2::Array{Int, 1}
#	- the second permutation
# Return Values:
#	Prod::Array{Int, 1}
#	- the permutation that is P1 * P2
# Notes:
# 	- P1 and P2 must be permutations of the same size

.. function:: sn_inverse(P)

::

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
# Return Values:
#	Inv::Array{Int, 1}
#	- the permutation that is the inverse of P

Group Element Constructors
--------------------------

.. function:: sn_p(N)

::

# Parameters:
#	N::Int
#	- the size of the permutation
# Return Values:
#	P::Array{Int, 1}
#	- a random permutation of N

.. function:: sn_cc(N, LB, UB)

::

# Parameters:
#	N::Int
#	- the size of the permutation
#	LB::Int
#	- the first position that is reassigned
#	UB::Int
#	- the last position that is reassigned
# Return Values:
#	CC::Array{Int, 1}
#	- the permutation of N that is the contiguous cycle [[LB, UB]]
#	- this is the permutation that sends LB to LB + 1, LB + 1 to LB + 2, ... ,  UB - 1 to UB, and UB to LB
# Notes
#	- 1 <= LB <= UB <= N

.. function:: sn_cc(N) 

::

# Parameters:
#	N::Int
#	- the size of the permutation
# Return Values:
#	CC::Array{Int, 1}
#	- a random contiguous cycle of N 

.. function:: sn_at(N, K)

::

# Parameters:
#	N::Int
#	- the size of the permutation
#	K::Int
#	- the position that is being reassigned
# Return Values:
#	AT::Array{Int, 1}
#	- the permutation of N that is the adjacent transposition (K, K+1)
#	- this is the permutation that sends K to K + 1 and K + 1 to K
# Notes:
#	- 1 <= K < N

.. function:: sn_at(N)

::

# Parameters:
#	N::Int
#	- the size of the permutation
# Return Values:
#	AT::Array{Int, 1}
#	- a random adjacent transposition of N

.. function:: sn_t(N, I, J)

::

# Parameters:
#	N::Int
#	- the size of the permutation
#	I::Int
#	- the first postition that is being reassigned
#	J::Int
#	- the second position that is being reassigned
# Return Values:
#	Tr::Array{Int, 1}
#	- the permutation of N that is the transposition (I, J)
#	- this is the permutation that sends I to J and J to I
# Notes:
#	- 1 <= I <= N
#	- 1 <= J <= N

.. function:: sn_t(N)

::

# Parameters:
#	N::Int
#	- the size of the permutation
# Return Values:
#	Tr::Array{Int, 1}
#	- a random transposition of N

Factorizations and Related Operations on the Left-Coset Tree
------------------------------------------------------------

.. function:: permutation_ccf(P)

::

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
# Return Values:
#	CCF::Array{Int, 1}
#	- the Contiguous Cycle Factoriztion of P
#	- P = product for i = 1:(N - 1) of sn_cc(N, CCF[i], N + 1 - i)

.. function:: ccf_index(CCF)

::

# Parameters:
#	CCF::Array{Int, 1}
#	- a contiguous cycle factorization of some permutation
# Return Values:
#	Index::Int
#	- the unique index that the permutation corresponding to CCF maps to

.. function:: permutation_index(P)

::

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
# Return Values:
#	Index::Int
#	- the unique index that P maps to

.. function:: index_ccf(N, Index)

::

# Parameters:
#	N::Int
#	- the size of the permutation that maps to Index
#	Index::Int
#	- the index of some permutation of N
# Return Values:
#	CCF::Array{Int, 1}
#	- the contiguous cycle factorization that corresponds to the permutation that maps to Index

.. function:: ccf_permutation(N::Int, Index::Int)

::

# Parameters:
#	CCF::Array{Int, 1}
#	- a contiguous cycle factorization of some permutation
# Return Values:
#	P::Array{Int, 1}
#	- the permutation that corresponds to CCF

.. function:: index_permutation(N, Index)

::

# Parameters:
#	N::Int
#	- the size of the permutation that maps to Index
#	Index::Int
#	- the index of some permutation of N
# Return Values:
#	P::Array{Int, 1}
#	- the permutation of N that maps to Index

.. function:: permutation_atf(P)

::

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
# Return Values
#	ATF::Array{Int, 1}
#	- the adjacent transposition factorization of P
#	- P = product for i = 1:length(ATF) of sn_at(N, ATF[i])

Young's Orthogonal Representation of a Permutation
----------------------------------------------------

.. function:: yor_permutation(P, YORnp)

::

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
#	YORnp::Array{SparseMatrixCSC, 1}
#	- YORnp[i] is Young's Orthogonal Representation for the adjacent transposition (i, i + 1) corresponding to the pth partition of n
# Return Values:
# 	RM::Array{Float64, 2}
#	- Young's Orthogonal Representation of P corresponding to the pth partition of n

.. |Sn| replace:: **S**\ :sub:`n` \

