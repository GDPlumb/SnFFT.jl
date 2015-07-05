Miscellaneous Functions
=======================

Partition Construction
----------------------
Although perhaps not computationally useful, *SnFFT* does export the function to construct the set of paritions of 1:N.  
Generally, this is useful for giving the output of other code easier to interpret.  

.. function:: partitions(N)

::

# Parameters:
#	N::Int 
#	- the problem size
# Return Values:
#	P::Array{Array{Array{Int, 1}, 1}, 1} (Partitions)
#	- P[n][p] contains the pth Partition of n
#	WI::Array{Int, 2} (Width Information)
#	- WI[n, w] contains the number of Paritions of n whose first element is less than or equal to w

Preference Matrices
-------------------
Constructs the preference matrix for a permutation.  

.. function:: preferencematrix(P)

::

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
# Return Values:
#	Q::Array{Int, 2}
#	- the preference matrix for P
#	- Q[i,j] = 1 if and only if j precedes i in P

Kendall Tau Distance
--------------------
Computes the Kendall Tau distance between two permutations using the permutation's preference matrices.  

.. function:: kendalldistance(Q1, Q2)

::

# Parameters:
#	Q1::Array{Int, 2}
#	- the preference matrix for the first permutation
#	Q2::Array{Int, 2}
#	- the preference matrix for the second permutation
# Return Values:
#	D::Int
#	- the Kendall Tau Distance between two permutations
#	- the the number of pairs (i, j) such that: P1[i] < P1[j] and P2[i] > P2[j]

Mallow's Distribution
---------------------
Constructs a Mallow's Distribution centered around a specified permutation with a given spreading factor.  

.. function:: mallowsdistribution(P, Gamma)

::

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
#	Gamma::Float64
#	- the spreading factor
# Return Values:
#	MD::Array{Float64, 1}
#	- the mallows distribution with spreading factor Gamma centered at P

Printing Methods
----------------
.. function:: permutation_string(Permutation)

::

# Parameters:
#	Permutation::Array{Int, 1}
#	- a permutation
# Return Values:
#	ST::String
#	- the string representation of permutation 

.. function:: partition_strign(Partition)

::

# Parameters:
#	Partition::Array{Int, 1}
#	- a partition
# Return Values:
#	ST::String
#	- the string representation of partition
