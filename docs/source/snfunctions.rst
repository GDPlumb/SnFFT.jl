Functions over the Symmetric Group
==================================

| SnFFT represents a function over |Sn| as an array of Float64 values.  
| Because this representation doesn't explicitly store the permutation that corresponds to each value of the function, SnFFT has a set of standards that it uses to define the correspondence between indices of this array and the permutations.  
| These standards will be explained for each of the three types of fast Fourier transforms that SnFFT implements. 

Dense Functions
---------------

| A dense function over |Sn| will have a length of n! because the dense fast Fourier transform doesn't rely on any prior knowledge of the function. 
| The dense fast Fourier transform assumes that the value stored at index *i* is the value associated with the permutation that permutation_index() maps to *i*.
| See example3() for more details.    


.. function:: snf(N, PA, VA)

::

# Parameters:
#	N::Int
#	- the problem size
#	PA::Array{Array{Int, 1}, 1}
#	- PA[i] is a Permutation of N
#	VA::Array{Float64, 1}
#	- VA[i] is the Value associated with PA[i]
# Return Values:
#	SNF::Array{Float64, 1}
#	- SNF[i] is the value associated with the permutation that permutation_index() maps to i
#	- this is the format for the SNF parameter of sn_fft()
# Notes:
#	- any permutation of N not represented in PA will be assigned a value of zero

Bandlimited Functions
---------------------

| A bandlimited function over |Sn| that is invariant at **S**\ :sub:`n - k` \ will have n!/(n - k)! blocks of identical values of length (n - k)! when it is represented in the format that the dense fast Fourier transform uses.  
| This representation both wastes space and makes the calculation of the fast Fourier transform much slower. 
| Consequently, SnFFT uses a representation of such a function that stores one value from each block. 
| The bandlimited fast Fourier transform assumes the the value stored at index *i* is the value associated with all of the permutations that permutation_index() maps to (i - 1) * (n - k)!  + 1 to i * (n - k)!.  
| See example7() for more details.  

.. function:: snf_bl(N, K, PA, VA)

::

# Parameters:
#	N::Int
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
#	PA::Array{Array{Int, 1}, 1}
#	- PA[i] is a Permutation of N
#	VA::Array{Float64, 1}
#	- VA[i] is the Value associated with PA[i]
# Return Values:
#	SNF::Array{Float64, 1}
#	- SNF[i] is the value associated with all of the permutations that permutation_index() maps to any value in the range ((i - 1) * factorial(N - K) + 1):(i * factorial(N - K))
#	- this is the format for the SNF parameter of sn_fft_bl()
# Notes:
#	- any homogenous coset that doesn't have a representative permutation in PA will be assigned a value of zero


Sparse Functions
----------------

| A sparse function over |Sn| is represented by two components.  
| The first is a set of values and the second is a set of indices.  
| The sparse fast Fourier transform assumes the the value at index *i* is the value associated with the permutation that permutation_index() maps the index at index *i*.  
| See example5() for more details.  

.. function:: snf_sp(N, PA, VA)

::

# Parameters:
#	N::Int
#	- the problem size
#	PA::Array{Array{Int, 1}, 1}
#	- PA[i] is a Permutation of N
#	VA::Array{Float64, 1}
#	- VA[i] is the Value associated with PA[i]
# Return Values:
#	SNF::Array{Float64, 1}
#	- SNF[i] is the value associated with the permutation that permutation_index() maps to NZL[i] 
#	- this is the format for the SNF parameter of sn_fft_sp()
#	NZL::Array{Int, 1}
#	- NZL must in increasing order
#	- this is the format for the NZL parameter of sn_fft_sp()
# Notes:
#	- the values in VA should be non-zero

.. |Sn| replace:: **S**\ :sub:`n` \

