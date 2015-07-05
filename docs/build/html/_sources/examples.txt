Examples
========
*SnFFT* comes with eight example functions that demonstrate some of the key properties of Young's Orthogonal Representations and the Fourier Transform of functions over |Sn|.  
Furthermore, they demonstrate the syntax used to call most of the functionality of the package. 
The examples are in the file `Examples.jl <https://github.com/GDPlumb/SnFFT.jl/blob/master/src/Examples.jl>`_.  

General Notes:
--------------
* Most of the example functions have two methods.  The first method has no parameters and will run the example with default values. The second has a full set of parameters and will be explained in each function's description. 
* A parameter name will appear in **bold** when it being referenced in the example's description.  
* *SnFFT* represents a `partition <http://en.wikipedia.org/wiki/Partition_(number_theory)>`_ in the following way:

::

# Let X be a Partition of N
#	X::Array{Int, 1}
#	X[i] > 0 for all i
#	sum(X) == N
#	X[i] >= X[j] when i < j

* *SnFFT* represents a `permutation <http://en.wikipedia.org/wiki/Permutation>`_ in the following way:

::

# Let X be a Permutation of N
#	X::Array{Int, 1}
#	length(X) == N
#	X[i] = j indicates that the item in position i is sent to position j

Example Function Explanations
-----------------------------

.. function:: example1(N, partition, permutation) 

::

# Parameters:
#	N::Int
#	- the problem size
#	partition::Array{Int, 1}
#	- a partition of N
#	permutation::Array{Int, 1}
#	- a permutation of N

This example finds the `Standard Tableau <http://en.wikipedia.org/wiki/Young_tableau>`_ corresponding to **partition**.  
It then calculates Young's Orthogonal Representation of **permutation** corresponding to **partition**.  

.. function:: example2(N, partition, p1, p2)

::

# Parameters:
#	N::Int
#	- the problem size
#	partition::Array{Int, 1}
#	- a partition of N
#	p1::Array{Int, 1}
#	- the first permutation of N
#	p2::Array{Int, 1}
#	- the second permutation of N

Let YOR1 and YOR2 be Young's Orthogonal Representations of **p1** and **p2** corresponding to **partition**.  
Let YORm be Young's Orthogonal Representations of **p1** x **p2** corresponding to **partition**. 
This example demonstrates that YORm = YOR1 x YOR2.  

.. function:: example3()

This example demonstrates the order of the permutations used by *SnFFT* to represent a function over |Sn|.  
It has no parameterized version.  

.. function:: example4(N)

::

# Parameters:
#	N::Int
#	- the problem size

This example constructs a random function over |Sn|.  
It then demonstrates how to calculate the (dense) fast Fourier transform and the inverse fast Fourier transform.
Finally, it shows that this process recovers the original function accurately.  

.. function:: example5(N, SC)

::

# Parameters:
#	N::Int
#	- the problem size
#	SC::Float64
#	- the portion of the function that is zero-valued

This example constructs a random sparse function over |Sn| and converts it to the format used to compute the spare fast Fourier transform. 
Next, it demonstrates how to compute the sparse fast Fourier transform.
Finally, it shows that this produces the same result as the dense fast Fourier transform.  

.. function:: example6(N, permutation)

::

# Parameters:
#	N::Int
#	- the problem size
#	Permutation::Array{Int, 1}
#	- a permutation of N

This example constructs a delta function over |Sn| that is centered on **permutation**. 
It then calculates the sparse fast Fourier transform of this function.  
Finally, it demonstrates that this produces the same results as computing Young's Orthogonal Representation for each partition of **N** corresponding to **permutation**.  

.. function:: example7(N, K)

::

# Parameters:
#	N::Int 
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K

This example constructs a random bandlimited function over |Sn|. 
To save space, the bandlimited function is compressed.  
It then constructs the equivalent non-compressed version.  
Next, it demonstrates how to compute the bandlimited fast Fourier transform using the compressed function and take the bandlimited inverse fast Fourier transform. 
Finally, it shows that the dense and bandlimited fast Fourier transforms produce the same result and that the inverse bandlimited fast Fourier transform recovers the original bandlimited function over |Sn|. 

.. function:: example8(N, M)

::

# Parameters:
#	N::Int 
#	- the problem size
#	M::Int
#	- the number of the top frequencies of the Fourier transform to use

This function constructs a random function over |Sn| and then demonstrates how to compute the paritial inverse fast Fourier transform.
It prints both the original and recovered function.  

.. |Sn| replace:: **S**\ :sub:`n` \

