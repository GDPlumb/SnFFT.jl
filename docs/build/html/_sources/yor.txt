Young's Orthogonal Representations
==================================
*SnFFT* uses |YOR| (**YOR**) to calculate the |FFT| of a function over |Sn|.     
In addition to |YOR|, the |FFT| needs to know the structure that determines the decomposition of |YOR| (**PT**). 
Additionally, the bandlimited |FFT| needs some information about whether or not a component is a zero-frequency component (**ZFI**).  
To make computing multiple fast Fourier transforms more efficient, **YOR**, **PT**, and **ZFI** are computed before calling the |FFT|.  
They only need to be computed once because they don't depend on the specific values of the function over |Sn|.  

Dense and Sparse Functions
--------------------------
Before computing a dense or sparse |FFT|, construct the necessary information with: 

.. code-block:: julia

    julia> RA, PT = yor(N)

.. function:: yor(N)

::

# Parameters
#	N::Int
#	- the problem size
# Return Values
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1} (Young's Orthogonal Representations)
#	- YOR[n][p][k] is Young's Orthogonal Representation for the Adjacent Transposition (K, K + 1) for the pth Partition of n
#	PT::Array{Array{Array{Int, 1}, 1}, 1} (Partition Tree)
#	- for each value, i, in PT[n][p], P[n][p] decomposes into P[n-1][i]
#	- length(PT[1]) = 0

Bandlimited Functions
---------------------
Before computing a bandlimited |FFT|, construct the necessary information with:

.. code-block:: julia

    julia> RA, PT, ZFI = yor_bl(N, K)
    
.. function:: yor_bl(N,K)

::

# Parameters:
#	N::Int 
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
# Return Values
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1} (Young's Orthogonal Representations)
#	- YOR[n][p][k] is Young's Orthogonal Representation for the Adjacent Transposition (K, K + 1) for the pth Partition of n that is needed for the bandlimited functionality
#	- length(YOR[n]) = 0 for n = 1:(N - K - 1)
#	- if p < ZFI[n], length(YOR[n][p] = 1) and YOR[n][p][1,1] contains the dimension of the full Young's Orthogonal Representation
#	PT::Array{Array{Array{Int, 1}, 1}, 1} (Partition Tree)
#	- for each value, i, in PT[n][j], P[n][j] decomposes into P[n-1][i]
#	- length(PT[n]) = 0 for n = 1:(N - K)
#	- length(PT[n][p]) = 0 for p <= ZFI[n]
#	ZFI::Array{Int, 1} (Zero Frequency Information)
#	- ZFI[n] = k if, for p<=k, P[n][p] is a zero frequency partition

.. |Sn| replace:: **S**\ :sub:`n` \
.. |YOR| replace:: Young's Orthogonal Representations
.. |FFT| replace:: fast Fourier transform

