Fast Fourier Transforms
=======================
*SnFFT* supports three types of fast Fourier transforms.
The dense |FFT| can use any type of function over |Sn|.  
The bandlimited |FFT| can use only bandlimited functions over |Sn|, but benefits greatly from the restriction.
The sparse |FFT| can run on any type of function over |Sn|, but becomes faster as the function becomes increasingly sparse. 

Dense Fast Fourier Transform
----------------------------
.. function:: sn_fft(N, SNF, YOR, PT)

::

# Parameters:
#	N::Int
#	- the problem size
#	SNF::Array{Foat64, 1}
#	- SNF[i] is the value associated with the Permutation that permutation_index() maps to i
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor()
# Return Values:
#	FFT::Array{Float64, 2}
#	- FFT is the Fast Fourier Transform of SNF

* See snf() for a detailed explanation of the **SNF** parameter
* **YOR** and **PT** are the outputs of yor(**N**)
* See the code for example4() for an example of the complete process to compute a dense |FFT|

Bandlimited Fast Fourier Transform
----------------------------------
.. function:: sn_fft_bl(N, K, SNF, YOR, PT, ZFI)

::

# Parameters:
#	N::Int
#	- the problem size N
#	K::Int
#	- the problem is homogenous at N-K
#	SNF::Array{Foat64, 1}
#	- SNF[i] is the value assigned to the ith homogenous subgroup of size N-K
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor_bl()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor_bl()
#	ZFI::Array{Int, 1}
#	- output3 from yor_bl()
# Return Values:
#	FFT::Array{Float64, 2}
#	- FFT is the Fast Fourier Transform of SNF

* See snf_bl() for a detailed explanation of the **SNF** parameter
* **YOR**, **PT**, and **ZFI** are the outputs of yor_bl(**N**, **K**)
* See the code for example7() for an example of the complete process to compute a bandlimited |FFT|

Sparse Fast Fourier Transform
-----------------------------
.. function:: sn_fft_sp(N, SNF, NZL, YOR, PT)

::

# Parameters:
#	N::Int
#	- the problem size is N
#	SNF::Array{Foat64, 1}
#	- SNF[i] is the value associated with the Permutation that permutation_index() maps to NZL[i]
#	NZL::Array{Int, 1}
#	- NZL[i] the set of NonZeroLocations for the sparse function over Sn
#	- NZL must be in increasing order
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor()
# Return Values:
#	FFT::Array{Float64, 2}
#	- FFT is the Fast Fourier Transform of SNF

* See snf_sp() for a detailed explanation of the **SNF** and **NZL** parameters
* **YOR** and **PT** are the outputs of yor(**N**)
* See the code for example5() for an example of the complete process to compute a dense |FFT|

.. |Sn| replace:: **S**\ :sub:`n` \
.. |FFT| replace:: fast Fourier transform

