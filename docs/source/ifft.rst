Inverse Fast Fourier Transform
==============================
*SnFFT* support three types of inverse fast Fourier transforms.  
The dense |iFFT| can take the inverse of the output of either the sparse or dense |FFT|.  
The bandlimited |iFFT| can only take the inverse of the output of the bandlimited |FFT|, but benefits greatly from the restriction.  
The partial |iFFT| can take the inverse of the output of either the sparse or dense |FFT|.  

Dense Inverse Fast Fourier Transform
------------------------------------
.. function:: sn_ifft(N, FFT, YOR, PT)

::

# Parameters:
#	N::Int
#	- the problem size
#	FFT::Array{Array{Float64, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft() or sn_fft_sp()
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor()
# Return Values:
#	SNF::Array{Float64, 1}
#	- the function over Sn that corresponds to FFT

* **FFT** is the output of sn_ftt() or sn_fft_sp()
* **YOR** and **PT** are the outputs of yor(**N**)
* See the code for example4() for an example of the complete process to compute a dense |iFFT|

Bandlimited Inverse Fast Fourier Transform
------------------------------------------
.. function:: sn_ifft_bl(N, K, FFT, YOR, PT, ZFI)

::

# Parameters:
#	N::Int
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
#	FFT::Array{Array{Float64, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft_bl()
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor_bl()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor_bl()
#	ZFI::Array{Int, 1}
#	- output3 from yor_bl()
# Return Values:
#	SNF::Array{Float64, 1}
#	- the bandlimited function over Sn that corresponds to FFT

* **FFT** is the output of sn_fft_bl()
* **YOR**, **PT**, and **ZFI** are the outputs of yor_bl(**N**, **K**)
* See the code for example7() for an example of the complete process to compute a bandlimited |iFFT|

Partial Inverse Fast Fourier Transform
---------------------------------------
.. function:: sn_ifft_p(N, M, FFT, YOR, PT)

::

# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the number of the top components of FFT to use
#	FFT::Array{Array{Float64, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft() or sn_fft_sp()
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor()
# Return Values:
#	SNF::Array{Float64, 1}
#	- the function over Sn that corresponds to FFT that has been smoothed to only use the top M componenets

* **FFT** is the output of sn_ftt() or sn_fft_sp()
* **YOR** and **PT** are the outputs of yor(**N**)
* See the code for example8() for an example of the complete process to compute a partial |iFFT|

.. |FFT| replace:: fast Fourier transform
.. |iFFT| replace:: inverse fast Fourier transform

