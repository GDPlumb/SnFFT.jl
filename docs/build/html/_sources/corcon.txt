Convolution and Correlation
===========================
*SnFFT* implements two functions to help analyze the Fourier transform of functions over |Sn|.  
They are correlation and convolution.  
It is important to keep track of the degree of bandlimitedness of the input Fourier transforms because the result will have the higer degree of bandlimitedness.  

.. function:: sn_convoluation(FFT1, FFT2)

::

# Parameters:
#	FFT1::Array{Array{Float64, 2}, 1}
#	- the first Fourier transform
#	FFT2::Array{Array{Float64, 2}, 1}
#	- the second Fourier transform
# Return Values:
#	Convolution::Array{Array{Float64, 2}, 1}
#	- the convolution of FFT1 and FFT2
# Notes:
#	- FFT1 and FFT2 have to be Fourier transforms of functions over the same Sn
#	- However, the don't have to have the same degree of bandlimitedness

.. function:: sn_correlation(FFT1, FFT2)

::

# Parameters:
#	FFT1::Array{Array{Float64, 2}, 1}
#	- the first Fourier transform
#	FFT2::Array{Array{Float64, 2}, 1}
#	- the second Fourier transform
# Return Values:
#	Correlation::Array{Array{Float64, 2}, 1}
#	- the correlation of FFT1 and FFT2
# Notes:
#	- FFT1 and FFT2 have to be Fourier transforms of functions over the same Sn
#	- However, the don't have to have the same degree of bandlimitedness

.. |Sn| replace:: **S**\ :sub:`n` \

