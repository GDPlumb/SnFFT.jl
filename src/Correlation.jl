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
function sn_correlation(FFT1::Array{Array{Float64, 2}, 1}, FFT2::Array{Array{Float64, 2}, 1})
	l1 = length(FFT1)
	l2 = length(FFT2)
	L = min(l1, l2)
	Correlation = Array(Array{Float64, 2}, L)
	while L != 0
		Correlation[L] = FFT1[l1] * transpose(FFT2[l2])
		l1 -= 1
		l2 -= 1
		L -= 1
	end
	return Correlation
end

