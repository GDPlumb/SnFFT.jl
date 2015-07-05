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
function snf(N::Int, PA::Array{Array{Int, 1}, 1}, VA::Array{Float64, 1})
	SNF = zeros(Float64, factorial(N))
	for i = 1:length(PA)
		Permutation = PA[i]
		Index = permutation_index(Permutation)
		SNF[Index] = VA[i]
	end
	return SNF
end

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
function snf_bl(N::Int, K::Int, PA::Array{Array{Int, 1}, 1}, VA::Array{Float64, 1})
	N_F = factorial(N)
	BS = factorial(N - K)
	SNF = zeros(Float64, int(N_F / BS))
	for i = 1:length(PA)
		Permutation = PA[i]
		Index = permutation_index(P)
		Index = ceil(Index / BS)
		SNF[Index] = VA[i]
	end
	return SNF
end

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
function snf_sp(N::Int, PA::Array{Array{Int, 1}, 1}, VA::Array{Float64, 1})
	L = length(PA)	
	SNF = Array(Float64, L)
	NZL = Array(Int, L)
	SNF[1] = VA[1]
	NZL[1] = permutation_index(PA[1])
	for i = 2:L
		val = VA[i]
		index = permutation_index(PA[i])
		j = 1		
		while j <= (i - 1)
			if index < NZL[j]
				break
			end
			j += 1
		end
		for k = i:-1:(j + 1)
			SNF[k] = SNF[k - 1]
			NZL[k] = NZL[k - 1]
		end
		SNF[j] = val
		NZL[j] = index
	end
	return SNF, NZL
end
