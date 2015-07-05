# Let X be a Permutation of N
#	X::Array{Int, 1}
#	length(X) == N
#	X[i] = j indicates that the item in position i is sent to position j

###
# Group Operations
###

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
function sn_multiply(P1::Array{Int, 1}, P2::Array{Int, 1})
	N = length(P1)   
	Prod = Array(Int, N)
	for i = 1:N
		Prod[i] = P1[P2[i]]
	end
	return Prod
end

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
# Return Values:
#	Inv::Array{Int, 1}
#	- the permutation that is the inverse of P	
function sn_inverse(P::Array{Int, 1})
	N = length(P)
	Inv = Array(Int, N)
	for i = 1:N
		Inv[P[i]] = i
	end
	return Inv
end

###
# Permutation Constructors
###

# Parameters:
#	N::Int
#   - the size of the permutation
# Return Values:
#	P::Array{Int, 1}
#	- a random permutation of N
function sn_p(N::Int)
	index = rand(1:factorial(N))
	P = index_permutation(N, index)
	return P
end

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
function sn_cc(N::Int, LB::Int, UB::Int)
	CC = Array(Int, N)
	for i = 1:(LB - 1)
		CC[i] = i
	end
	for i = LB:(UB - 1)
		CC[i] = i + 1
	end
	CC[UB] = LB
	for i = (UB + 1):N
		CC[i] = i
	end
	return CC
end

# Parameters:
#	N::Int
#	- the size of the permutation
# Return Values:
#	CC::Array{Int, 1}
#	- a random contiguous cycle of N
function sn_cc(N::Int)
	lb = rand(1:N)
	ub = rand(lb:N)
	CC = sn_cc(N, lb, ub)
	return CC
end

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
function sn_at(N::Int, K::Int)
	AT = Array(Int, N)
	for i = 1:(K - 1)
		AT[i] = i
	end
	AT[K] = K + 1
	AT[K + 1] = K
	for i = (K + 2):N
		AT[i] = i
	end
	return AT
end

# Parameters:
#	N::Int
#	- the size of the permutation
# Return Values:
#	AT::Array{Int, 1}
#	- a random adjacent transposition of N
function sn_at(N::Int)
	k = rand(1:(N - 1))
	AT = sn_at(N, k)
	return AT
end

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
function sn_t(N::Int, I::Int, J::Int)
	Tr = Array(Int, N)
	for i = 1:N
		Tr[i] = i
	end
	Tr[I] = J
	Tr[J] = I
	return Tr
end

# Parameters:
#	N::Int
#	- the size of the permutation
# Return Values:
#	Tr::Array{Int, 1}
#	- a random transposition of N
function sn_t(N::Int)
	i = rand(1:N)
	j = rand(1:N)
	Tr = sn_t(N, i, j)
	return Tr
end

###
# Factorizations (and related operations) on the Left Coset Tree
###

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
# Return Values:
#	CCF::Array{Int, 1}
#	- the Contiguous Cycle Factoriztion of P
#	- P = product for i = 1:(N - 1) of sn_cc(N, CCF[i], N + 1 - i)
function permutation_ccf(P::Array{Int, 1})
	N = length(P)
	CCF = Array(Int, N - 1)
	i = 1
	for j = N:-1:2
		CCF[i] = P[j]
		cc = sn_cc(N, P[j], j)
		cc_inv = sn_inverse(cc)
		P = sn_multiply(cc_inv, P)
		i += 1
	end
	return CCF
end

# Parameters:
#	CCF::Array{Int, 1}
#	- a contiguous cycle factorization of some permutation
# Return Values:
#	Index::Int
#	- the unique index that the permutation corresponding to CCF maps to
function ccf_index(CCF::Array{Int, 1})
	N = length(CCF) + 1
	Index = 1
	for i = 1:(N - 1)
		N -= 1  	    
		if CCF[i] != 1
			Index += (CCF[i] - 1) * factorial(N)
		end
	end
	return Index
end

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
# Return Values:
#	Index::Int
#	- the unique index that P maps to
function permutation_index(P::Array{Int, 1})
	ccf = permutation_ccf(P)
	index = ccf_index(ccf)
	return index
end

# Parameters:
#	N::Int
#	- the size of the permutation that maps to Index
#	Index::Int
#	- the index of some permutation of N
# Return Values:
#	CCF::Array{Int, 1}
#	- the contiguous cycle factorization that corresponds to the permutation that maps to Index
function index_ccf(N::Int, Index::Int)
	CCF = Array(Int, N - 1)
	Index -= 1
	for i = 1:(N - 1)
		q = floor(Index / factorial(N - i))
		Index -= q * factorial(N - i)
		CCF[i] = q + 1
	end
	return CCF
end

# Parameters:
#	CCF::Array{Int, 1}
#	- a contiguous cycle factorization of some permutation
# Return Values:
#	P::Array{Int, 1}
#	- the permutation that corresponds to CCF
function ccf_permutation(CCF::Array{Int, 1})
	N = length(CCF) + 1
	P = Array(Int, N)
	for i = 1:N
		P[i] = i
	end
	for i = 1:(N - 1)
		cc = sn_cc(N, CCF[i], N + 1 - i)
		P = sn_multiply(P, cc)
	end
	return P
end

# Parameters:
#	N::Int
#	- the size of the permutation that maps to Index
#	Index::Int
#	- the index of some permutation of N
# Return Values:
#	P::Array{Int, 1}
#	- the permutation of N that maps to Index
function index_permutation(N::Int, Index::Int)
	ccf = index_ccf(N, Index)
	permutation = ccf_permutation(ccf)
	return permutation
end

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
# Return Values
#	ATF::Array{Int, 1}
#	- the adjacent transposition factorization of P
#	- P = product for i = 1:length(ATF) of sn_at(N, ATF[i])
function permutation_atf(P::Array{Int, 1})
	N = length(P)
	CCF = permutation_ccf(P)
	dATF = Array(Array{Int, 1}, length(CCF))
	L = 0
	for i = 1:length(CCF)
		P = CCF[i]
		D = N - P
		P -= 1
		L += D
		N -= 1
		dATFs = Array(Int, D)
		for j = 1:D
			dATFs[j] = P + j
		end
		dATF[i] = dATFs
	end
	ATF = Array(Int, L)
	index = 1
	for i = 1:length(CCF)
		for j = 1:length(dATF[i])
			ATF[index] = dATF[i][j]
			index += 1
		end
	end
	return ATF
end

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
#	YORnp::Array{SparseMatrixCSC, 1}
#	- YORnp[i] is Young's Orthogonal Representation for the adjacent transposition (i, i + 1) corresponding to the pth partition of n
# Return Values:
# 	RM::Array{Float64, 2}
#	- Young's Orthogonal Representation of P corresponding to the pth partition of n
function yor_permutation(P::Array{Int, 1}, YORnp::Array{SparseMatrixCSC, 1})
	ATF = permutation_atf(P)
	if length(ATF) == 0
		Dim = size(YORnp[1], 1)
		RM = eye(Dim)
		return RM
	else
		RM = copy(full(YORnp[ATF[1]]))
		for i = 2:length(ATF)
			RM = RM * full(YORnp[ATF[i]])
		end
		return RM
	end
end
