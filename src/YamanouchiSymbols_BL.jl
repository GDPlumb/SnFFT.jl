# Parameters:
#	N::Int 
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
#	P::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 from partitions_bl()
#	ZFI::Array{Int, 1}
#	- output2 from partitions_bl()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 from partition_tree_bl()
# Return Values:
#	YS::Array{Array{Array{Array{Int, 1}, 1}, 1}, 1} (Yamanouchi Symbols)
#	- YS[n][p][i] is the ith Yamanouchi Symbol of the pth Partition of n is required for the bandlimited functionality
#	- length(YS[n]) = 0 for n = 1:(N - K - 1)
function ys_symbols_bl(N::Int, K::Int,  P::Array{Array{Array{Int, 1}, 1}, 1}, ZFI::Array{Int, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1})
	YS = Array(Array{Array{Array{Int, 1}, 1}, 1}, N)
	for n = 1:(N - K - 1)
		YS[n] = Array(Array{Array{Int, 1}, 1}, 0)
	end
	YSn = Array(Array{Array{Int, 1}, 1}, length(P[N - K]))
	for p = 1:length(P[N - K])
		Pnp = P[N - K][p]
		YSn[p] = ys_partition(N - K, Pnp)
	end
	YS[N - K] = YSn
	for n = (N - K + 1):N
		Pn = P[n]
		PTn = PT[n]
		Pn_L = length(Pn)
		YSn = Array(Array{Array{Int, 1}, 1}, Pn_L)
		for p = 1:ZFI[n]
			Pnp = Pn[p]
			YSn[p] = ys_partition(n, Pnp)
		end
		for p = (ZFI[n] + 1):Pn_L
			Pnp = Pn[p]
			Pnp_L = length(Pnp)
			PTnp = PTn[p]
			PTnp_L = length(PTnp)
			c = 1
			row = 1
			YSnp = Array(Array{Int, 1}, degree(n, Pnp))
			for d = 1:PTnp_L
				while row < Pnp_L && Pnp[row] <= Pnp[row + 1]
					row += 1
				end
				YSd = YS[n - 1][PTnp[d]]
				YSd_L = length(YSd)
				for i = 1:YSd_L
					YSnpc = Array(Int, n)
					YSdi = YSd[i]
					for j = 1:(n - 1)
						YSnpc[j] = YSdi[j]
					end
					YSnpc[n] = row
					YSnp[c] = YSnpc 
					c += 1
				end
				row += 1
			end
			YSn[p] = YSnp
		end
		YS[n] = YSn
	end
	return YS
end

# Parameters
#	N::Int
#	- N is the size of P
#	P::Array{Int, 1}
#	- P is a Partition of N
# Return Values
#	YS_P::Array{Array{Int, 1}, 1} (Yamanouchi Sybol for the Partition P)
#	- YS_P[i] is the ith Yamanouchi Symbol of the Partition P
function ys_partition(N::Int, P::Array{Int, 1})
	if length(P) == 1
		YS_P = Array(Array{Int, 1}, 1)
		YS_P[1] = ones(Int, N)
		return YS_P
	end	
	NST = degree(N, P)	
	YSymbol = Array(Int, N)
	PA = copy(P)
	YS_P = Array(Array{Int, 1}, NST)
	C = Counter(NST)
	ys_fillin(N, YSymbol, PA, YS_P, C)
	return YS_P	
end

# Parameters:
#	n::Int
#	- the next element to be placed into YSymbol
#	YSymbol::Array{Int, 1}
#	- the partially completed Yamanouchi Symbol that is being filled in
#	PA::Array{Int, 1}
#	- PA[i] is the number of Positions Available in row i of the Standard Tableau corresponding to YSymbol
#	YS_P::Array{Array{Int, 1}, 1}
#	- the array that stores the Yamanouchi symbols as they are completed
#	C::Counter
#	- A wrapper for an int that allows the recursive calls to know the next empty index in YS_P
# Return Values:
#	NONE
#	- once the recursion is finished, YS_P is full 
# Notes:
#	- This shouldn't be called on its own, use its wrapper function ys_partition()
function ys_fillin(n::Int, YSymbol::Array{Int, 1}, PA::Array{Int, 1}, YS_P::Array{Array{Int, 1}, 1}, C::Counter)
	PA_L = length(PA)	
	if PA[PA_L] > 0
		YSymbol_new = copy(YSymbol)
		YSymbol_new[n] = PA_L
		PA_new = copy(PA)
		PA_new[PA_L] -= 1
		ys_fillin(n - 1, YSymbol_new, PA_new, YS_P, C)
	end
	for i = (PA_L - 1):(-1):2 
		if PA[i] > PA[i + 1]
			YSymbol_new = copy(YSymbol)
			YSymbol_new[n] = i
			PA_new = copy(PA)
			PA_new[i] -= 1
			ys_fillin(n - 1, YSymbol_new, PA_new, YS_P, C)
		end
	end
	if PA[1] > PA[2]
		if PA[1] == 1
			YSymbol_new = copy(YSymbol)
			YSymbol_new[1] = 1
			YS_P[C.N] = YSymbol_new
			C.N -= 1
		else
 			YSymbol_new = copy(YSymbol)
			YSymbol_new[n] = 1
			PA_new = copy(PA)
			PA_new[1] -= 1
			ys_fillin(n - 1, YSymbol_new, PA_new, YS_P, C)
		end
	end
end

