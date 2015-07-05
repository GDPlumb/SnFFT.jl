# Let X be a Yamanouchi Symbol corresponding to a Partition of N whose length is R
#	X::Array{Int, 1}
#	length(X) = N
#	maximum(X) = R
#	X[i] = r means that i appears in row r of the Standard Tableau represented by X

# Parameters:
#	N::Int
#	- the problem size
#	P::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 from partitions()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 from partition_tree()
# Return Values
#	YS::Array{Array{Array{Array{Int, 1}, 1}, 1}, 1} (Yamanouchi Symbols)
#	- YS[n][p][i] is the ith Yamanouchi Symbol of the pth Partition of n
function ys_symbols(N::Int, P::Array{Array{Array{Int, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1})
	YS = Array(Array{Array{Array{Int, 1}, 1}, 1}, N)
	YSn = Array(Array{Array{Int, 1}, 1}, 1)
	YSnp = Array(Array{Int, 1}, 1)
	YSnp[1] = [1]
	YSn[1] = YSnp
	YS[1] = YSn
	for n = 2:N
		Pn = P[n]
		PTn = PT[n]
		Pn_L = length(Pn)
		YSn = Array(Array{Array{Int, 1}, 1}, Pn_L)
		for p = 1:Pn_L
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

# Parameters:
#	N::Int
#	- the size of the Partition that the elements of YSymbols correspond to
#	R::Int  
#	- the length of the Partition that the elements of YSymbols correspond to
#	YSymbols::Array{Array{Int, 1}, 1}
#	- the array of Yamanouchi Symbols corresponding to a particular Partition
# Return Values
#	DA::Array{Int, 2} (Distance Array)
#	- output1 from ys_distance()
#	IA::Array{Int, 2} (Index Array)
#	- output1 from ys_indices()
function ys_information(N::Int, R::Int, YSymbols::Array{Array{Int, 1}, 1})
	CA = ys_content(N, R, YSymbols)
	DA = ys_distance(CA)
	IA = ys_indices(N, YSymbols, DA)
	return DA, IA
end

#Parameters:
#	N::Int
#	- the size of the Partition that the elements of YSymbols correspond to
#	R::Int  
#	- the length of the Partition that the elements of YSymbols correspond to
#	YSymbols::Array{Array{Int, 1}, 1}
#	- the array of Yamanouchi Symbols corresponding to a particular Partition
#Return Value
#	CA::Array{Int, 2} (Content Array)
#	- CA[i, n] is the Content of n for the ith Standard Tableau
function ys_content(N::Int, R::Int, YSymbols::Array{Array{Int, 1}, 1})	
	CA = Array(Int, length(YSymbols), N)
	for i = 1:length(YSymbols)
		YSymbol = YSymbols[i]
		PF = zeros(Int, R)
		for n = 1:N
			row = YSymbol[n]
			PF[row] += 1
			CA[i,n] = PF[row] - row
		end
	end
	return CA
end

# Parameters
#	CA::Array{Int, 2}
#	- output1 from ys_content()
# Return Values
#	DA::Array{Int, 2} (Distance Array)
#	- DA[i, k] is the Distance of the ith Standard Tableau for the Adjacent Transposition (k, k + 1)
function ys_distance(CA::Array{Int, 2})
	L, N = size(CA)
	DA = Array(Int, L, N - 1)
	for i = 1:L
		for n = 1:(N - 1)
			DA[i, n] = CA[i, n + 1] - CA[i, n]
		end
	end
	return DA
end

#Parameters
#	N::Int
#	- the size of the Partition that the elements of YSymbols correspond to
#	YSymbols::Array{Array{Int, 1}, 1}
#	- the array of Yamanouchi Symbols corresponding to a particular Partition
#	DA::Array{Int, 2}
#	- output1 from ys_distance()
#Return Values
#	IA::Array{Int, 2} (Index Array)
#	- IA[i, k] contains the index of the Standard Tableau that is the result of applying the Adjacent Transposition (k, k + 1) to the Standard Tableau at index i in YSymbols
#	- IA[i, k] = 0 when applying the Adjacent Transposition (K, K + 1) to the Standard Tableau at index i doesn't result in a Standard Tableau 
#Notes 
#	- there are block patterns that appear in this that could be used to speed up this search
#	- there may be ways to avoid searching all together
function ys_indices(N::Int, YSymbols::Array{Array{Int, 1}, 1}, DA::Array{Int, 2})
	L = length(YSymbols)
	IA = zeros(Int, L, N - 1)
	PF = zeros(Int, L, N - 1)
	for i = 1:L
		YSymbol = YSymbols[i]
		k = 1
		ti = i + 1
		while k < N 
			D = DA[i, k]
			if D != 1 && D != -1 && PF[i, k] == 0
				while ti <= L
					if istranspose(YSymbol, YSymbols[ti], k) == 1
						IA[i, k] = ti
						IA[ti, k] = i
						PF[i, k] = 1
						PF[ti, k] = 1
						break
					end
					ti += 1
				end
			end
			k += 1
		end
	end
	return IA
end

# Parameters
#	YSymbol1::Array{Int, 1}
#	- the starting Yamanouchi Symbol
#	YSmobl2::Array{Int, 2}
#	- the potential resulting Yamanouchi Symbol
#	K::Int
#	- represents the Adjacent Transposition (K, K + 1)
# Return Values
#	- 1, if YSymbol2 is the result of applying the Adjacent Transposition (K, K + 1) to YSymbol1
#	- 0, otherwise
function istranspose(YSymbol1::Array{Int, 1}, YSymbol2::Array{Int, 1}, K::Int)
	for i = 1:(K - 1)
		if YSymbol1[i] != YSymbol2[i]
			return 0
		end
	end
	if YSymbol1[K] != YSymbol2[K + 1] || YSymbol1[K + 1] != YSymbol2[K]
		return 0
	end
	for i = (K + 2):length(YSymbol1)
		if YSymbol1[i] != YSymbol2[i]
			return 0
		end
	end
	return 1
end

