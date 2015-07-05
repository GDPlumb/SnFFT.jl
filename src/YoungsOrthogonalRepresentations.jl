# Parameters
#	N::Int
#	- the problem size
# Return Values
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1} (Young's Orthogonal Representations)
#	- YOR[n][p][k] is Young's Orthogonal Representation for the Adjacent Transposition (K, K + 1) for the pth Partition of n
#	PT::Array{Array{Array{Int, 1}, 1}, 1} (Partition Tree)
#	- output1 from partition_tree()
function yor(N::Int)
	P, WI = partitions(N)
	PT = partition_tree(N, P, WI)
	YS = ys_symbols(N, P, PT)
	YOR = Array(Array{Array{SparseMatrixCSC, 1}, 1}, N)
	for n = 1:N
		YSn = YS[n]
		Pn = P[n]
		Pn_L = length(Pn)
		YORn = Array(Array{SparseMatrixCSC, 1}, Pn_L)
		for p = 1:Pn_L
			YSnp = YSn[p]
			Pnp = Pn[p]
			R = length(Pnp)
			YORn[p] = yor_p(n, R, YSnp)
		end
		YOR[n] = YORn
	end
	return YOR, PT
end

# Parameters:
#	N::Int
#	- the size of P
#	R::Int
#	- the length of P
#	P::Array{Int, 1}
#	- P is a Partition of N
# Return Values:
#	YORp::Array{SparseMatrixCSC, 1}
#	- YORp[k] is the Young's Orthogonal Representation of P for the Adjacent Transposistion (k, k + 1)
function yor_p(N::Int, R::Int, YSymbols::Array{Array{Int, 1}, 1})
	YORp = Array(SparseMatrixCSC, N - 1)
	L = length(YSymbols)
	DA, IA = ys_information(N, R, YSymbols)
	for k = 1:(N - 1)
		YORp[k] = yor_pk(DA, IA, L, k)
	end
	return YORp
end

#Parameters
#	DA::Array{Int, 2}  
#	- output1 from ys_information()
#	IA::Array{Int, 2}
#	- outpu2 from ys_information()
#	L::Int
#	- the number of Yamanouchi Symbols for the current Partition
#	K::Int
#	- represents the Adjacent Transposition (K, K + 1)
#Return Values
#	YORpk::SparseMatrixCSC 
#	- Young's Orthongal Representation of the Adjacent Transposition (K, K + 1) for the Partition p
function yor_pk(DA::Array{Int, 2}, IA::Array{Int, 2}, L::Int, k::Int)
	n = L
	for i = 1:L
		if IA[i,k] != 0
			n += 1
		end
	end
	colptr = Array(Int, L + 1)
	colptr[1] = 1
	rowval = Array(Int, n)
	nzval = Array(Float64, n)
	cp = 1
	for i = 1:L
		D = 1/DA[i, k]
		I = IA[i,k]
		if I == 0
			colptr[i + 1] = colptr[i] + 1
			rowval[cp] = i
			nzval[cp] = D
			cp += 1
		else
			colptr[i + 1] = colptr[i] + 2
			if I > i
				rowval[cp] = i
				nzval[cp] = D
				cp += 1
				rowval[cp] = I
				nzval[cp] = sqrt(1 - D * D)
				cp += 1
			else
				rowval[cp] = I
				nzval[cp] = sqrt(1 - D * D)
				cp += 1
				rowval[cp] = i
				nzval[cp] = D
				cp += 1
			end
		end
	end
	YORpk = SparseMatrixCSC(L, L, colptr, rowval, nzval) 
	return YORpk
end

