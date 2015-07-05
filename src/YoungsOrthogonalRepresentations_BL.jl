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
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 of partition_tree_bl()
#	ZFI::Array{Int, 1}
#	- output2 for paritions_bl()
function yor_bl(N::Int, K::Int)
	P, ZFI, WI = partitions_bl(N, K)
	PT = partition_tree_bl(N, K, P, ZFI, WI)
	YS = ys_symbols_bl(N, K, P, ZFI, PT)
	YOR = Array(Array{Array{SparseMatrixCSC, 1}, 1}, N)
	for n = 1:(N - K - 1)
		YOR[n] = Array(Array{SparseMatrixCSC, 1}, 0)
	end
	for n = (N - K):N
		YSn = YS[n]
		Pn = P[n]
		Pn_L = length(Pn)
		YORn = Array(Array{SparseMatrixCSC, 1}, Pn_L)
		for p = 1:ZFI[n]
			YORnp_sparse = spzeros(1, 1)
			YORnp_sparse[1,1] = length(YSn[p])
			YORnp = Array(SparseMatrixCSC, 1)
			YORnp[1] = YORnp_sparse
			YORn[p] = YORnp
		end		
		for p = (ZFI[n] + 1):Pn_L
			YSnp = YSn[p]
			Pnp = Pn[p]
			R = length(Pnp)
			YORn[p] = yor_p(n, R, YSnp)
		end
		YOR[n] = YORn
	end
	return YOR, PT, ZFI
end

