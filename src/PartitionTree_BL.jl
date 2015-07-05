# Parameters:
#	N::Int 
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
#	P::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 from partitions_bl()
#	ZFI::Array{Int, 1}
#	- output2 from partitions_bl()
#	WI::Array{Int, 2}
#	- output3 from partitions_bl()
# Return Values:
#	PT::Array{Array{Array{Int, 1}, 1}, 1} (Partition Tree)
#	- for each value, i, in PT[n][j], P[n][j] decomposes into P[n-1][i]
#	- length(PT[n]) = 0 for n = 1:(N - K)
#	- length(PT[n][p]) = 0 for p <= ZFI[n]
function partition_tree_bl(N::Int, K::Int, P::Array{Array{Array{Int, 1}, 1}, 1}, ZFI::Array{Int, 1}, WI::Array{Int, 2})
	PT = Array(Array{Array{Int, 1}, 1}, N)
	for n = 1:(N - K)
		PT[n] = Array(Array{Int, 1}, 0)
	end
	for n = (N - K + 1):N
		Pd = P[n - 1]
		Pn = P[n]
		PTn = Array(Array{Int, 1} , length(Pn))
		for p = 1:ZFI[n]
			PTn[p] = Array(Int, 0)
		end
		for p = (ZFI[n] + 1):length(Pn)
			PDA = pda(Pn[p])
			if length(PDA) == 1 || PDA[1][1] == PDA[2][1]
				lb = 1
				if PDA[1][1] != 1
					lb = WI[n - 1, PDA[1][1] - 1] + 1
				end
				PTn[p] = dia(PDA, Pd, lb)
			else
				lb1 = 1
				if PDA[1][1] != 1
					lb1 = WI[n - 1, PDA[1][1] - 1] + 1
				end
				lb2 = WI[n - 1, PDA[2][1] - 1]
				PTn[p] = dia(PDA, Pd, lb1, lb2)
			end
		end
		PT[n] = PTn
	end
	return PT
end

