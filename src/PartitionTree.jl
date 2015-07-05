# Parameters:
#	N::Int
#	- the problem size
#	P::Array{Array{Array{Int, 1}, 1}, 1}
#	- output1 from partitions()
#	WI::Array{Int, 2}
#	- output2 from partitions()
# Return Values:
#	PT::Array{Array{Array{Int, 1}, 1}, 1} (Partition Tree)
#	- for each value, i, in PT[n][p], P[n][p] decomposes into P[n-1][i]
#	- length(PT[1]) = 0
function partition_tree(N::Int, P::Array{Array{Array{Int, 1}, 1}, 1}, WI::Array{Int, 2})
	PT = Array(Array{Array{Int, 1}, 1}, N)
	PT[1] = Array(Array{Int, 1}, 0)
	for n = 2:N
		Pd = P[n - 1]
		Pn = P[n]
		PTn = Array(Array{Int, 1} , length(Pn))
		for p = 1:length(Pn)
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

# Parameters:
#	P::Array{Int, 1}
#	- P is a Partition
# Return Values:
#	PDA::Array{Array{Int, 1}, 1} (Partition Decomposition Array)
#	- PDA[i] is the the ith Partition that P decomposes into 
function pda(P::Array{Int, 1})
	P_L = length(P)    
	num = 1
	for i = 1:(P_L - 1)
		if P[i] > P[i + 1]
			num += 1
		end
	end
	PDA = Array(Array{Int, 1}, num)	
	c = 1
	for i = 1:(P_L - 1)
		if P[i] > P[i + 1]
			D = copy(P)
			D[i] -= 1
			PDA[c] = D
			c += 1
		end
	end
	if P[P_L] == 1
		D = copy(P[1:(P_L - 1)])
		PDA[c] = D
	else
		D = copy(P)
		D[P_L] -= 1
		PDA[c] = D
	end
	return PDA
end

# Parameters:
#	PDA::Array{Array{Int, 1}, 1}
#	- output1 from pda()
#	Pd::Array{Array{Int, 1}, 1}
#	- Pd is array containing all the partitions of the same size as those in PDA
#	lb::Int
#	- the lower bound for the index where the first element of PDA will be found
# Return Values:
#	DIA::Array{Int, 1} (Decomposition Index Array)
#	- PDA[i] = Pd[DIA[i]]
function dia(PDA::Array{Array{Int, 1}, 1}, Pd::Array{Array{Int, 1}, 1}, lb::Int)
	PDA_L = length(PDA)
	DIA = Array(Int, PDA_L)
	c = 1
	p = lb
	while true
		if Pd[p] == PDA[c]
			DIA[c] = p
			c += 1
			if c > PDA_L
				return DIA
			end
		end
		p += 1
	end
end

# Parameters:
#	PDA::Array{Array{Int, 1}, 1}
#	- output1 from pda()
#	Pd::Array{Array{Int, 1}, 1}
#	- Pd is array containing all the partitions of the same size as those in PDA
#	lb1::Int
#	- the lower bound for the index where the first element of PDA will be found
#	lb2::Int
#	- the lower bound for the index where the second element of PDA will be found
# Return Values:
#	DIA::Array{Int, 1}
#	- PDA[i] = Pd[DIA[i]]
function dia(PDA::Array{Array{Int, 1}, 1}, Pd::Array{Array{Int, 1}, 1}, lb1::Int, lb2::Int)
	PDA_L = length(PDA)
	DIA = Array(Int, PDA_L)
	p = lb1
	while true
		if Pd[p] == PDA[1]
 			DIA[1] = p
			break
		end
		p += 1
	end
	p = lb2
	c = 2
	while true
		if Pd[p] == PDA[c]
			DIA[c] = p
			c += 1
			if c > PDA_L
				return DIA
			end
		end
		p += 1
	end
end

