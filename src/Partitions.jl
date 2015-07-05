# Let X be a Partition of N
#	X::Array{Int, 1}
#	X[i] > 0 for all i
#	sum(X) == N
#	X[i] >= X[j] when i < j

# Parameters:
#	N::Int 
#	- the problem size
# Return Values:
#	P::Array{Array{Array{Int, 1}, 1}, 1} (Partitions)
#	- P[n][p] contains the pth Partition of n
#	WI::Array{Int, 2} (Width Information)
#	- WI[n, w] contains the number of Paritions of n whose first element is less than or equal to w
function partitions(N::Int)
	WI = zeros(Int, N, N) 
	for n = 1:N
		WI[n, 1] = 1
	end
	for w = 2:N
		WI[w, w] = 1
		for n = (w + 1):N
			num = 0
			for i = 1:min(n - w, w)
				num += WI[n - w, i]
			end
			WI[n, w] = num
		end
	end	
	for n = 2:N
		for w = 2:n
			WI[n, w] += WI[n, w - 1]
		end
	end	
	P = Array(Array{Array{Int, 1},1}, N)
	Pn = Array(Array{Int, 1}, 1)
	Pn[1] = [1]
	P[1] = Pn
	for n = 2:N
		Pn = Array(Array{Int, 1}, WI[n, n])
		Pn[1] = ones(Int, n)
		i = 2
		for w = 2:(n - 1)
			for spi = 1:WI[n - w, min(n - w, w)] 
				Pn[i] = [w; P[n - w][spi]]
				i += 1
			end
		end
		Pn[i] = [n]
		P[n] = Pn
	end
	return P, WI
end

