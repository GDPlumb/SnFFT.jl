# Parameters:
#	N::Int 
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
# Return Values:
#	P::Array{Array{Array{Int, 1}, 1}, 1} (Partitions)
#	- P[n][p] contains the pth Partition of n that is required for the bandlimited functionality
#	- This is equivalent to the Partitions that correspond to non-zero frequencies in FFT and their immediate descendents
#	- length(P[n]) = 0 for n = 1:(N-K-1)
#	ZFI::Array{Int, 1} (Zero Frequency Information)
#	- ZFI[n] = k if, for p<=k, P[n][p] is a zero frequency partition
#	WI::Array{Int, 2} (Width Information)
#	- WI[n, w] contains the number of Partitions of n whose first element is less than or equal to w
#	- This only counts the Partitions that are in P and not all Partitions
function partitions_bl(N::Int, K::Int)
	P_K, WI_K = partitions(K)
	WI = zeros(Int, N, N)
	for n = (N - K):N
		WI[n,n] = 1
	end	
	for n = (N - K):(N - 1)
		for w = (N - K - 1):(n - 1)
			WI[n,w] = WI_K[n - w, min(n - w, w)]
		end
	end
	for w = (N - K):(N - 1)
		WI[N,w] = WI_K[N - w, min(N - w, w)]
	end
	for n = (N - K):N
		for w = (N - K):n
			WI[n, w] += WI[n, w - 1]
		end
	end
	P = Array(Array{Array{Int, 1}, 1}, N)
	ZFI = Array(Int, N)
	for n = 1:(N - K - 1)
		P[n] = Array(Array{Int, 1}, 0)
		ZFI[n] = 0
	end
	for n = (N - K):(N - 1)
		Pn = Array(Array{Int, 1}, WI[n,n])
		i = 1
		for spi = 1:WI_K[n - N + K + 1, min(n - N + K + 1, N - K - 1)]
			Pn[i] = [N - K - 1; P_K[n - N + K + 1][spi]]
			i += 1
		end
		ZFI[n] = i - 1
		for w = (N - K):(n - 1)
			for spi = 1:WI_K[n - w, min(n - w, w)] 
				Pn[i] = [w; P_K[n - w][spi]]
				i += 1
			end
		end
		Pn[i] = [n]
		P[n] = Pn
	end
	Pn = Array(Array{Int, 1}, WI[N, N])
	i = 1
	for w = (N - K):(N - 1)
		for spi = 1:WI_K[N - w, min(N - w, w)] 
			Pn[i] = [w; P_K[N - w][spi]]
			i += 1
		end
	end
	Pn[i] = [N]
	P[N] = Pn
	ZFI[N] = 0
	return P, ZFI, WI
end

