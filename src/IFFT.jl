# Parameters:
#	N::Int
#	- the problem size
#	FFT::Array{Array{Float64, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft() or sn_fft_sp()
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor()
# Return Values:
#	SNF::Array{Float64, 1}
#	- the function over Sn that corresponds to FFT
function sn_ifft(N::Int, FFT::Array{Array{Float64, 2}, 1},  YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1})
	np = nprocs()
	if np == 1 || N < 10
		SNF = Array(Float64, factorial(N))
		compute_ifft(N, SNF, FFT, YOR, PT, Counter(1))
		return SNF
	else
		pSNFA = Array(Array{Float64, 1}, N)
		RR_FFT = Array(RemoteRef, np)
		RR_YOR = Array(RemoteRef, np)
		RR_PT = Array(RemoteRef, np)
		for p = 1:np
			if p != myid()
				RR_FFT[p] = RemoteRef(p)
				put!(RR_FFT[p], FFT)
				RR_YOR[p] = RemoteRef(p)
				put!(RR_YOR[p], YOR)
				RR_PT[p] = RemoteRef(p)
				put!(RR_PT[p], PT)
			end
		end
		i = 1
		nextidx() = (idx = i; i += 1; idx)
		@sync begin
			for p = 1:np
				if p != myid()
					@async begin
						while true
							idx = nextidx()
							if idx > N
								break
							end
							pSNFA[idx] = remotecall_fetch(p, compute_sifft_remote, N, idx, RR_FFT[p], RR_YOR[p], RR_PT[p])
						end
					end
				end
			end
		end
		SNF = Array(Float64, factorial(N))
		BS = factorial(N - 1)
		i = 1
		for n = 1:N
			pSNF = pSNFA[n]
			for si = 1:BS
				SNF[i] = pSNF[si]
				i += 1
			end
		end
		return SNF
	end
end

# Parameters:
#	N::Int
#	- the problem size
#	SNF::Array{Float64, 1}
#	- the function over Sn that is being calculated
#	FFT::Array{Array{Float64, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft() or sn_fft_sp()
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor()
#	C::Counter
#	- a wrapper for an Int that is used to coordinate the recursive defining of the elements of SNF
# Return Values:
#	None
#	- once the recursion is finished, SNF is full
function compute_ifft(N::Int, SNF::Array{Float64, 1}, FFT::Array{Array{Float64, 2}, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}, C::Counter)
	if N == 2
		YORn = YOR[2]
		sFFT = 0.5 * (full(YORn[1][1]) * FFT[1] + full(YORn[2][1]) * FFT[2])
		SNF[C.N] = sFFT[1,1]
		C.N += 1
		sFFT = 0.5 * (FFT[1] + FFT[2])
		SNF[C.N] = sFFT[1,1]
		C.N += 1
	else
		YORn = YOR[N]
		NPn = length(YORn)
		YORd = YOR[N - 1]
		NPd = length(YORd)
		PTn = PT[N]
		sFFT = Array(Array{Float64, 2}, NPd)
		for n = 1:N
			for p = 1:NPd
				Dim = size(YORd[p][1], 1)
				sFFT[p] = zeros(Float64, Dim, Dim)
			end
			for p = 1:NPn
				update_sfft!(N, n, sFFT, FFT[p], YORn[p], YORd, PTn[p])
			end
			compute_ifft(N - 1, SNF, sFFT, YOR, PT, C)
		end
	end
end

# Parameters:
#	N::Int
#	- the problem size
#	n::Int
#	- determines which left-sided coset we are defining the FFT for
#	sFFT::Array{Array{Float64, 2}, 1}
#	- the FFT of the nth left-sided coset that we are defining
#	FFTp::Array{Float64, 2}
#	- the pth component of the FFT of size N
#	YORnp::Array{SparseMatrixCSC, 1}, 1}
#	- Youngs Orthogonal Representations for the pth Partition of N
#	YORd::Array{Array{SparseMatrixCSC, 1}, 1}
#	- Youngs Orthogonal Reprsentations for the Partitions of N - 1
#	PTnp::Array{Array{Array{Int, 1}, 1}, 1}
#	- The decomposition indices for the pth Partition of N
# Return Values:
#	None
#	- sFFT is fully defined once all of the components of FFT have been used
function update_sfft!(N::Int, n::Int, sFFT::Array{Array{Float64, 2}, 1}, FFTp::Array{Float64, 2}, YORnp::Array{SparseMatrixCSC, 1}, YORd::Array{Array{SparseMatrixCSC, 1}, 1}, PTnp::Array{Int, 1})
	Dim = size(YORnp[1], 1)
	M = eye(Float64, Dim, Dim)
	for ccn = n:(N - 1)
		M = M * YORnp[ccn]
	end
	M = transpose(M)
	M = M * FFTp
	lb = 1
	for d = 1:length(PTnp)
		index = PTnp[d]
		sDim = size(YORd[index][1], 1)
		ub = lb + sDim - 1
		sFFT[index] += (Dim / (sDim * N)) * M[lb:ub, lb:ub]
		lb += sDim
	end
end

# Parameters:
#	N::Int
#	- the problem size
#	n::Int
#	- determines which left-sided coset we are defining the function for
#	FFT::Array{Array{Float64, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft() or sn_fft_sp()
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor()
# Return Values:
#	pSNF::Array{Float64, 1}
#	- the function on the nth left-sided coset of Sn
function compute_sifft(N::Int, n::Int, FFT::Array{Array{Float64, 2}, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1},  PT::Array{Array{Array{Int, 1}, 1}, 1})
	pSNF = Array(Float64, factorial(N - 1))
	C = Counter(1)
	YORn = YOR[N]
	NPn = length(YORn)
	YORd = YOR[N - 1]
	NPd = length(YORd)
	PTn = PT[N]
	sFFT = Array(Array{Float64, 2}, NPd)
	for p = 1:NPd
		Dim = size(YORd[p][1], 1)
		sFFT[p] = zeros(Float64, Dim, Dim)
	end
	for p = 1:NPn
		update_sfft!(N, n, sFFT, FFT[p], YORn[p], YORd, PTn[p])
	end
	compute_ifft(N - 1, pSNF, sFFT, YOR, PT, C)
	return pSNF
end

function compute_sifft_remote(N::Int, n::Int, RR_FFT::RemoteRef, RR_YOR::RemoteRef, RR_PT::RemoteRef)
	FFT = fetch(RR_FFT)
	YOR = fetch(RR_YOR)
	PT = fetch(RR_PT)
	pSNF = compute_sifft(N, n, FFT, YOR, PT)
	return pSNF
end 

