# Parameters:
#	N::Int
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
#	FFT::Array{Array{Float64, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft_bl()
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor_bl()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor_bl()
#	ZFI::Array{Int, 1}
#	- output3 from yor_bl()
# Return Values:
#	SNF::Array{Float64, 1}
#	- the bandlimited function over Sn that corresponds to FFT
function sn_ifft_bl(N::Int, K::Int, FFT::Array{Array{Float64, 2}, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}, ZFI::Array{Int, 1})
	np = nprocs()
	if np == 1 || N < 10 || K == 1 || K == 2 && N < 60 || K == 3 && N < 20
		SNF = Array(Float64, round(Int, factorial(N) / factorial(N - K)))
		compute_ifft_bl(N, K, SNF, FFT, YOR, PT, ZFI, Counter(1))
		return SNF
	else
		pSNFA = Array(Array{Float64, 1}, N)
		RR_FFT = Array(RemoteRef, np)
		RR_YOR = Array(RemoteRef, np)
		RR_PT = Array(RemoteRef, np)
		RR_ZFI = Array(RemoteRef, np)
		for p = 1:np
			if p != myid()
				RR_FFT[p] = RemoteRef(p)
				put!(RR_FFT[p], FFT)
				RR_YOR[p] = RemoteRef(p)
				put!(RR_YOR[p], YOR)
				RR_PT[p] = RemoteRef(p)
				put!(RR_PT[p], PT)
				RR_ZFI[p] = RemoteRef(p)
				put!(RR_ZFI[p], ZFI)
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
							pSNFA[idx] = remotecall_fetch(p, compute_sifft_bl_remote, N, K, idx, RR_FFT[p], RR_YOR[p], RR_PT[p], RR_ZFI[p])
						end
					end
				end
			end
		end
		SNF = Array(Float64, round(Int, factorial(N) / factorial(N - K)))
		BS = round(Int, factorial(N - 1) / factorial(N - K))
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
#	K::Int
#	- the problem is homogenous at N-K
#	SNF::Array{Float64, 1}
#	- the bandlimited function over Sn that is being calculated
#	FFT::Array{Array{Float64, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft_bl()
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor_bl()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor_bl()
#	ZFI::Array{Int, 1}
#	- output3 from yor_bl()
#	C::Counter
#	- a wrapper for an Int that is used to coordinate the recursive defining of the elements of SNF
# Return Values:
#	None
#	- once the recursion is finished, SNF is ful
function compute_ifft_bl(N::Int, K::Int, SNF::Array{Float64, 1}, FFT::Array{Array{Float64, 2}, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}, ZFI::Array{Int, 1}, C::Counter)
	if K == 0
		V = FFT[end][1, 1]
		V /= factorial(N)
		SNF[C.N] = V
		C.N += 1
	else
		YORn = YOR[N]
		NPn = length(YORn)
		YORd = YOR[N - 1]
		NPd = length(YORd)
		PTn = PT[N]
		ZFIn = ZFI[N]
		ZFId = ZFI[N - 1]
		sFFT = Array(Array{Float64, 2}, NPd)
		for p = 1:ZFId
			sFFT[p] = zeros(Float64, 1, 1)
		end
		for n = 1:N
			for p = (ZFId + 1):NPd
				Dim = size(YORd[p][1], 1)
				sFFT[p] = zeros(Float64, Dim, Dim)
			end
			for p = (ZFIn + 1):NPn
				update_sfft_bl!(N, n, sFFT, FFT[p], YORn[p], YORd, PTn[p], ZFId)
			end
			compute_ifft_bl(N - 1, K - 1, SNF, sFFT, YOR, PT, ZFI, C)
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
#	- Youngs Orthogonal Representations for the pth Partition of N as defined by yor_bl()
#	YORd::Array{Array{SparseMatrixCSC, 1}, 1}
#	- Youngs Orthogonal Reprsentations for the Partitions of N - 1 as defined by yor_bl()
#	PTnp::Array{Array{Array{Int, 1}, 1}, 1}
#	- The decomposition indices for the pth Partition of N as defined by yor_bl()
#	ZFId::Int
#	- if p <= ZFId, sFFT[p] is a zero-frequency component
# Return Values:
#	None
#	- sFFT is fully defined once all of the components of FFT have been used
function update_sfft_bl!(N::Int, n::Int, sFFT::Array{Array{Float64, 2}, 1}, FFTp::Array{Float64, 2}, YORnp::Array{SparseMatrixCSC, 1}, YORd::Array{Array{SparseMatrixCSC, 1}, 1}, PTnp::Array{Int, 1}, ZFId::Int)
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
		if index > ZFId
			sDim = size(YORd[index][1], 1)
			ub = lb + sDim - 1
			sFFT[index] += (Dim / (sDim * N)) * M[lb:ub, lb:ub]
			lb += sDim
		else
			lb += round(Int, YORd[index][1][1, 1])
		end
	end
end

# Parameters:
#	N::Int
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
#	n::Int
#	- determines which left-sided coset we are defining the function for
#	FFT::Array{Array{Float64, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft_bl()
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor_bl()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor_bl()
#	ZFI::Array{Int, 1}
#	- output3 from yor_bl()
# Return Values:
#	-SNF::Array{Float64, 1}
#	- the bandlimited function on the nth left-sided coset of Sn
function compute_sifft_bl(N::Int, K::Int, n::Int, FFT::Array{Array{Float64, 2}, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1},  PT::Array{Array{Array{Int, 1}, 1}, 1}, ZFI::Array{Int, 1})
	pSNF = Array(Float64, round(Int, factorial(N - 1) / factorial(N - K)))
	C = Counter(1)
	YORn = YOR[N]
	NPn = length(YORn)
	YORd = YOR[N - 1]
	NPd = length(YORd)
	PTn = PT[N]
	ZFIn = ZFI[N]
	ZFId = ZFI[N - 1]
	sFFT = Array(Array{Float64, 2}, NPd)
	for p = 1:ZFId
		sFFT[p] = zeros(Float64, 1, 1)
	end
	for p = (ZFId + 1):NPd
		Dim = size(YORd[p][1], 1)
		sFFT[p] = zeros(Float64, Dim, Dim)
	end
	for p = (ZFIn + 1):NPn
		update_sfft_bl!(N, n, sFFT, FFT[p], YORn[p], YORd, PTn[p], ZFId)
	end
	compute_ifft_bl(N - 1, K - 1, pSNF, sFFT, YOR, PT, ZFI, C)
	return pSNF
end

function compute_sifft_bl_remote(N::Int, K::Int, n::Int, RR_FFT::RemoteRef, RR_YOR::RemoteRef, RR_PT::RemoteRef, RR_ZFI::RemoteRef)
	FFT = fetch(RR_FFT)
	YOR = fetch(RR_YOR)
	PT = fetch(RR_PT)
	ZFI = fetch(RR_ZFI)
	pSNF = compute_sifft_bl(N, K, n, FFT, YOR, PT, ZFI)
	return pSNF
end

