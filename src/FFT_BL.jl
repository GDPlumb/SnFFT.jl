# Parameters:
#	N::Int
#	- the problem size N
#	K::Int
#	- the problem is homogenous at N-K
#	SNF::Array{Foat64, 1}
#	- SNF[i] is the value assigned to the ith homogenous subgroup of size N-K
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor_bl()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor_bl()
#	ZFI::Array{Int, 1}
#	- output3 from yor_bl()
# Return Values:
#	FFT::Array{Float64, 2}
#	- FFT is the Fast Fourier Transform of SNF
function sn_fft_bl(N::Int, K::Int, SNF::Array{Float64, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}, ZFI::Array{Int, 1})
	np = nprocs()
	if np == 1 || N < 10 || K == 1 || K == 2 && N < 60 || K == 3 && N < 20
		return compute_fft_bl(N, K, SNF, YOR, PT, ZFI, Counter(1))
	else 
		sFFT = Array(Array{Array{Float64, 2}, 1}, N)
		RR_YOR = Array(RemoteRef, np)
		RR_PT = Array(RemoteRef, np)
		RR_ZFI = Array(RemoteRef, np)
		for p = 1:np
			if p != myid()
				RR_YOR[p] = RemoteRef(p)
				put!(RR_YOR[p], YOR)
				RR_PT[p] = RemoteRef(p)
				put!(RR_PT[p], PT)
				RR_ZFI[p] = RemoteRef(p)
				put!(RR_ZFI[p], ZFI)
			end
		end
		BS = round(Int, length(SNF) / N)
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
							lb = (idx - 1) * BS + 1
							ub = lb + BS - 1
							sSNF = SNF[lb:ub]
							sFFT[idx] = remotecall_fetch(p, compute_fft_bl_remote, N - 1, K - 1, sSNF, RR_YOR[p], RR_PT[p], RR_ZFI[p], Counter(1))
						end
					end
				end
			end
		end
		YORn = YOR[N]
		NP = length(YORn)
		PTn = PT[N]
		ZFId = ZFI[N - 1]
		np = nprocs()
		RR_sFFT = Array(RemoteRef, np)
		for p = 1:np
			if p != myid()
				RR_sFFT[p] = RemoteRef(p)
				put!(RR_sFFT[p], sFFT)
			end
		end
		FFT = Array(Array{Float64, 2}, NP)
		i = 1
		nextidx() = (idx = i; i += 1; idx)
		@sync begin
			for p = 1:np
				if p != myid()
					@async begin
						while true
							idx = nextidx()
							if idx > NP
								break
							end
							FFT[idx] = remotecall_fetch(p, fc_bl_remote, N, idx, RR_YOR[p], RR_PT[p], RR_sFFT[p], ZFId)
						end
					end
				end
			end
		end
		return FFT
	end
end

# Parameters:
#	N::Int
#	- the problem size N
#	K::Int
#	- the problem is homogenous at N-K
#	SNF::Array{Foat64, 1}
#	- SNF[i] is the value assigned to the ith homogenous subgroup of size N-K
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor_bl()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor_bl()
#	ZFI::Array{Int, 1}
#	- output3 from yor_bl()
#	C::Counter()
#	- a wrapper for an Int that is used to coordinate the recursive access to the elements of SNF
# Return Values:
#	FFT::Array{Float64, 2}
#	- FFT is the Fast Fourier Transform of SNF
function compute_fft_bl(N::Int, K::Int, SNF::Array{Float64, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}, ZFI::Array{Int, 1}, C::Counter)
	sFFT = Array(Array{Array{Float64, 2}, 1}, N)
	if K != 1
		for n = 1:N
			sFFT[n] = compute_fft_bl(N - 1, K - 1, SNF, YOR, PT, ZFI, C)
		end
	else
		for n = 1:N
			V = SNF[C.N]
			sFFT[n] = compute_fft_iv(N - 1, V, YOR[N - 1])
			C.N += 1 
		end
	end
	FFT = combine_sfft_bl(N, YOR[N], PT[N], sFFT, ZFI[N], ZFI[N - 1])
	return FFT
end

function compute_fft_bl_remote(N::Int, K::Int, SNF::Array{Float64, 1}, RR_YOR::RemoteRef, RR_PT::RemoteRef, RR_ZFI::RemoteRef, C::Counter)
	YOR = fetch(RR_YOR)
	PT = fetch(RR_PT)
	ZFI = fetch(RR_ZFI)
	FFT = compute_fft_bl(N, K, SNF, YOR, PT, ZFI, C)
	return FFT
end

# Parameters:
#	N::Int
#	- the size of the homogenous (invariant) subgroup
#	V::Float64
#	- the value associated with this homogenous subgroup
#	YORn::Array{SparseMatrixCSC, 1}
#	- Young's Orthogonal Representations for the Partitions of N as defined by yor_bl()
function compute_fft_iv(N::Int, V::Float64, YORn::Array{Array{SparseMatrixCSC, 1}, 1})
	NP = length(YORn)
	FFT = Array(Array{Float64, 2}, NP)
	for p = 1:(NP - 1)
		FFT[p] = full(YORn[p][1])
	end
	NZ = Array(Float64, 1, 1)
	NZ[1,1] = factorial(Int128(N)) * V
	FFT[NP] = NZ
	return FFT
end

# Parameters:
#	N::Int
#	- the size of the FFT being calculated is N
#	YORn::Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- Young's Orthogonal Representations for the Partitions of N as defined by yor_bl()
#	PTn::Array{Array{Int, 1}, 1}
#	- the decomposition indices for the Partitions of N as defined by partition_tree_bl()
#	sFFT::Array{Array{Array{Float64, 2}, 1}, 1}
#	- the array of N FFT's of size N-1 that is used to compute this FFT of size N
# Return Values:
#	FFT::Array{Array{Float64, 2}, 1}
#	- FFT is the Fast Fourier Transform for this group
function combine_sfft_bl(N::Int, YORn::Array{Array{SparseMatrixCSC, 1}, 1}, PTn::Array{Array{Int, 1}, 1}, sFFT::Array{Array{Array{Float64, 2}, 1}, 1}, ZFIn::Int, ZFId::Int)
	NP = length(YORn)
	FFT = Array(Array{Float64, 2}, NP)
	for p = 1:ZFIn
		FFT[p] = full(YORn[p][1])
	end
	for p = (ZFIn + 1):NP
		YORnp = YORn[p]
		PTnp = PTn[p]
		FC = fc_bl(N, YORnp, PTnp, sFFT, ZFId)
		FFT[p] = FC
	end
	return FFT
end

# Parameters:
#	N::Int
#	- the Fourier Component that is being calculated corresponds to a Partition, P, of N
#	YORnp::Array{SparseMatrixCSC, 1}
#	- Young's Orthogonal Representations for P as defined by yor_bl()
#	PTnp::Array{Int, 1}
#	- the indices for the Partitions that P decomposes into as defined by pt_bl()
#	sFFT::Array{Array{Array{Float64, 2}, 1}, 1}
#	- the array of N FFT's of size N-1 that is used to compute this FFT of size N
# Return Values:
#	FC::Array{Float64, 2}
#	- FC is the Fourier Coefficient corresponding to the Partition P
function fc_bl(N::Int, YORnp::Array{SparseMatrixCSC, 1}, PTnp::Array{Int, 1}, sFFT::Array{Array{Array{Float64, 2}, 1}, 1}, ZFId::Int)
	Dim = size(YORnp[1], 1)
	FC = dsm_bl(Dim, sFFT[N], PTnp, ZFId)
	CCM = eye(Dim)
	for n = (N - 1):-1:1
		CCM = YORnp[n] * CCM
		DSM = dsm_bl(Dim, sFFT[n], PTnp, ZFId)
		FC_n = CCM * DSM
		FC += FC_n
	end
	return FC
end

function fc_bl_remote(N::Int, p::Int, RR_YOR::RemoteRef, RR_PT::RemoteRef, RR_sFFT::RemoteRef, ZFId::Int)
	YOR = fetch(RR_YOR)
	PT = fetch(RR_PT)	
	sFFT = 	fetch(RR_sFFT)
	FC = fc_bl(N, YOR[N][p], PT[N][p], sFFT, ZFId)
	return FC
end
	
# Parameters:
#	Dim::Int
#	- the size of the DSM that is being caculated
#	sFFTn::Array{Array{Float64, 2}, 1}
#	- one of the elements of sFFT from fc_bl()
#	PTnp::Array{Int, 1}
#	- same as in fc_bl()
# Return Values:
#	DSM::Array{Float64, 2} (Direct Sum Matrix)
#	- the direct sum of the coefficients of sFFTn that correspond to the Partitions indicated by PTnp 
function dsm_bl(Dim::Int, sFFTn::Array{Array{Float64, 2}, 1}, PTnp::Array{Int, 1}, ZFId::Int)
	DSM = zeros(Float64, Dim, Dim)
	offset = 0
	for i = 1:length(PTnp)
		index = PTnp[i]
		if index > ZFId
			FC = sFFTn[index]
			sDim = size(FC, 1)
			for r = 1:sDim
				for c = 1:sDim
					v = FC[r, c]
					if v != 0.0
						DSM[offset + r, offset + c] = v
					end
				end
			end
			offset += sDim
		else
			FC = sFFTn[index]
			offset += round(Int, FC[1, 1])
		end
	end
	return DSM
end
