# Parameters:
#	N::Int
#	- the problem size is N
#	SNF::Array{Foat64, 1}
#	- SNF[i] is the value associated with the Permutation that permutation_index() maps to NZL[i]
#	NZL::Array{Int, 1}
#	- NZL[i] the set of NonZeroLocations for the sparse function over Sn
#	- NZL must be in increasing order
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor()
# Return Values:
#	FFT::Array{Float64, 2}
#	- FFT is the Fast Fourier Transform of SNF
function sn_fft_sp(N::Int, SNF::Array{Float64, 1}, NZL::Array{Int, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1})
	np = nprocs()
	if np == 1 || N < 10
		return compute_fft_sp(N, SNF, YOR, PT, NZL, 1, Counter(1))
	else
		sFFT = Array(Array{Array{Float64, 2}, 1}, N)
		RR_YOR = Array(RemoteRef, np)
		RR_PT = Array(RemoteRef, np)
		for p = 1:np
			if p != myid()
				RR_YOR[p] = RemoteRef(p)
				put!(RR_YOR[p], YOR)
				RR_PT[p] = RemoteRef(p)
				put!(RR_PT[p], PT)
			end
		end
		BS = factorial(N - 1)
		sSNF = Array(Array{Float64, 1}, N)
		sNZL = Array(Array{Int, 1}, N)
		SM = zeros(Int, N)
		si = 1
		ei = 1
		for n = 1:N
			if si > length(NZL)
				sSNF[n] = Array(Float64, 0)
				sNZL[n] = Array(Int, 0)
			else
				ni = NZL[ei]
				lb = (n - 1) * BS + 1
				ub = lb + BS - 1
				while lb <= ni && ni <= ub
					SM[n] = 1
					ei += 1
					if ei > length(NZL)
						break
					end
					ni = NZL[ei]
				end
				sSNF[n] = SNF[si:(ei - 1)]
				sNZL[n] = NZL[si:(ei - 1)]
				si = ei
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
							if SM[idx] == 1
								lb = (idx - 1) * BS + 1
								sFFT[idx] = remotecall_fetch(p, compute_fft_sp_remote, N - 1, sSNF[idx], RR_YOR[p], RR_PT[p], sNZL[idx], lb, Counter(1))
							else
								sFFT[idx] = Array(Array{Float64, 2}, 0)
							end
						end
					end
				end
			end
		end
		YORn = YOR[N]
		NP = length(YORn)
		PTn = PT[N]
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
							FFT[idx] = remotecall_fetch(p, fc_sp_remote, N, idx, RR_YOR[p], RR_PT[p], RR_sFFT[p], SM)
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
#	- the problem size is N
#	SNF::Array{Foat64, 1}
#	- SNF[i] is the value associated with the Permutation that permutation_index() maps to NZL[i]
#	NZL::Array{Int, 1}
#	- NZL[i] the set of NonZeroLocations for the sparse function over Sn
#	- NZL must be in increasing order
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor()
#	LB::Int
#	- the LowerBound on the value of permutation_index() for Permutations in this subgroup
#	C::Counter
#	- a wrapper for an Int that is used to coordinate the recursive access to the elements of SNF
# Return Values:
#	FFT::Array{Float64, 2}
#	- FFT is the Fast Fourier Transform of SNF
function compute_fft_sp(N::Int, SNF::Array{Float64, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}, NZL::Array{Int, 1}, LB::Int, C::Counter)
	if N == 1
		sFFT = Array(Array{Float64, 2}, 1)
		sFFTi = Array(Float64, 1, 1)
		sFFTi[1, 1] = SNF[C.N]
		sFFT[1] = sFFTi
		C.N += 1
		return sFFT
	end

	sFFT = Array(Array{Array{Float64, 2}, 1}, N)
	SM = zeros(Int, N)
	lbi = factorial(N - 1) - 1
	for n = 1:N
		if C.N > length(NZL)
			sFFT[n] = Array(Array{Float64, 2}, 0)
		else
			index = NZL[C.N]
			ub = LB + lbi
			if LB <= index && index <= ub
				SM[n] = 1
				sFFT[n] = compute_fft_sp(N - 1, SNF, YOR, PT, NZL, LB, C)
			else
				sFFT[n] = Array(Array{Float64, 2}, 0)
			end
			LB = ub + 1
		end
	end
	FFT = combine_sfft_sp(N, YOR[N], PT[N], sFFT, SM)
	return FFT
end

function compute_fft_sp_remote(N::Int, SNF::Array{Float64, 1}, RR_YOR::RemoteRef, RR_PT::RemoteRef, NZL::Array{Int, 1}, LB::Int, C::Counter)
	YOR = fetch(RR_YOR)
	PT = fetch(RR_PT)
	FFT = compute_fft_sp(N, SNF, YOR, PT, NZL, LB, C)
	return FFT
end

# Parameters:
#	N::Int
#	- the size of the FFT being calculated is N
#	YORn::Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- Young's Orthogonal Representations for the Partitions of N
#	PTn::Array{Array{Int, 1}, 1}
#	- the decomposition indices for the Partitions of N
#	sFFT::Array{Array{Array{Float64, 2}, 1}, 1}
#	- the array of N FFT's of size N-1 that is used to compute this FFT of size N
#	SM::Array{Int, 1}
#	- SM[n] == 0 if sFFT[n] corresponded to a homogenous subgroup with a value of 0
#	- SM[n] == 1 otherwise
# Return Values:
#	FFT::Array{Array{Float64, 2}, 1}
#	- FFT is the Fast Fourier Transform for this group
function combine_sfft_sp(N::Int, YORn::Array{Array{SparseMatrixCSC, 1}, 1}, PTn::Array{Array{Int, 1}, 1}, sFFT::Array{Array{Array{Float64, 2}, 1}, 1}, SM::Array{Int, 1})
	NP = length(YORn)
	FFT = Array(Array{Float64, 2}, NP)
	for p = 1:NP
		YORnp = YORn[p]
		PTnp = PTn[p]
		FC = fc_sp(N, YORnp, PTnp, sFFT, SM)
		FFT[p] = FC
	end
	return FFT
end

# Parameters:
#	N::Int
#	- the Fourier Component that is being calculated corresponds to a Partition, P, of N
#	YORnp::Array{SparseMatrixCSC, 1}
#	- Young's Orthogonal Representations for P
#	PTnp::Array{Int, 1}
#	- the indices for the Partitions that P decomposes into
#	sFFT::Array{Array{Array{Float64, 2}, 1}, 1}
#	- the array of N FFT's of size N-1 that is used to compute this FFT of size N
#	SM::Array{Int, 1}
#	- SM[n] == 0 if sFFT[n] corresponded to a homogenous subgroup with a value of 0
#	- SM[n] == 1 otherwise
# Return Values:
#	FC::Array{Float64, 2}
#	- FC is the Fourier Coefficient corresponding to the Partition P
function fc_sp(N::Int, YORnp::Array{SparseMatrixCSC, 1}, PTnp::Array{Int, 1}, sFFT::Array{Array{Array{Float64, 2}, 1}, 1}, SM::Array{Int, 1})
	Dim = size(YORnp[1], 1)
	if SM[N] == 1
		FC = dsm(Dim, sFFT[N], PTnp)
	else
		FC = zeros(Float64, Dim, Dim)
	end
	CCM = eye(Dim)
	for n = (N - 1):-1:1
		CCM = YORnp[n] * CCM
		if SM[n] == 1
			DSM = dsm(Dim, sFFT[n], PTnp)
			FC_n = CCM * DSM
			FC += FC_n
		end
	end
	return FC
end

function fc_sp_remote(N::Int, p::Int, RR_YOR::RemoteRef, RR_PT::RemoteRef, RR_sFFT::RemoteRef, SM::Array{Int, 1})
	YOR = fetch(RR_YOR)
	PT = fetch(RR_PT)	
	sFFT = 	fetch(RR_sFFT) 
	FC = fc_sp(N, YOR[N][p], PT[N][p], sFFT, SM)
	return FC
end
