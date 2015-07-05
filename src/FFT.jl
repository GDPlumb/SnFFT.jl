# Parameters:
#	N::Int
#	- the problem size
#	SNF::Array{Foat64, 1}
#	- SNF[i] is the value associated with the Permutation that permutation_index() maps to i
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from yor()
# Return Values:
#	FFT::Array{Float64, 2}
#	- FFT is the Fast Fourier Transform of SNF
function sn_fft(N::Int, SNF::Array{Float64, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}) 
	np = nprocs()
	if np == 1 || N < 10
		C = Counter(1)
		return compute_fft(N, SNF, YOR, PT, C)
	else #This is a modified version of Julia's pmap()
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
							SNF_subgroup = SNF[lb:ub]
							sFFT[idx] = remotecall_fetch(p, compute_fft_remote, N - 1, SNF_subgroup, RR_YOR[p], RR_PT[p], Counter(1))
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
							FFT[idx] = remotecall_fetch(p, fc_remote, N, idx, RR_YOR[p], RR_PT[p], RR_sFFT[p])
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
#	- the size of the FFT being calculated is N
#	SNF::Array{Foat64, 1}
#	- SNF[i] is the value associated with the Permutation that permutation_index() maps to i
#	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1() from yor()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2() from yor()
#	C::Counter
#	- a wrapper for an Int that is used to coordinate the recursive access to the elements of SNF
# Return Values:
#	- FFT is the Fast Fourier Transfrom of SNF
function compute_fft(N::Int, SNF::Array{Float64, 1}, YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}, C::Counter)
	if N == 1
		sFFT = Array(Array{Float64, 2}, 1)
		sFFTi = Array(Float64, 1, 1)
		sFFTi[1, 1] = SNF[C.N]
		sFFT[1] = sFFTi
		C.N += 1
		return sFFT
	end
	sFFT = Array(Array{Array{Float64, 2}, 1}, N)
	for n = 1:N
		sFFT[n] = compute_fft(N - 1, SNF, YOR, PT, C)
	end
	FFT = combine_sfft(N, YOR[N], PT[N], sFFT)
	return FFT
end

function compute_fft_remote(N::Int, SNF::Array{Float64, 1}, RR_YOR::RemoteRef, RR_PT::RemoteRef, C::Counter)
	YOR = fetch(RR_YOR)
	PT = fetch(RR_PT)
	FFT =  compute_fft(N, SNF, YOR, PT, C)
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
# Return Values:
#	FFT::Array{Array{Float64, 2}, 1}
#	- FFT is the Fast Fourier Transform for this group
function combine_sfft(N::Int, YORn::Array{Array{SparseMatrixCSC, 1}, 1}, PTn::Array{Array{Int, 1}, 1}, sFFT::Array{Array{Array{Float64, 2}, 1}, 1})
	NP = length(YORn)
	FFT = Array(Array{Float64, 2}, NP)
	for p = 1:NP
		YORnp = YORn[p]
		PTnp = PTn[p]
		FC = fc(N, YORnp, PTnp, sFFT)
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
# Return Values:
#	FC::Array{Float64, 2}
#	- FC is the Fourier Coefficient corresponding to the Partition P
function fc(N::Int, YORnp::Array{SparseMatrixCSC, 1}, PTnp::Array{Int, 1}, sFFT::Array{Array{Array{Float64, 2}, 1}, 1})
	Dim = size(YORnp[1], 1)
	FC = dsm(Dim, sFFT[N], PTnp)
	CCM = eye(Dim)
	for n = (N - 1):-1:1
		CCM = YORnp[n] * CCM
		DSM = dsm(Dim, sFFT[n], PTnp)
		FC_n = CCM * DSM
		FC += FC_n
	end
	return FC
end
	
function fc_remote(N::Int, p::Int, RR_YOR::RemoteRef, RR_PT::RemoteRef, RR_sFFT::RemoteRef)
	YOR = fetch(RR_YOR)
	PT = fetch(RR_PT)	
	sFFT = 	fetch(RR_sFFT) 
	FC = fc(N, YOR[N][p], PT[N][p], sFFT)
	return FC
end

# Parameters:
#	Dim::Int
#	- the size of the DSM that is being caculated
#	sFFTn::Array{Array{Float64, 2}, 1}
#	- one of the elements of sFFT from fc()
#	PTnp::Array{Int, 1}
#	- same as in fc()
# Return Values:
#	DSM::Array{Float64, 2} (Direct Sum Matrix) 
#	- the direct sum of the coefficients of sFFTn that correspond to the Partitions indicated by PTnp 
function dsm(Dim::Int, sFFTn::Array{Array{Float64, 2}, 1}, PTnp::Array{Int, 1})
	DSM = zeros(Float64, Dim, Dim)
	offset = 0
	for i = 1:length(PTnp)
		FC = sFFTn[PTnp[i]]
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
	end
	return DSM
end

