#Finds the representation of a permutation for a particular partition
function example1()
	example1(5,[3, 2],[2,3,1,4,5])
end

# Parameters:
#	N::Int
#	- the problem size
#	partition::Array{Int, 1}
#	- a partition of N
#	permutation::Array{Int, 1}
#	- a permutation of N
function example1(N::Int, partition::Array{Int, 1}, permutation::Array{Int, 1})

	P, WI = partitions(N) #Creates the partitions of i for i is 1 through N
	p = 1 #Find the index of the specified partitions
	while P[N][p] != partition
		p += 1
	end

	YS = ys_partition(N, partition) #Creates the Yamanouchi Symbols of the specified partition
	degree = length(YS)

	YOR, PT = yor(N) #Creates Young's Orthogonal Representations and the Pointer Tree for partition decomposistion
	YORnp = YOR[N][p] #Get the Representation of each Adjacent Transposition of the specified partition

	RM = yor_permutation(permutation, YORnp) #Creates the Representation Matrix of the specified permutation for the specified partition

	println("Finds the representation of a permutation for a particular partition")
	println("")
	println("Permutation: ")
	println(permutation_string(permutation))
	println("")
	println("Partition: ")
	println(partition_string(partition))
	println("")
	println("Degree: ")
	println(degree)
	println("")
	println("Standard Tableaux: ")
	for i = 1:degree #Converts the Yamanouchi Symbols into Standard Tableau and prints them
		print_ys(YS[i])
		println("")
	end
	println("Representation matrix of the permutation for this partition: ")
	print(full(round(RM,5)))
end

#If a and b are permutations and c = a * b, demonstrates that the representation of c is the representation of a multiplied by the representation of b
function example2()
	example2(5,[3,2],[1,2,4,3,5],[2,3,1,4,5])
end

# Parameters:
#	N::Int
#	- the problem size
#	partition::Array{Int, 1}
#	- a partition of N
#	p1::Array{Int, 1}
#	- the first permutation of N
#	p2::Array{Int, 1}
#	- the second permutation of N
function example2(N::Int, partition::Array{Int, 1}, p1::Array{Int, 1}, p2::Array{Int, 1})

	P, WI = partitions(N) #Creates the partitions of i for i is 1 through N
	p = 1 #Find the index of the specified partitions
	while P[N][p] != partition
		p += 1
	end

	pm = sn_multiply(p1, p2) #Find the product of the two permutations

	YOR, PT = yor(N) #Creates Young's Orthogonal Representations and the Pointer Tree for partition decomposistion
	YORnp = YOR[N][p] #Get the Representation of each Adjacent Transposition of the specified partition

	RM1 = yor_permutation(p1, YORnp) #Find the representation of the first permutation
	RM2 = yor_permutation(p2, YORnp) #Find the representation of the second permutation
	RMm = yor_permutation(pm, YORnp) #Find the representation of the product of the permutations
	error = RM1 * RM2 - RMm #Find the error in the process

	println("If a and b are permutations and c = a * b, demonstrates that the representation of c is the representation of a multiplied by the representation of b")
	println("")
	ST1 = permutation_string(p1)
	println("Representation of ", ST1)
	println(full(round(RM1, 5)))
	println("")
	ST2 = permutation_string(p2)
	println("Representation of ", ST2)
	println(full(round(RM2, 5)))
	println("")
	STm = permutation_string(pm)
	println(ST1," * ", ST2, " is ", STm)
	println("")
	println("Representation of ", STm)
	println(full(round(RMm, 5)))
	println("")
	println("Maximum error: ")
	println(maximum(abs(error)))
end

#Demonstrates how to put the values of a function on Sn in the order used to compute the FFT
function example3()

	PA = Array(Array{Int, 1}, 24) #Create an array to hold all of the permutations of 1,2,3,4
	PA[1] = [1,2,3,4]
	PA[2] = [1,2,4,3]
	PA[3] = [1,3,2,4]
	PA[4] = [1,3,4,2]
	PA[5] = [1,4,2,3]
	PA[6] = [1,4,3,2]
	PA[7] = [2,1,3,4]
	PA[8] = [2,1,4,3]
	PA[9] = [2,3,1,4]
	PA[10] = [2,3,4,1]
	PA[11] = [2,4,1,3]
	PA[12] = [2,4,3,1]
	PA[13] = [3,1,2,4]
	PA[14] = [3,1,4,2]
	PA[15] = [3,2,1,4]
	PA[16] = [3,2,4,1]
	PA[17] = [3,4,1,2]
	PA[18] = [3,4,2,1]
	PA[19] = [4,1,2,3]
	PA[20] = [4,1,3,2]
	PA[21] = [4,2,1,3]
	PA[22] = [4,2,3,1]
	PA[23] = [4,3,1,2]
	PA[24] = [4,3,2,1]

	VA = Array(Float64, 24) #Create the values of the function on Sn
	for i = 1:24
		VA[i] = Float64(i)
	end

	SNF = snf(4,PA,VA) #Put the values into the order used to compute the FFT

	println("Demonstrates how to put the values of a function on Sn in the order used to compute the FFT")
	println("")
	println("Start with any ordering: ")
	println("Permutation   Function value for that Permutation")
	for i = 1:24
		ST = permutation_string(PA[i])
		ST = string(ST, "  ", VA[i])
		println(ST)
	end
	println("")
	println("Is sorted into this order for SnFFT to use: ")
	println("Permutation   Function value for that Permutation")
	for i = 1:24
		ni = round(Int, SNF[i])
		ST = permutation_string(PA[ni])
		ST = string(ST, "  ", VA[ni])
		println(ST)
	end
end

#Demonstrates how to compute the FFT, the IFFT, and shows that the IFFT recovers the initial function
function example4()
	example4(4)
end

# Parameters:
#	N::Int
#	- the problem size
function example4(N::Int)

	YOR, PT = yor(N) #Creates Young's Orthogonal Representations and the Pointer Tree for Partition decomposistion
	SNF = 2 * rand(factorial(N)) #Create a random function on Sn
	FFT = sn_fft(N, SNF, YOR, PT) #Find the FFT
	iSNF = sn_ifft(N, FFT, YOR, PT) #Find the inverse FFT
	dif = SNF - iSNF #Find the error in the process

	P, WI = partitions(N) #Create the Partitions of 1:N
	Pn = P[N] #Get the partitions that act as the labels for the irreducibles


	println("Demonstrates how to compute the FFT, the IFFT, and shows that the IFFT recovers the initial function")
	println("")
	println("A random function on Sn: ")
	for i = 1:length(SNF)
		println(SNF[i])
	end
	println("")
	println("The Fourier Transform: ")
	for i = 1:length(FFT)
		println(partition_string(Pn[i]))
		println(round(FFT[i], 4))
		println("")
	end
	println("Maximum error in recovering the original function on Sn: ")
	println(maximum(abs(dif)))
	end

#Demonstrates how to compute a sparse FFT and shows that it produces the same results as the normal FFT
function example5()
	example5(4,0.8)
end

# Parameters:
#	N::Int
#	- the problem size
#	SC::Float64
#	- the portion of the function that is zero-valued
function example5(N::Int, SC::Float64)

	#create a sparse function on Sn
	SNF = rand(factorial(N))
	L = 0
	for(i = 1:length(SNF))
		if(SNF[i] < SC) #determine if the value will be zero
			SNF[i] = 0.0
		else #if non-zero, assign it a random value
			L += 1
			SNF[i] = 2 * rand()
		end
	end
	sSNF = Array(Float64, L)
	NZL = Array(Int, L)
	index = 1
	for i = 1:length(SNF)
		if SNF[i] != 0.0
			sSNF[index] = SNF[i]
			NZL[index] = i
			index += 1
		end
	end

	YOR, PT = yor(N) #Creates Young's Orthogonal Representations and the Pointer Tree for Partition decomposistion

	nFFT = sn_fft(N, SNF, YOR, PT) #Compute the FFT normally
	sFFT = sn_fft_sp(N, sSNF, NZL, YOR, PT) #Compute the sparse FFT

	P, WI = partitions(N) #Create the Partitions of 1:N
	Pn = P[N] #Get the partitions that act as the labels for the irreducibles

	println("Demonstrates how to compute a sparse FFT and shows that it produces the same results as the normal FFT")
	println("")
	println("A random sparse function on Sn")
	for i = 1:length(SNF)
		println(SNF[i])
	end
	println("")
	println("The sparse Fourier Transform: ")
	for i = 1:length(sFFT)
		println(partition_string(Pn[i]))
		println(round(sFFT[i], 4))
		println("")
	end
	println("Comparing the results of the two methods of computing the FFT using '==' evaluates to: ")
	println(sFFT == nFFT)
end

#Demonstrates the property of the FFT where the components of the FFT are the representations when the FFT is computed on a delta function
function example6()
	example6(5,[1,3,2,4,5])
end

# Parameters:
#	N::Int
#	- the problem size
#	Permutation::Array{Int, 1}
#	- a permutation of N
function example6(N::Int, Permutation::Array{Int, 1})

	YOR, PT = yor(N) #Creates Young's Orthogonal Representations and the Pointer Tree for Partition decomposistion
	YORn = YOR[N]

	#Create a delta function on Sn
	SNF = [1.0]
	NZL = [permutation_index(Permutation)]

	FFT = sn_fft_sp(N, SNF, NZL, YOR, PT) #Find the FFT, because delta functions are almost completely sparse the sparse FFT is used

	error = 0.0
	for i = 1:length(FFT)
		YORp = yor_permutation(Permutation, YORn[i]) #Find the representation for the selected permutation for this partition

		#Find the error
		error_p = maximum(abs(YORp - FFT[i]))
		if error_p > error
			error = error_p
		end
	end

	P, WI = partitions(N) #Create the Partitions of 1:N
	Pn = P[N] #Get the partitions that act as the labels for the irreducibles

	println("Demonstrates the property of the FFT where the components of the FFT are the representations when the FFT is computed on a delta function")
	println("")
	println("The permutation for which the funciton is non zero: ", permutation_string(Permutation))
	println("")    
	println("The Fourier transform: ")
	for i = 1:length(FFT)
		println(partition_string(Pn[i]))
		println(round(FFT[i], 4))
		println("")
	end
	println("Maximum absolute difference between the Fourier components and the YOR representations for the permutation: ")
	println(error)
end

#Creates a random Band Limited funtion and returns both a compressed version (that FFT_BL() can use) and a full version (that FFT() can use)
#computing a bandlimited FFT and iFFT
function example7()
	example7(4,2)
end

# Parameters:
#	N::Int 
#	- the problem size
#	K::Int
#	- the problem is homogenous at N-K
function example7(N::Int, K::Int)

	blSNF = 2 * rand(round(Int, factorial(N) / factorial(N - K)))
	nSNF = Array(Float64, factorial(N))
	index = 1
	for i = 1:length(blSNF)
		v = blSNF[i]
		for j = 1:factorial(N - K)
			nSNF[index] = v
			index += 1
		end
	end

	YOR, PT = yor(N)  #Creates Young's Orthogonal Representations and the Pointer Tree for Partition decomposistion
	nFFT = sn_fft(N, nSNF, YOR, PT) #compute the normal FFT

	YOR, PT, ZFI = yor_bl(N, K)  #Creates Young's Orthogonal Representations and the Pointer Tree for Partition decomposistion for bandlimited functions
	blFFT = sn_fft_bl(N, K, blSNF, YOR, PT, ZFI) #compute the band limited FFT

	iblSNF = sn_ifft_bl(N, K, blFFT, YOR, PT, ZFI) #compute the band limited iFFT

	P, ZFI, WI = partitions_bl(N,K) #Create the Partitions of 1:N
	Pn = P[N] #Get the partitions that act as the labels for the irreducibles

	error = 0.0
	shift = length(nFFT) - length(blFFT)
	for i = 1:length(blFFT) 
		#Find the error
		error_p = maximum(abs(nFFT[i + shift] - blFFT[i]))
		if error_p > error
			error = error_p
		end
	end

	println("Creates a random Band Limited funtion and returns both a compressed version (that FFT_BL() can use) and a full version (that FFT() can use) computing a bandlimited FFT and iFFT")
	println("")
	println("The full representation of a bandlimited function on Sn:")
	for i = 1:length(nSNF)
		println(nSNF[i])
	end
	println("")
	println("The compressed representation of the same function used by SnFFT:")
	for i = 1:length(blSNF)    
		println(blSNF[i])
	end
	println("")
	println("The band limited Fourier transform: ")
	for i = 1:length(blFFT)
		println(partition_string(Pn[i]))
		println(round(blFFT[i], 4))
		println("")
	end

	println("The maximum absolute difference between the full and bandlimited FFT: ")
	println(error)
	println("")
	println("The absolute error in recovering the band limited function on Sn from the inverse bandlimited FFT: ")
	println(maximum(abs(iblSNF - blSNF)))
end

#Calculating a partial iFFT
function example8()
	example8(4,3)
end

# Parameters:
#	N::Int 
#	- the problem size
#	M::Int
#	- the number of the top frequencies of the Fourier transform to use
function example8(N::Int, M::Int)

	SNF = rand(factorial(N)) #create a random function on Sn
	YOR, PT = yor(N) #create YOR matrics and the Pointer Tree
	FFT = sn_fft(N, SNF, YOR, PT) #calculate the normal FFT (the sparse FFT also works)
	piSNF = sn_ifft_p(N, M, FFT, YOR, PT) #calculate the partial iFFT

	println("Demonstrates how to calculate a partial iFFT")
	println("")
	println("The original random function on Sn:")
	for i = 1:length(SNF)    
		println(SNF[i])
	end
	println("")
	println("The recovered function from the partial iFFT considering the lowest ", M, " frequencies: ")
	for i = 1:length(piSNF)    
		println(piSNF[i])
	end
end

#Uses sparcl (an R library) to perform sparse Hierarchical Clustering
#The output is saved into your home directory
function example_clustering()
	example_clustering(2,5,1.0,7)
end

# Parameters:
#	C::Int
#	- number of clusters
#	S::Int
#	- sizes of the clusters
#	F::Float64
#	- spreading factors for the distributions
#	N::Int
#	- the problem size
# Return Values:
# None - saves the data to CSV files and runs the clustering script
# Notes: For Julia V0.4
# - the TmStruct and related code can be replaced with "now()"
# - the appropriate version number needs to be used to define $loc
function example_clustering(C::Int,S::Int,F::Float64,N::Int)

	mkdir = "mkdir"
	Time = now() 
	dir = string("ClusteringExample-",Time)
	R = "R"	
	CMD = "CMD"
	BATCH = "BATCH"
	loc = string(homedir(),"/.julia/v0.4/SnFFT/src/sparcl_script.R") 

	run(`$mkdir $dir`)
	cd(dir)
	clustering_setup(C,S,F,N)
	run(`$R $CMD $BATCH $loc`)
end

#########################################
# Functionality needed for the examples #
#########################################

###
# ClusteringExample
###

#Sets up the data for the clustering
function clustering_setup(C::Int64, S::Int, F::Float64, N::Int64)

	#Construct the distribution Centers	
	Centers = Array(Array{Int64, 1}, C)
	for i = 1:C	
		Centers[i] = sn_p(N)
		println(Centers[i])
	end

	#Construct the Distributions
	Permutations = Array(Array{Int64, 1}, C * S)
	Distributions = Array(Array{Float64, 1}, C * S)
	k = 1
	for i = 1:C
		center = Centers[i]
		for j = 1:S
			rt = sn_t(N)
			perm = sn_multiply(rt, center)
			Permutations[k] = perm
			Distributions[k] = mallowsdistribution(perm, F)
			println("Distribution ", k)
			k += 1
		end
	end 

	#Calculate the Fourier Transforms
	FTA = Array(Array{Array{Float64, 2}, 1}, C * S)
	YOR, PT = yor(N)	
	for i = 1:(C * S)
		FTA[i] = sn_fft(N, Distributions[i], YOR, PT)
		println("FFT ", i)
	end

	X = Array(Float64, C * S, factorial(N))
	for c = 1:(C * S)
		k = 1
		FT = FTA[c]
		for p = 1:length(FT)
			FC = FT[p]
			d = size(FC, 1)
			for i = 1:d
				for j = 1:d
					X[c, k] = FC[i,j]
					k += 1
				end
			end
		end
	end

	writecsv("centers.csv",Centers)
	writecsv("perms.csv",Permutations)
	writecsv("input-x.csv",X)
end

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
#	Gamma::Float64
#	- the spreading factor
# Return Values:
#	MD::Array{Float64, 1}
#	- the mallows distribution with spreading factor Gamma centered at P
function mallowsdistribution(P::Array{Int64, 1}, Gamma::Float64)
	N = length(P)
	N_F = factorial(N)
	Q = preferencematrix(P)
	MD = Array(Float64, N_F)
	Gamma_neg = -1 * Gamma
	V = (1 - exp(Gamma_neg)) ^ N
	for j = 1:N
		V /= (1 - exp(j * Gamma_neg))
	end
	for i = 1:N_F
		Pi = index_permutation(N, i)
		Qi = preferencematrix(Pi) #Find the Preference Matrix of that permutation
		MD[i] = V * exp(Gamma_neg * kendalldistance(Qi, Q)) #Find the value of the Mallow Distribution for that permutation
	end
	return MD
end

# Parameters:
#	P::Array{Int, 1}
#	- a permutation
# Return Values:
#	Q::Array{Int, 2}
#	- the preference matrix for P
#	- Q[i,j] = 1 if and only if j precedes i in P
function preferencematrix(P::Array{Int64, 1})
	N = length(P)
	Q = zeros(Int64, N, N)
	for i = 1:N
		for j = 1:N
			p = P[j]
			if p == i
				break
			else
				Q[i, p] = 1
			end
		end
	end
	return Q
end

# Parameters:
#	Q1::Array{Int, 2}
#	- the preference matrix for the first permutation
#	Q2::Array{Int, 2}
#	- the preference matrix for the second permutation
# Return Values:
#	D::Int
#	- the Kendall Tau Distance between two permutations
#	- the the number of pairs (i, j) such that: P1[i] < P1[j] and P2[i] > P2[j]
function kendalldistance(Q1::Array{Int64, 2}, Q2::Array{Int64, 2})
	N = size(Q1, 1)
	D = 0
	for i = 1:N
		for j = 1:N
			D += abs(Q1[i,j] - Q2[i,j])
		end
	end
	D /= 2
	D = round(Int, D)
	return D
end


###
# Printing methods
###

# Parameters:
#	Permutation::Array{Int, 1}
#	- a permutation
# Return Values:
#	ST::String
#	- the string representation of permutation
function permutation_string(Permutation::Array{Int, 1})
	ST = "[ "
	for i = 1:length(Permutation)
		ST = string(ST, Permutation[i], " ")
	end
	ST = string(ST, "] ")
	return ST
end

# Parameters:
#	Partition::Array{Int, 1}
#	- a partition
# Return Values:
#	ST::String
#	- the string representation of partition
function partition_string(Partition::Array{Int, 1})
	ST = "( "
	for i = 1:length(Partition)
		ST = string(ST, Partition[i], " ")
	end
	ST = string(ST, ") ")
	return ST
end

# Parameters:
#	YS::Array{Int, 1}
#	- a yamanouchi symbol
# Return Values:
#	None - prints YS
function print_ys(YS::Array{Int, 1})
	N = length(YS)
	numRows = maximum(YS)
	ST = Array(AbstractString, numRows)
	IA = ones(Int, numRows)
	for i = 1:numRows
		ST[i] = ""
	end
	for n = 1:N	
		row = YS[n]
		NI = IA[row]
		ST[row] = string(ST[row]," ", n)
		NI += 1
		IA[row] = NI
	end
	for i = 1:numRows
		println(ST[i])
	end
end

