module SnFFT

<<<<<<< HEAD
#load the files included
include("Convolution.jl")
include("Correlation.jl")
include("Counter.jl")
include("Element.jl")
include("Examples.jl")
include("FFT.jl")
include("FFT_BL.jl")
include("FFT_SP.jl")
include("Function.jl")
include("Hook.jl")
include("IFFT.jl")
include("IFFT_BL.jl")
include("IFFT_P.jl")
include("Partitions.jl")
include("Partitions_BL.jl")
include("PartitionTree.jl")
include("PartitionTree_BL.jl")
include("YamanouchiSymbols.jl")
include("YamanouchiSymbols_BL.jl")
include("YoungsOrthogonalRepresentations.jl")
include("YoungsOrthogonalRepresentations_BL.jl")

export
	#core functions
	yor,
	yor_bl,
	sn_fft,
	sn_fft_bl,
	sn_fft_sp,
	sn_ifft,
	sn_ifft_bl,
	sn_ifft_p,
	sn_convolution,
	sn_correlation,
	
	#examples
	example1,
	example2,
	example3,
	example4,
	example5,
	example6,
	example7,
	example8,
	example_clustering, 

	#functions from SnElement.jl
	sn_multiply,
	sn_inverse,
	sn_p,
	sn_cc,
	sn_at,
	sn_t,
	permutation_ccf,
	ccf_index,
	permutation_index,
	index_ccf,
	ccf_permutation,
	index_permutation,
	permutation_atf,
	yor_permutation,

	#functions fromSnFunction.jl
	snf,
	snf_bl,
	snf_sp, 

	#other fucntions Examples.jl
	mallowsdistribution,
	preferencematrix,
	kendalldistance,
	permutation_string,
	partition_string,
	
	#other
	partitions

end 

=======
# package code goes here

end # module
>>>>>>> c2921e6a62dd7f9c4ff629c01803d0b504b1f3fa
