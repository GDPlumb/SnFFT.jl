#POTENTIAL ISSUE - line 22 may throw an Inexact Error for large problems

# Parameters:
#	N::Int
#	- N is the size of P
#	P::Array{Int, 1}
#	- P is a Partition
#Return Value:
#	NST::Int
#	- the number of Standard Tableau of the Young Diagram of P
function degree(N::Int, P::Array{Int, 1})
	if N > 20
		NST = 2432902008176640000.0
		for i = 1:length(P)
			for j = 1:P[i]
				NST /= hook_length(P, i, j)
			end
		end
		for n = 21:N
			NST *= n
		end
		NST = int(NST)
		return NST
	else
		NST = factorial(N)
		NST /= hook_product(P)
		NST = int(NST)
		return NST
	end
end

# Parameters:
#	P::Array{Int, 1}
#	- P is a Partition
# Return Values:
#	HP::Int
#	- the Hookproduct of the Young Diagram of P
function hook_product(P::Array{Int, 1})
	HP = 1
	for i = 1:length(P)
		for j = 1:P[i]
			HP *= hook_length(P, i, j)
		end
	end
	return HP
end

# Parameters:
#	P::Array{Int, 1}
#	- P is a Partition
#	i::Int 
#	- the row index
#	j::Int
#	- the column index
# Return Value
#	HL::Int
#	- the HookLength of the box in the ith row and jth colomn of the Young Diagram of P
function hook_length(P::Array{Int, 1}, i::Int, j::Int)
	R = P[i] - j
	B = 0
	i += 1
	while i <= length(P) && j <= P[i]
		i += 1
		B += 1
	end
	HL = 1 + R + B
	return HL
end

