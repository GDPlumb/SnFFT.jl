#Installs SnFFT using Julia's package manager and then runs the examples

#Open Terminal and change directories to the location of this file
#Start an interactive julia session
#Run:
#	require("InitialSetUp.jl")

Pkg.clone("git://github.com/GDPlumb/SnFFT.jl")
using SnFFT
example1()
println("")
example2()
println("")
example3()
println("")
example4()
println("")
example5()
println("")
example6()
println("")
example7()
println("")
example8()
