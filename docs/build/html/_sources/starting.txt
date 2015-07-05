Getting Started
===============

Installation
------------

Method 1
^^^^^^^^
Download the setup script `InitialSetup.jl <https://github.com/GDPlumb/SnFFT.jl/blob/master/InitialSetUp.jl>`_.  
Then open a Julia session and run:

.. code-block:: julia

    julia> require("InitialSetUp.jl")

This script will install and load the SnFFT library and then run the examples.  
    
Method 2
^^^^^^^^
SnFFT can also be installed using Julia's pakage manager as follows:

.. code-block:: julia

    julia> Pkg.clone("git://github.com/GDPlumb/SnFFT.jl")
    

Set Up
------
After the installation, the user can load the library with:

.. code-block:: julia

    julia> using SnFFT
    
By default, lower level functions in the source code are not made available during installation.  
These functions are often wrapper or helper functions that may not be relevant to most users.  
However, these internal functions can be made available by modifying the file `SnFFT.jl <https://github.com/GDPlumb/SnFFT.jl/blob/master/src/SnFFT.jl>`_.  

Parallelism
-----------
*SnFFT* allows the user to compute fast Fourier transforms and inverse fast Fourier transforms in parralel with no change to their code. 
On startup, simply run:

.. code-block:: julia
	
	julia> addprocs(p)
	julia> using SnFFT
	
This will and *p* worker processes and load the SnFFT library onto them.  
Afterwards, parts of *SnFFT* will automatically run in parallel, if their heuristics think say that it will be beneficial.  

