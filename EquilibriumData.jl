include("./Wolff.jl")
using LinearAlgebra
using Random

i = parse(Int64,ARGS[1])
Lattices  = [8,10,14,20,24]    #Different systems avaible

βc = 0.693        #Critical temperature
Ttime = 1000      #Max time 
N_samples = 1000  #Number of simulations performed 
L = Lattices[i]   #Lattice lentgh 
L2 = L*L          #Lattice area
L3 = L2*L         #Lattice volume

neighbours = GetNeighbours(L,L2,L3)  #Matrix with the indexes of the 6 adjacent spins of each in the lattice

Equilibrium(βc, Ttime, N_samples, L, L2, L3, neighbours)