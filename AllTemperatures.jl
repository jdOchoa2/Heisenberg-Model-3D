include("./Wolff.jl")
using LinearAlgebra
using Random

L = parse(Int64,ARGS[1])  #Lattice size
n_samples = parse(Int64,ARGS[2]) 
L2 = L*L   #Lattice Area
L3 = L*L2  #Lattice Volume
τ = 10
#Random.seed!(51)
neighbours = GetNeighbours(L,L2,L3)  #Matrix with the indexes of the 6 adjacent spins of each in the lattice
Temperatures = [kk for kk in 0.7:.1:2.1 ]
βs = 1 ./ Temperatures;

for kk in βs                         #Loop over all tempetures
    state = [RandomSpin() for i in 1:L3] #The 3D lattice is store as a 1D-Array
    Proceed(state, kk, n_samples, L, L2, L3, neighbours,τ)
end