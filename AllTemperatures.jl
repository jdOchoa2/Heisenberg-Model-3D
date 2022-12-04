include("./Wolff.jl")
using LinearAlgebra
using Random

i = parse(Int64,ARGS[1]) 
n_samples = 800
Lattices  = [8,10,12,14,16]    #Different systems avaible
CorrTimes = [3,3,3,3,3]
L  = Lattices[i]     #Lattice size
L2 = L*L             #Lattice Area
L3 = L*L2            #Lattice Volume
τ  = CorrTimes[i]    #Correlation times
τ0 = 3               #Equilibrium times
neighbours = GetNeighbours(L,L2,L3)  #Matrix with the indexes of the 6 adjacent spins of each in the lattice
# Temperatures
Ti = 1.2
Tf = 2.0
Tn = 40 
Ts = (Tf-Ti)/Tn
Temperatures = [kk for kk in Ti:Ts:Tf]
for kk in Temperatures                  #Loop over all tempetures
    save = 0; if (kk == Ti) || (kk == Tf) || (kk == 1.44)
        save = 1 end
    state = [RandomSpin() for i in 1:L3] #The 3D lattice is store as a 1D-Array
    Proceed(state, 1/kk, n_samples, L, L2, L3, neighbours, τ, τ0,save)
end