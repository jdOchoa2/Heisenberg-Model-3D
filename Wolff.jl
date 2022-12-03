function RandomSpin()
    """Creates a random unit vector from an homogenous distribution
    on the unit sphere"""
    θ = acos(1-rand(0:1e-6:2))
    ϕ = rand(0:1e-6:2π)
    return [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
end

function AcceptanceProbability(Spin1, Spin2, n, β)
    """Computes the probability of adding Spin2 - neighbour of
    Spin1 - to the cluster. It uses the temperature β and the 
    vector of reflection n"""
    return 1 - exp(min(0,-2β*transpose(Spin1)*n*transpose(Spin2)*n))
end

function ReflectSpin(S,n,state)
    """Reflectes the spin state[S] around the plane defined by the
    normal vector n"""
    state[S] -=  2(transpose(state[S])*n)*n
end

function GetNeighbours(L,L2,L3)
    """Gets the coordinates of the 6 neighbours of every spin in the state. It uses 
    periodic boundery conditions in a 3D lattice of length L (area L2, volume L3). 
    The L3 × 6 matriz it returns contains the neighbours as folows: 
    [Up-Down-Left-Right-Front-Back]"""
    Neighbours = zeros(Int16,L3,6)
    for x in 1:L3
        Neighbours[x,:] = [L2*((x-1)÷L2)+(x-1-L+L2)%L2+1,
            L2*((x-1)÷L2)+(x-1+L)%L2+1,
            L*((x-1)÷L)+(x+L-2)%L+1,
            L*((x-1)÷L)+(x)%L+1,
            (x-L2+L3-1)%L3+1,
            (x+L2-1)%L3+1]
    end
    return Neighbours
end

function Grow_Reflect(S, Cluster, n, β, L, L2, L3, state, neighbours)
    """Checks every neighbour of spin state[S] and adds them to the Cluster given
        a certain probability (refair to AcceptanceProbability). Finally, it 
        reflects the spin around the normal vector n (refair to ReflectSpin)"""
    #Checks neighbours
    for Sn in neighbours[S,:]
        if Sn ∉ Cluster && rand(0:1e-15:1) < AcceptanceProbability(state[S],state[Sn],n,β)
            #Adds neighbour to Cluster
            push!(Cluster,Sn)
            #Checks neighbours of neighbour
            Grow_Reflect(Sn, Cluster, n, β, L, L2, L3, state, neighbours)
        end
    end
    #Reflects the spin
    ReflectSpin(S,n, state)
end

function NewState(state, L, L2, L3, β, neighbours)
    """Chooses a randon spin from the state with lattice sice L, and also a
    random normal vector. From there it builds the cluster using Wolff algorithm
    (see Grow_Reflect)"""
    n = RandomSpin()     #Initial random Normal Vector
    S0 = rand(1:1:L3)    #Initial random Spin
    Cluster = [S0]       #Stores the indexes of the spins added to the cluster
    Grow_Reflect(S0, Cluster, n, β, L, L2, L3, state, neighbours) #Builds cluster and flips for new state
    return length(Cluster)
end

function Magnetization(state)
    """Returns the magnetization as the norm of the sum of all spins in the grid"""
    return norm(sum(state))
end

function Energy(state,L3,neighbours)
    """Returns the energy of a state"""
    energy = 0.0
    #Goes through every spin in the lattice
    for ii in 1:L3  
        Sp = [0,0,0]
        #Only checks 3 out the 6 neighbourds to avoid double counting of bounderies
        for Sn in neighbours[ii,[2,4,6]] 
           Sp += state[Sn]
        end
        energy -= dot(Sp,state[ii])
    end
    return energy
end

function Proceed(state, β, n_samples, L, L2, L3, neighbours,τ)
    """For a single value of β in a lattice of size L, creates a sequence of 
    states and meassures the observables for the de-correlated configurations.
    Then it prints the desired quantities."""
    Mmean,M2mean,M4mean,Emean,E2mean,p = zeros(6)
    progress = 0      #Samples taken
    while progress < n_samples
        M = Magnetization(state); E = Energy(state,L3,neighbours)
        Mmean += M; M2mean += M*M; M4mean += M^4 
        Emean += E; E2mean += E*E
        t = 0         #Current virtual time in MCSS
        while t < τ 
            mcss = 0  #Number of spins flipped
            while mcss < L3 
                mcss += NewState(state, L, L2, L3, β, neighbours) 
            end 
        t += 1 
        end
        p += NewState(state, L, L2, L3, β, neighbours) #Additional step just to get the cluster size
        progress += 1 
    end
    Mmean /= n_samples; M2mean /= n_samples; M4mean /= n_samples
    Emean /= n_samples; E2mean /= n_samples; p /= n_samples
    χ  = β*(M2mean-Mmean*Mmean)
    Cv = β*β*(E2mean-Mmean*Mmean)
    U  = 1. - 1. /3. *(M4mean/(M2mean*M2mean))
    print(β," ",Mmean," ",M2mean," ",M4mean," ",Emean," ",E2mean," ",p," ",χ," ",Cv," ",U,"\n")
end