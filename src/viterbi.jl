using Distributions
function viterbi{T<:Real}(y::Array{T,1}, p::Number, mu::Array{Float64,1}, C::Number)
    nstates = length(mu)
    A = prepA(p,nstates)
    x, T2,T1 = viterbi_2(y,A,mu,C)
    return x
end

function viterbi{T<:Real}(y::Array{T,1}, p::Number, mu::Array{Float64,2}, C::Number,chunk::Number)
    nstates = size(mu,1)
    ntemps = size(mu,2)
    A = prepA(p,nstates)
    Q = createIndex(ntemps,nstates)
    chunks = collect(1:chunk:length(y))
    X = zeros(Int16,length(y),ntemps)
    #make sure we pick up the last chunk as well
    if chunks[end] != length(y)
        push!(chunks,length(y))
    end
    #process array in chunks
    i = 1
    while chunks[i+1] < length(y)
        println(string("Sorting chunk ",chunks[i],":",chunks[i+1]-1))
        x,T2,T1 = viterbi_3(y[chunks[i]:chunks[i+1]],A,mu,C)
        #check to see whether the last point is noise
        #if not, go backwards until we hit noise, then resort the non-noise part in the next chunk
        if maximum(x) == 1
            i+= 1
            continue #pure noise; go to the next chunk
        end
        j = 0 
        while x[end-j] == 1 && j < chunk
            j+=1
        end
        chunks[i+1] -= j
        if chunks[i+1] <= chunks[i]
            i+=1
            continue #gone through the whole chunk
        end
        X[chunks[i]:chunks[i+1],:] = unfoldIndex(x[1:(end-j)],Q)
        #rewind next chunk
        i+=1
    end
    return X
end

function viterbi{T<:Real}(y::Array{T,1}, p::Number, mu::Array{Float64,2}, C::Number)
    nstates = size(mu,1)
    ntemps = size(mu,2)
    A = prepA(p,nstates)
    Q = createIndex(ntemps,nstates)
    x,T2,T1 = viterbi_3(y,A,mu,C)
    X = unfoldIndex(x,Q) 
    return X
end

function viterbi(y::Array{Int16,1}, A::Array{Float64,2}, mu::Array{Float64,1}, C::Float64)
    #straight forward implementation of the Viterbi algorithm
    #assume gaussian emission probabilities
    nstates = size(A,1)
    nstates == size(mu,1) || throw(ArgumentError("size(B,2) should be equal to the number of states"))
    nobs = size(y,1)
    x = zeros(Int16, nobs)
    T1 = log(zeros(Float64, nstates, nobs))
    T2 = zeros(Int16, nstates, nobs)
    #define the emitting distributions; assume the same covariance matrix
    Pn = [Normal(mu[i],sqrt(C))  for i=1:size(mu,1)]
    #define initial elements
    #for i=1:nstates
    #    T1[i,1] = log(P[i])+logpdf(Pn[i],y[1])
    #end
    T1[1,1] = 0
    for i=2:nobs
        for j=1:nstates    
            tm = -Inf
            km = 0
            q = logpdf(Pn[j],y[i]) #emission probability for symbol i from state j
            for k=1:nstates
                #probability of transitioning from state k to state j and emitting observation i
                t = T1[k,i-1]+A[k,j]+q
                #j==1 && k==60 && println(string("1->60 $t"))
                #j==1 && k==1 && println(string("1->1 $t"))
                if t > tm
                    tm = t
                    km = k
                end
            end
            T1[j,i] = tm
            T2[j,i] = km
        end
    end
    #define the last state
    #mx,x[end] = findmax(T1[:,end])
    x[end] = 1
    #run backward
    for i=nobs:-1:2
        x[i-1] = T2[x[i],i]
    end
    return x,T2, T1
end

function viterbi_2{T<:Real}(y::Array{T,1}, A::Array{Float64,2}, mu::Array{Float64,1}, C::Float64)
    #implementation tailor made for ring transition probabilities
    #assume gaussian emission probabilities
    nstates = size(A,1)
    nstates == size(mu,1) || throw(ArgumentError("size(B,2) should be equal to the number of states"))
    nobs = size(y,1)
    x = zeros(Int16, nobs)
    T1 = log(zeros(Float64, nstates, nobs))
    T2 = zeros(Int16, nstates, nobs)
    #define the emitting distributions; assume the same covariance matrix
    Pn = [Normal(mu[i],sqrt(C))  for i=1:size(mu,1)]
    T1[1,1] = 0
    for i=2:nobs
        #for state j=1, only two transitions are poossible; 1->1 and nstates->1
        j = 1
        q = logpdf(Pn[j],y[i]) #emission probability for symbol i from state j
        t1_1 = T1[1,i-1]+A[1,j]+q
        tl_1 = T1[60,i-1]+A[60,j]+q
        if t1_1 > tl_1
            T1[j,i] = t1_1
            T2[j,i] = 1
        else
            T1[j,i] = tl_1
            T2[j,i] = 60 
        end
        #remaining transitions are deterministic
        for j=2:nstates
            q = logpdf(Pn[j],y[i]) #emission probability for symbol i from state j
            tm = T1[j-1,i-1]+A[j-1,j]+q
            T1[j,i] = tm
            T2[j,i] = j-1
        end
    end
    #define the last state
    x[end] = 1
    #run backward
    for i=nobs:-1:2
        x[i-1] = T2[x[i],i]
    end
    return x,T2, T1
end

function prepA(p,n)
    A = log(zeros(n,n))
    A[1,1] = log(1.-p)
    A[1,2] = log(p)
    for i=2:n-1
        A[i,i+1] = 0
    end
    A[n,1] = 0
    return A
end

function viterbi_3{T<:Real}(y::Array{T,1}, A::Array{Float64,2}, mu::Array{Float64,2}, C::Float64)
    #implementation tailor made for ring transition probabilities; accepts multiple neurons
    #assume gaussian emission probabilities
    nstates = size(A,1)
    ntemps = size(mu,2) #number of templates
    nstates == size(mu,1) || throw(ArgumentError("size(B,2) should be equal to the number of states"))
    ntotstates = 1 + ntemps*(nstates-1) + div(ntemps*(ntemps-1),2)*(nstates-1)*(nstates-1)
    nobs = size(y,1)
    x = zeros(Int16, nobs)
    T1 = log(zeros(Float64, ntotstates, nobs))
    T2 = zeros(Int16, ntotstates, nobs)
    #define the emitting distributions; assume the same covariance matrix
    Pn = Array(Distribution, nstates,ntemps)
    S = sqrt(C)
    for i=1:ntemps
        for j=1:nstates
            Pn[j,i] = Normal(mu[j,i],S)
        end
    end
    LL = zeros(nstates,ntemps)
    #get the mapping from individual states for each neuron to a global state index
    joint_states = createIndex(ntemps,nstates)
    pairs = zeros(Int16,ntemps,ntemps)
    ppq = 1 #variable to track which pair we are looking at
    for k1=1:ntemps-1
        for k2=k1+1:ntemps
            pairs[k1,k2] = ppq
            ppq += 1
        end
    end
    pairs = pairs + pairs'
    T1[1,1] = 0
    for i=2:nobs
        #loop through each pair
        for k1=1:ntemps-1
            for k2=k1+1:ntemps
                #both neurons are in state 1
                #for state j=1, only two transitions are poossible; 1->1 and nstates->1
                ppq = pairs[k1,k2]
                j1 = 1
                j2 = 1
                #4 possibilities to consider as each neuron could have come into state 1 either by staying
                #in the first state or coming in from the last state.
                jstate1 = joint_states[ppq,j1,j2] #get the joint state
                jstate2 = joint_states[ppq, nstates, nstates]
                jstate3 = joint_states[ppq, 1, nstates]
                jstate4 = joint_states[ppq, nstates, 1]
                q1 = logpdf(Pn[1,k1],y[i]) #emission probability for symbol i from state j
                q2 = logpdf(Pn[1,k2],y[i]) #emission probability for symbol i from state j
                LL[1,k1] = q1
                LL[1,k2] = q2
                t1 = T1[jstate1,i-1]+A[1,1]+q1 + q2 + A[1,1] #neuron 1 and 2 from state 1
                t2 = T1[jstate2,i-1]+A[nstates,1]+q1 + q2 + A[nstates,1] #neuron 1 and 2 from last state
                t3 = T1[jstate3,i-1]+A[1,1]+q1 + q2 + A[nstates,1] #neuron 1 from 1, neuron 2 from last state
                t4 = T1[jstate4,i-1]+A[nstates,1]+q1 + q2 + A[1,1] #neuron 1 from last state, neuron 2 from state 1
                mx,xi = findmax([t1,t2,t3,t4])
                T1[jstate1,i] = mx
                T2[jstate1,i] = [jstate1,jstate2,jstate3,jstate4][xi]
                #neuron 2 in state 1 neuron 1 in other states
                for j1=2:nstates
                    q1 = logpdf(Pn[j1,k1],y[i]) #emission probability for symbol i from state j
                    LL[j1,k1] = q1
                    jstate1 = joint_states[ppq,j1-1,1] #get the joint state
                    jstate2 = joint_states[ppq, j1-1, nstates]
                    jstatef = joint_states[ppq,j1,1]
                    tm1 = T1[jstate1,i-1]+A[1,1]+q2 + A[j1-1,j1] + q1
                    tm2 = T1[jstate2,i-1]+A[nstates,1]+q2 +  A[j1-1,j1] + q1
                    if tm1 > tm2
                        T1[jstatef,i] = tm1 #neuron 1 transitioned from state 1, neuron 2 transitions to state j2
                        T2[jstatef,i] = jstate1
                    else
                        T1[jstatef,i] = tm2 #neuron 1 transitioned from state 60
                        T2[jstatef,i] = jstate2
                    end
                end
                #neuron 1 in state 1
                #remaining transitions are deterministic
                q1 = LL[1,k1] #emission probability for symbol i from state j
                for j2=2:nstates    
                    q2 = logpdf(Pn[j2,k2],y[i]) #emission probability for symbol i from state j
                    jstate1 = joint_states[ppq,1,j2-1] #get the joint state
                    jstate2 = joint_states[ppq, nstates, j2-1]
                    jstatef = joint_states[ppq,1,j2]
                    tm1 = T1[jstate1,i-1]+A[j2-1,j2]+q2 + A[1,1] + q1
                    tm2 = T1[jstate2,i-1]+A[j2-1,j2]+q2 +  A[nstates,1] + q1
                    if tm1 > tm2
                        T1[jstatef,i] = tm1 #neuron 1 transitioned from state 1, neuron 2 transitions to state j2
                        T2[jstatef,i] = jstate1
                    else
                        T1[jstatef,i] = tm2 #neuron 1 transitioned from state 60
                        T2[jstatef,i] = jstate2
                    end
                    #neuron 1 in other states as well
                    a2 = A[j2-1,j2] + q2 #only compute once for the following loop
                    for j1=2:nstates
                        q1 = LL[j1,k1]
                        jstatep = joint_states[ppq,j1-1,j2-1]
                        jstatef = joint_states[ppq,j1,j2]
                        tm = T1[jstatep,i-1] + A[j1-1,j1] + q1 + a2
                        T1[jstatef,i] = tm
                        T2[jstatef,i] = jstatep
                    end
                end
            end
        end
    end
    #define the last state
    mx, x[end] = findmax(T1[:,end])
    #run backward
    for i=nobs:-1:2
        x[i-1] = T2[x[i],i]
    end
    return x,T2, T1
end

function createIndex(N,K)
    #create a linear index for N rings, each with K states
    npairs = div(N*(N-1),2)
    Q = zeros(Int64,npairs,K,K)
    KK = (K-1)*(K-1)
    NK = N*(K-1)+1
    NP = NK + KK
    i = 1
    for k1=1:N-1
        for k2=k1+1:N
            Q[i,2:end,2:end] = reshape([1:KK]+NK + (i-1)*NP, (1,K-1,K-1))
            Q[i,1,:] = [1:K] + (k1-1)*(K-1)
            Q[i,2:end,1] = [2:K] + (k2-1)*(K-1)
            i += 1
        end
    end
    return Q
end

function unfoldIndex(states,stateMatrix)
    npairs,nstates,nstates = size(stateMatrix)
    ncells = div(int(1+sqrt(1+4*2*npairs)),2) #get the number of cells
    #cheat
    pairs = zeros(Int16,ncells,ncells)
    k = 1
    for k1=1:ncells-1
        for k2=k1+1:ncells
            pairs[k1,k2] = k
            k+=1
        end
    end
    #pairs = pairs + pairs'
    X = zeros(Int64,length(states),ncells)
    for i=1:length(states)
        pair,state1,state2 = ind2sub((npairs,nstates,nstates), find(stateMatrix.==states[i]))
        cell1,cell2 = ind2sub((ncells,ncells),find(pairs.==pair[1]))
        X[i,cell1] = state1[1]
        X[i,cell2] = state2[1]
    end
    return X-1
end
