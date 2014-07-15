using Distributions

function forward(V::Array{Float64,1},A::Array{Float64,2}, mu::Array{Float64,1}, sigma::Real)
    T = length(V)
    nstates = size(A,1)
    a = zeros(nstates, T)
    a[1,1] = 1
    LP = Array(Normal, nstates)
    for i=1:nstates
        LP[i] = Normal(mu[i],sigma)
    end
    aa = 0.0
    Z = 0.0
    for i=2:T
        Z = 0.0
        for j=1:nstates
            b = pdf(LP[j],V[i])
            aa = 0.0
            for k=1:nstates
                aa += a[k,i-1]*A[k,j]
            end
            a[j,i] = b*aa
            Z += a[j,i]
        end
        for j=1:nstates
            a[j,i] /= Z
        end
    end
    return a
end

function backward(V::Array{Float64,1},A::Array{Float64,2}, mu::Array{Float64,1}, sigma::Real)
    T = length(V)
    nstates = size(A,1)
    a = zeros(nstates, T)
    a[1,T] = 1
    LP = Array(Normal, nstates)
    for i=1:nstates
        LP[i] = Normal(mu[i],sigma)
    end
    aa = 0.0
    Z = 0.0
    for i=T-1:-1:1
        Z = 0.0
        for j=1:nstates
            aa = 0.0
            for k=1:nstates
                b = pdf(LP[k],V[i+1])
                aa += a[k,i+1]*A[j,k]
            end
            a[j,i] = aa 
            Z += a[j,i]
        end
        for j=1:nstates
            a[j,i] /= Z
        end
    end
    return a
end

function gamma(a,b,A,mu,sigma,x)
    nstates = size(A,1)
    n = size(a,2)
    G = zeros(nstates,nstates,n)
    LP = Array(Normal, nstates)
    for i=1:nstates
        LP[i] = Normal(mu[i],sigma)
    end
    for t=2:n
        for j=1:nstates
            for i=1:nstates
                G[i,j,t] = a[i,t-1]*A[i,j]*b[j,t]*pdf(LP[j],x[t])
            end
        end
    end
    return G
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

function test()
    A = prepA(1e-9)
end
