module TestViterbi
    import Base.Test
    include("../src/viterbi.jl")
    
    create_template(nstates::Integer) = create_template(nstates, 1.0, 0.8, 0.2)
    function create_template(nstates::Integer, a::Real, b::Real, c::Real)
        x = linspace(0,1.5,nstates)
        y = a*sin(2*pi*x).*exp(-(b-x).^2/c)
        return y
    end
    
    function create_signal(N::Integer, sigma::Number,nspikes::Integer, nstates::Integer) 
        T = create_template(nstates)
        return create_signal(N,sigma,nspikes,T)
    end

    function create_signal(N::Integer, sigma::Number, nspikes::Integer, template::Array{Float64,1})
        nstates = length(template)
        S = rand(Normal(0,sigma),N)
        spike_idx = rand(1:N,nspikes)
        for i=1:nspikes
            S[spike_idx[i]:spike_idx[i] + nstates-1] += template 
        end
        return S
    end

    function viterbi(N::Integer, sigma::Number, nspikes::Integer, nstates::Integer)
        T = create_template(nstates)
        srand(1234)
        S = create_signal(N, sigma, nspikes, T)
        X = viterbi(S,nspikes/N, T, sigma*sigma)
        Q = T[X] #reconstruction
        return S,Q
    end

    function viterbi()
        S,Q = viterbi(10000,0.1,10,60)
        err = (Q-S)
        err = sum(err.*err)
        Base.Test.@test_approx_eq err 109.70934022342433
        println("Test passed. sse = $err")
    end

    function viterbi2(N::Integer, sigma::Number, nspikes::Integer, nstates::Integer)
        T1 = create_template(nstates)
        T2 = create_template(nstates, 3, 0.6, 0.1)
        srand(1235)
        S1 = create_signal(N, sqrt(2)*sigma, div(nspikes,2), T1)
        S2 = create_signal(N, sqrt(2)*sigma, div(nspikes,2), T2)
        S = S1 + S2
        T = cat(2,T1,T2)
        X = viterbi(S,nspikes/N, T, sigma*sigma)
        Q = T1[X[:,1]+1] + T2[X[:,2]+1] #reconstruction
        return S,Q
    end

    function viterbi2()
        S,Q = viterbi2(10000,0.1,10,60)
        err = S-Q
        err = sum(err.*err)
        Base.Test.@test_approx_eq err 464.267547552014
        println("Test passed. sse = $err")
    end
end
