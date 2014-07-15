
create_template(nstates::Integer) = create_template(nstates, 1.0, 0.8, 0.2)
function create_template(nstates::Integer, a::Real, b::Real, c::Real)
    x = linspace(0,1.5,nstates)
    y = a*sin(2*pi*x).*exp(-(b-x).^2/c)
    return y
end

create_signal() = create_signal(10000,0.1,10,60)
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
