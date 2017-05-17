function get_threshold(X::Array{Float64,1})
  n = length(X)
  μ = mean(X)
  σ = std(X)::Float64 #because of type instability in std
  l = μ - 6*σ 
  u = μ + 6*σ 
  σ2 = 0.0
  for x in X
    if l < x < u
      σ2 +=  x*x
    end
  end
  σ2 = sqrt(σ2/n - μ*μ)
  μ, σ2
end

function extract_spikes(X::Array{Float64,1};μ::Float64=NaN, σ::Float64=NaN,nq::Tuple{Int64, Int64}=(20,40))
  if isnan(σ) || isnan(μ)
    μ, σ = get_threshold(X)
  end
  n = length(X)
  l = μ - 6*σ
  u = μ + 6*σ
  pidx = Array{Int64}(0)
  xmin = Inf
  xmax = -Inf
  j = 1
  in_peak = false
  _idx = 0
  while j <= n
    #check whether points are outside μ ± σ
    #fast forward points that are within noise levels
    x = X[j]
    if x < l
      if x < xmin
        xmin = x
        _idx = j
      end
      in_peak = true
    elseif x > u
      if x > xmax
        xmax = x
        _idx = j
      end
      in_peak = true
    else
      if in_peak == true #just exited a peak; reset
        push!(pidx, _idx)
      end
      xmin = Inf
      xmax = -Inf
      in_peak = false
    end
    j += 1
  end
  waveforms = extract_spikes(X, pidx,nq)
  pidx, waveforms
end

function extract_spikes(X::Array{Float64,1}, idx::Array{Int64,1}, n::Tuple{Int64, Int64}=(20,40))
  nt = sum(n)
  nx = length(idx)
  nb,na = n
  waveforms = zeros(nt, nx)
  for i in 1:nx
    _idx = idx[i]
    for j in -nb:na-1
      waveforms[j+nb+1,i] = X[_idx + j]
    end
  end
  waveforms
end

