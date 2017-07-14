function highpass_filter(X::Array{Float64,1},sampling_rate::Float64, cutoff::Float64=300.0;method=Butterworth(4))
  filtfilt(digitalfilter(Highpass(cutoff;fs=sampling_rate), method), X)
end

function get_threshold(X::Array{Float64,1},θ=6.0)
  n = length(X)
  μ = mean(X)
  σ = std(X)::Float64 #because of type instability in std
  l = μ - θ*σ
  u = μ + θ*σ
  σ2 = 0.0
  for x in X
    if l < x < u
      σ2 +=  x*x
    end
  end
  σ2 = sqrt(σ2/n - μ*μ)
  μ, σ2
end

function extract_spikes(X::Array{Float64,1};μ::Float64=NaN, σ::Float64=NaN,nq::Tuple{Int64, Int64}=(20,40),θ=6.0,only_negative=true)
  if isnan(σ) || isnan(μ)
    μ, σ = get_threshold(X,θ)
  end
  n = length(X)
  l = μ - θ*σ
  u = μ + θ*σ
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
  elseif (x > u ) && !only_negative
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
  pidx, waveforms = extract_spikes(X, pidx,nq)
  pidx, waveforms
end

function extract_spikes(X::Array{Float64,1}, idx::Array{Int64,1}, n::Tuple{Int64, Int64}=(20,40))
  nt = sum(n)
  nb,na = n
  idxl = find(x->x>nb, idx)
  nx = length(idxl)
  waveforms = zeros(nt, nx)
  for i in 1:nx
     _idx = idx[idxl[i]]
    for j in -nb:na-1
      waveforms[j+nb+1,i] = X[_idx + j]
    end
  end
  idx[idxl], waveforms
end

