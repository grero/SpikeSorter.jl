import Grid

@doc meta("Return the width in miliseconds of the spike", returns = (Float64,) ) ->
function spike_width(spike::Array{Float64,1},samplingrate::Real)
	spikei = Grid.InterpGrid(spike,Grid.BCnil, Grid.InterpQuadratic)
	dt = length(spike)/samplingrate
	#upsample
	spikes = spikei[linspace(1,length(spike),1000)]
	idx = find(spikes .< 0.5*minimum(spikes))
	return dt*(idx[end]	- idx[1])
end

function spike_width(TF::TemplateFile,samplingrate::Real)
	#find the largest
	w = zeros(size(TF.templates,3))
	for i in 1:size(TF.templates,3)
		idx = indmin(TF.templates[:,:,i][:]) 
		pt,ch = ind2sub((size(TF.templates,1), size(TF.templates,2)),idx)
		w[i] = spike_width(TF.templates[:,ch,i],samplingrate)
	end
	w
end

spike_width(A) = spike_width(A, 40000)

function spike_dv_ratio(spike::Array{Float64,1},samplingrate::Real)
	spikei = Grid.InterpGrid(spike,Grid.BCnil, Grid.InterpQuadratic)
	dt = length(spike)/samplingrate
	#upsample
	spikes = spikei[linspace(1,length(spike),1000)]
	dv = diff(spikes)
	abs(maximum(dv)/minimum(dv))
end

spike_dv_ratio(spike::Array{Float64,1}) = spike_dv_ratio(spike,40000)

@compat function get_features(session::AbstractString)
	features = Array(Features,0)
	get_features!(features,session)
	features
end

@compat function get_features!(features::Array{Features,1},session::AbstractString)
	templatefiles = split(readchomp(`find . -name "$(session)*templates*.hdf5"`),"\n")	
	return get_features!(features,templatefiles)
end

@compat function get_features{T<:AbstractString}(templatefiles::Array{T,1})
	features = Array(Features,0)
	get_features!(features,templatefiles)
	features
end

@compat function get_features!{T<:AbstractString}(features::Array{Features,1},templatefiles::Array{T,1})
	w = Float64[]
	isi = Float64[]
	for tf in templatefiles
		m = match(r"([[:alpha:][:digit:]_]*)_templatesg([[:digit:]]*)", tf)
		dd,pp = splitdir(tf)
		session = m.captures[1]
		group = m.captures[2]
		TF = get_valid_templates(tf)
		if TF == nothing
			continue
        end
		if !isempty(dd)
			f = "$(dd)/$(session)g$(group)Spiketrains.mat"
		else
			f = "$(session)g$(group)Spiketrains.mat"
		end
		if !isfile(f)
			continue
		end
		DD = MAT.matread(f)
		if length(DD) == size(TF.templates,3)
			_w = spike_width(TF)
			for c in 1:size(TF.templates,3)
				cc = @sprintf "%02d" c
				_isi = diff(DD["cluster$(cc)s"][:])
				spikes_in_bursts = (sum(_isi .< 5)+1)/length(_isi+1)
				m_isi = percentile(_isi,5)
				dvr = spike_dv_ratio(TF.templates[:,:,c][:])
				push!(features, Features(_w[c], spikes_in_bursts,m_isi,dvr))
			end
		end
	end
	return features 
end
