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

function get_features(session::String)
	templatefiles = split(readchomp(`find . -name "$(session)*templates*.hdf5"`),"\n")	
	return get_features(templatefiles)
end

function get_features{T<:String}(templatefiles::Array{T,1})
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
			append!(w, spike_width(TF))
			for c in 1:size(TF.templates,3)
				cc = @sprintf "%02d" c
				push!(isi, percentile(diff(DD["cluster$(cc)s"][:]),5))
			end
		end
	end
	return w, isi
end
