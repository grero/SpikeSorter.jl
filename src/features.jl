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
