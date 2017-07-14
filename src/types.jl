import HDF5

type SpikeWaveforms
    waveforms::Array{Float64,3}
    timestamps::Array{Float64,1}
end

abstract SpikeSorterResult
abstract FeatureModel

type PCAFeatureModel <: FeatureModel
    features::Array{Float64,2}
    model::PCA
end

transform(m::PCAFeatureModel, y) = transform(m.model, y)

type FeatureSorter{T1 <: FeatureModel, T2} <: SpikeSorterResult
    waveforms::SpikeWaveforms
    featuremodel::T1
    clusterid::Array{Int64,1}
    clustermodel::T2
end

function SpikeWaveforms(waveforms::Array{Float64,2}, timestamps::Array{Float64,1})
    npoints,nspikes = size(waveforms)
    _waveforms = zeros(npoints, 1, nspikes)
    _waveforms[:,1,:] = waveforms
    SpikeWaveforms(_waveforms, timestamps)
end

type TemplateFile
	templates::Array{Float64,3}
	cinv::Array{Float64,2}
	ntemplates::Int64
	nchannels::Int64
end

type SortedFile
	templates::Array{Float64,3}
	cinv::Array{Float64,2}
	ntemplates::Int64
	nchannels::Int64
	mlseq::Array{UInt8,2}
end

type Features
	spike_width::Float64
	spikes_in_bursts::Float64
	low_isi::Float64
	dv_ratio::Float64
end

@compat function TemplateFile(fname::AbstractString;verbose::Integer=0)
	if !HDF5.ishdf5(fname)
		verbose > 0 && println("$fname is not a valid hdf5 file")
		return nothing
	end
	ff = HDF5.h5open(fname) 
	if !("spikeForms" in names(ff))
		verbose > 0 && println("File does not contain any spike forms")
		close(ff)
		return nothing
	end
	templates = HDF5.read(ff, "spikeForms")
	nspikes = size(templates,3)
	nchannels = size(templates,2)
	if nchannels > 1
		cinv = HDF5.read(ff, "cinv")
	else
		cinv = zeros(1,1)
		cinv[1,1] = HDF5.read(ff, "cinv")[:][1]
	end
	close(ff)
	return TemplateFile(templates, cinv,nspikes,nchannels)
end

