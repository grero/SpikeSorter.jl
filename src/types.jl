import HDF5

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
	mlseq::Array{Uint8,2}
end

type Features
	spike_width::Float64
	spikes_in_bursts::Float64
	low_isi::Float64
end

function TemplateFile(fname::String;verbose::Integer=0)
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

