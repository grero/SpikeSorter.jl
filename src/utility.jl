import HDF5
import Base.isvalid

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

function is_TemplateFile(fname::String)
	if !HDF5.ishdf5(fname)	
		return false
	end
	ff = HDF5.h5open(fname) 
	_names = names(ff)
	close(ff)
	if !("spikeForms" in _names)
		return false
	end
	return true
end

function is_valid_template_file(fname::String)
	return get_valid_templates(fname) != nothing
end

function get_valid_templates(fname::String)
	tf = TemplateFile(fname)
	if tf == nothing
		return nothing
	end
	return get_valid_templates(tf)
end

function get_valid_templates(tf::TemplateFile)
	vidx = zeros(Bool,tf.ntemplates)
	for i=1:tf.ntemplates
		mm = tf.templates[:,:,i]
		ii = indmax(-mm) #get the minimum point
		vidx[i]= ii in div(size(mm,1),3):div(2*size(mm,1),3) #minimum should occur in the middle 3rd of the spike
		vidx[i] = vidx[i] && -minimum(mm) >= maximum(mm) 
	end
	if any(vidx)
		return TemplateFile(tf.templates[:,:,vidx], tf.cinv, sum(vidx), tf.nchannels)
	else
		return nothing
	end
end


