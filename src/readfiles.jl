import MAT

function loadWaveformsFile(fname::String;time_conversion::Real=0.001)
	fid = open(fname,"r")
	waveforms = Array(Int16,0)
	timestamps = Array(Uint64,0)
	try
		headersize = int(read(fid,Uint32))
		num_spikes = int(read(fid,Uint32))
		num_chs = int(read(fid,Uint8))
		sampling_rate = int(read(fid,Uint32))
		timepts = int(read(fid,Uint32))
		seek(fid,headersize)
		append!(waveforms,read(fid,Int16,num_chs*timepts*num_spikes))
		waveforms = reshape(waveforms,(num_chs,timepts,num_spikes))
		seek(fid,headersize + num_spikes*num_chs*timepts*sizeof(Int16))
		append!(timestamps,read(fid,Uint64,num_spikes))
	finally
		close(fid)
	end
	return float(waveforms),float(timestamps)*time_conversion
end

function getspiketrains(;session::String="",groups::Array{Int64,1}=Array(Int64,0),checkArtifact::Bool=true,verbose::Integer=0)
	if isempty(groups)
		#get the groups from the waveforms files in the current directory
		files = split(chomp(readall(`find . -name "*waveforms.bin"`)))
	else
		files = []
	end
	SS = pmap(f->begin
					m = match(r"([A-Za-z0-9_]*)g([0-9]*)",f)
					session = m.captures[1]
					group = int(m.captures[2])
					sptrains = getspiketrains(session,group;checkArtifact=checkArtifact,verbose=verbose)
					sptrains
			end,files)
	return merge(SS...)
end

function getspiketrains(session::String, group::Integer;checkArtifact::Bool=true,verbose::Integer=0)
	verbose > 0 && println("Processing group $group ...")
	wfname = @sprintf "%sg%04dwaveforms.bin" session group
	overlap_fname = @sprintf "%sg%04dwaveforms.overlap" session group
	cut_fname = @sprintf "%sg%04dwaveforms.cut" session group

	waveforms,timestamps = loadWaveformsFile(wfname)
	D = MAT.matread(overlap_fname)
	cids = map(int,readlines(open(cut_fname,"r")))
	clusters = unique(cids)
	clusters = clusters[clusters.>0]
	sptrains = Dict{String,Array{Float64,1}}()
	c2 = 1
	for c in clusters
		spidx = (cids.==c) | ((cids.==-1)&(D["overlaps"][:,c].==1))
		cidx = cids.==c
		if checkArtifact
			mm = mean(waveforms[:,:,bool(cidx[:])],3)
			ii = indmax(-mm) #get the minimum point
			goodcluster = ii in div(size(mm,2),3):div(2*size(mm,2),3)
			goodcluster = goodcluster && -minimum(mm) >= maximum(mm)
			if goodcluster
				cell = "g$(group)c$c2"
				sptrains[cell] = timestamps[spidx[:]]
				c2 += 1
			else
				verbose > 0 && println("\tCluster $c classified as artifact and removed")
			end
		else
			cell = "g$(group)c$c"
			sptrains[cell] = timestamps[spidx[:]]
		end
	end
	return sptrains
end
