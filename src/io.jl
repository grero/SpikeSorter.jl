import Base.write

function Base.write(io::IO, waveforms::SpikeWaveforms)
    #write a header
    headersize = 64
    nbytes = 0
    npoints, nchannels, nspikes = size(waveforms.waveforms) 
    nbytes += write(io, npoints) 
    nbytes += write(io, nchannels)
    nbytes += write(io, nspikes)
    seek(io, headersize)

    #write timestamps, followed by waveforms
    nbytes += write(io, waveforms.timestamps)
    nbytes += write(io, waveforms.waveforms)
    nbytes
end

function Base.read(io::IO, ::Type{SpikeWaveforms})
    #read the header
    headersize = 64
    npoints = read(io, Int64)
    nchannels = read(io, Int64)
    nspikes = read(io, Int64)
    seek(io, headersize)

    #read timestamps
    timestamps = read(io, Float64, nspikes)
    waveforms = read(io, Float64, (npoints, nchannels, nspikes))
    SpikeWaveforms(waveforms, timestamps)
end

