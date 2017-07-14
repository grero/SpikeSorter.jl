module SpikeSorter
using Compat
using DSP
using MultivariateStats
using StatsBase
import StatsBase.transform, StatsBase.fit

include("types.jl")
include("utility.jl")
include("hmm.jl")
include("viterbi.jl")
include("readfiles.jl")
include("features.jl")
include("extraction.jl")
include("io.jl")

#import GUICheck
#if GUICheck.hasgui()
#	include("plot.jl")
#end

export SpikeWaveforms

end # module
