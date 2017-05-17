module SpikeSorter
using Compat
using DSP
include("types.jl")
include("utility.jl")
include("hmm.jl")
include("viterbi.jl")
include("readfiles.jl")
include("features.jl")
include("extraction.jl")
import GUICheck
if GUICheck.hasgui()
	include("plot.jl")
end


end # module
