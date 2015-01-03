module SpikeSorter
include("types.jl")
include("utility.jl")
include("hmm.jl")
include("viterbi.jl")
include("readfiles.jl")
import GUICheck
if GUICheck.hasgui()
	include("plot.jl")
end


end # module
