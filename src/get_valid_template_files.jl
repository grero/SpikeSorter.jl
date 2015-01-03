#!/usr/bin/env julia
import SpikeSorter
gfiles = filter(SpikeSorter.is_valid_template_file,ARGS)
println(join(gfiles,"\n"))

