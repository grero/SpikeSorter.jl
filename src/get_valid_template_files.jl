#!/usr/bin/env julia
using ArgParse
import SpikeSorter
#set up arguments
s = ArgParseSettings()
@add_arg_table s begin
	"--session"
		help = "The name of the session to be checked"
	"--channel"
		help = "The channel to be processed"
		arg_type = Int
	"--chunk"
		help = "The chunk for process"
		arg_type = Int
	"--dir"
		help = "Directory to run from"
	"files"
		action = :store_arg
		nargs = '*'
		arg_type = String
		help = "Files to be processed"
end

D = parse_args(ARGS,s)
if isempty(D["files"])
	ss = D["session"]
	if !isempty(D["channel"])
		ch  = @sprintf "%04d" D["channel"]
	else
		ch = "*"
	end
	if !isempty(D["chunk"])
		chunk = @sprintf "%04d" D["chunk"]
	else
		chunk = "*"
	end
	cwd = get(D, "dir", ".")
	_files = readchomp(`find $cwd -name "$(ss)_templatesg$(ch).$(chunk).hdf5" -depth 1`)
	if isempty(_files)
		exit(1)
	end
	files = split(_files, "\n")

else
	#files = convert(Array{ASCIIString,1},D["files"])
	files = D["files"]
end
gfiles = filter(SpikeSorter.is_valid_template_file,files)
if isempty(gfiles)
	exit(1)
end
println(join(gfiles,"\n"))
exit(0)
