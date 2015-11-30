import Winston
#import Visualizer.visualize
import HDF5
import Colors

@compat function visualize{T<:AbstractString}(fnames::Array{T,1})
	p = Winston.FramedPlot()
	Winston.setattr(p.x2, "draw_axis", false)
	Winston.setattr(p.y2, "draw_axis", false)
	Winston.setattr(p.y1, "tickdir", +1)
	Winston.setattr(p.x1, "tickdir", +1)
	function func3(p,fnames,i)
		if 0 < i <= length(fnames)
			tf = TemplateFile(fnames[i])
			if tf != nothing
				plot(p, tf)
				_title = replace(fnames[i], "_", "\\_")
				Winston.setattr(p, "title", _title)
			end
		end
	end
	visualize(fnames,800,600,"Spike forms", func3, p)
end

@compat function plot(fname::AbstractString)
	p = Winston.FramedPlot()
	Winston.setattr(p.x2, "draw_axis", false)
	Winston.setattr(p.y2, "draw_axis", false)
	Winston.setattr(p.y1, "tickdir", +1)
	Winston.setattr(p.x1, "tickdir", +1)
	plot(p,fname)
	p
end

@compat function plot(p::Winston.FramedPlot,fname::AbstractString)
	tf = TemplateFile(fname)	
	plot(p, tf)
	p
end

function plot(p::Winston.FramedPlot, tf::TemplateFile)
	if tf == nothing
		return p
	end
	colors = Colors.distinguishable_colors(tf.ntemplates)
	t = 1:size(tf.templates,1)
	for i=1:tf.ntemplates
		Winston.add(p, Winston.Curve(t,tf.templates[:,:,i][:];color=colors[i]))
	end
	p
end

function plot(F::Array{Features,1};kvs...)
	p = Winston.FramedPlot()
	plot(p, F;kvs...)
	p
end

function plot(p::Winston.FramedPlot, F::Array{Features,1};feature1::Symbol=:spike_width, feature2::Symbol=:spikes_in_bursts)
	f1 = zeros(length(F))
	f2 = zeros(length(F))
	if !(feature1 in names(F[1]) && feature2 in names(F[1]))
		ArgumentError("Invalid feature name")
	end
	for i in 1:length(f1)
		f1[i] = getfield(F[i],feature1)
		f2[i] = getfield(F[i],feature2)
	end
	#throw away outliers
	pf1 = percentile(f1, [1,99])	
	pf2 = percentile(f2, [1,99])	
	fidx = (f1.>pf1[1])&(f1.<pf1[2])&(f2.>pf2[1])&(f2.<pf2[2])
	f1 = f1[fidx]
	f2 = f2[fidx]

	Winston.add(p, Winston.Points(f1,f2))
	Winston.setattr(p.x2, "draw_axis",false)
	Winston.setattr(p.y2, "draw_axis",false)
	Winston.setattr(p.x1, "tickdir", 1)
	Winston.setattr(p.y1, "tickdir", 1)
	Winston.setattr(p, "xlabel", replace(string(feature1),"_", " "))
	Winston.setattr(p, "ylabel", replace(string(feature2),"_", " "))
	p
end

