import Winston
import Visualizer.visualize
import HDF5
import Color

function visualize{T<:String}(fnames::Array{T,1})
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

function plot(fname::String)
	p = Winston.FramedPlot()
	Winston.setattr(p.x2, "draw_axis", false)
	Winston.setattr(p.y2, "draw_axis", false)
	Winston.setattr(p.y1, "tickdir", +1)
	Winston.setattr(p.x1, "tickdir", +1)
	plot(p,fname)
	p
end

function plot(p::Winston.FramedPlot,fname::String)
	tf = TemplateFile(fname)	
	plot(p, tf)
	p
end

function plot(p::Winston.FramedPlot, tf::TemplateFile)
	if tf == nothing
		return p
	end
	colors = Color.distinguishable_colors(tf.ntemplates)
	t = 1:size(tf.templates,1)
	for i=1:tf.ntemplates
		Winston.add(p, Winston.Curve(t,tf.templates[:,:,i][:];color=colors[i]))
	end
	p
end


