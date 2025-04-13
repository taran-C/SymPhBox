module SymPhBox

using SymPh, SymPh.Maths 
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread
using Plots
using Gtk4
using Printf

include("models.jl")
include("initial_conditions.jl")

#Building the window
ENV["JULIA_NUM_THREADS"] = (2,1)
b = GtkBuilder("interface/symphbox.ui")
win = b["main"]
main_container = b["main_container"]

params = Dict{String, Any}((
	 "stepping" => false,
	 "playing" => false,
	 "canvas" => nothing,

	 #Current state of run config
	 "model" => nothing,
	 "var" => nothing, #plotted + edited var
	 "ics" => nothing,

	 #Display params
	 "palette" => :phase,#:balance,
	 ))

#FROM GR.jl examples (gtk4_plots_ex)
function _plot(ctx, w, h)
	ENV["GKS_WSTYPE"] = "142"
	ENV["GKSconid"] = @sprintf("%lu", UInt64(ctx.ptr))

	gr(show=true)
	heatmap(getproperty(params["model"].state, params["var"]), size = (w,h), c=params["palette"])
end

function do_step(widget)
	if !params["stepping"]
		params["stepping"] = true
		Threads.@spawn begin
			step!(params["model"]; n=10)
			draw(params["canvas"])
			params["stepping"] = false
		end
	end
end

#TODO send a signal when finished emitting to allow waiting safely and not crashing
function play(widget)
	if !params["playing"]
		params["playing"] = true
		params["stepping"] = true
		Threads.@spawn begin
			while params["playing"]
				step!(params["model"]; n=3)
				draw(params["canvas"])
				sleep(1/60) #TODO adjust framerate
			end
		end
	end

end

function stop(widget)
	params["playing"] = false
	params["stepping"] = false
end

function change_ics(widget, thing)
	stop(widget)

	ic = selected_string(widget)

	if ic != params["ics"]
		
		params["ics"] = ic

		if ic == "Dipole"
			set_dipole(params["var"], params["model"])
		elseif ic == "Tripole"
			set_tripole(params["var"], params["model"])
		else
			println("Not implemented")
		end
		
		draw(params["canvas"])
	end
end

#Creating a model and setting it up TODO maybe only generate model once clicked on a validation button
function create_model(widget, thing)
	stop(widget)

	#Model
	model_string = selected_string(b["model_dropdown"])
	if model_string == "EulerPsi"
		model_func = get_eulerpsi
	else
		println("Not Implemented")
	end

	#Integrator
	integ_string = selected_string(b["integrator_dropdown"])
	if integ_string == "RK3"
		step_func = rk3step!
	elseif integ_string == "EulerForward"
		step_func = euler_forwardstep!
	else
		println("Not Implemented")
	end

	#Interpolation
	interp_string = selected_string(b["interp_dropdown"])
	if interp_string == "Weno5"
		interp_func = Arrays.weno
	elseif interp_string == "Upwind5"
		interp_func = Arrays.upwind
	else
		println("Not Implemented")
	end

	#Creating the model
	model, var = model_func(step_func, interp_func)
	params["model"] = model
	params["var"] = var
	
	#Initial conditions
	params["ics"] = nothing
	change_ics(b["ic_dropdown"], thing)

	#(re)draw canvas
	draw(params["canvas"])
end

#Defining and adding drawing area
canvas = GtkCanvas()
params["canvas"] = canvas
@guarded draw(canvas) do widget
	w, h = width(canvas), height(canvas)
	#Create a Cairo context for drawing
	ctx = getgc(canvas)

	_plot(ctx, w, h)
end
main_container[2] = canvas

#Control Buttons
step_button = b["step_button"]
play_button = b["play_button"]
stop_button = b["stop_button"]

signal_connect(do_step, step_button, "clicked")
signal_connect(play, play_button, "clicked")
signal_connect(stop, stop_button, "clicked")

#Configuration Controls
model_dropdown = b["model_dropdown"]
integrator_dropdown = b["integrator_dropdown"]
interp_dropdown = b["interp_dropdown"]
signal_connect(create_model, model_dropdown, "notify::selected")
signal_connect(create_model, integrator_dropdown, "notify::selected")
signal_connect(create_model, interp_dropdown, "notify::selected")

ic_dropdown = b["ic_dropdown"]
signal_connect(change_ics, ic_dropdown, "notify::selected")

#Creating the initial model
create_model(nothing, nothing)

#Displaying everything
show(win)

#Allow non REPL usage, slow for X reason
if !isinteractive()
	cond = Condition()
	signal_connect(win, :destroy) do widget
		notify(cond)
	end

	@async Gtk4.GLib.glib_main()

	wait(cond)
end

end
