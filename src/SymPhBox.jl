module SymPhBox

using SymPh, SymPh.Maths 
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread
using Plots
using Gtk4
using Printf

include("models.jl")

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
	 "model_func" => nothing,
	 "step_func" => nothing,
	 "interp_func" => nothing,
	 "var" => nothing, #plotted + edited var
	 "ics" => nothing,
	 "icid" => nothing,

	 #Display params
	 "palette" => :balance,
	 ))

#FROM GR.jl examples (gtk4_plots_ex)
function _plot(ctx, w, h)
	ENV["GKS_WSTYPE"] = "142"
	ENV["GKSconid"] = @sprintf("%lu", UInt64(ctx.ptr))

	gr(show=true)
	model = params["model"]
	mesh = model.mesh
	heatmap(getproperty(model.state, params["var"])[mesh.nh+1:mesh.nx-mesh.nh, mesh.nh+1:mesh.ny-mesh.nh], 
		size = (w,h), 
		c=params["palette"],
		aspect_ratio = 1)
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
	dd = b["ic_dropdown"]
	
	if length(Gtk4.model(dd)) >0
		ic = selected_string(dd)
	
		params["ics"][ic](params["model"])
		
		draw(params["canvas"])
	end
end

function change_model(widget, thing)
	#Model
	model_string = selected_string(b["model_dropdown"])
	if model_string == "EulerPsi"
		params["model_func"] = get_eulerpsi
	elseif model_string == "RSW"
		params["model_func"] = get_rsw
	else
		println("Not Implemented")
	end
end

function change_integrator(widget, thing)
	#Integrator
	integ_string = selected_string(b["integrator_dropdown"])
	if integ_string == "RK3"
		params["step_func"] = rk3step!
	elseif integ_string == "EulerForward"
		params["step_func"] = euler_forwardstep!
	else
		println("Not Implemented")
	end
end

function change_interpolation(widget, thing)
	#Interpolation
	interp_string = selected_string(b["interp_dropdown"])
	if interp_string == "Weno5"
		params["interp_func"] = Arrays.weno
	elseif interp_string == "Upwind5"
		params["interp_func"] = Arrays.upwind
	else
		println("Not Implemented")
	end
end

function change_palette(widget, thing)
	palette_string = selected_string(b["cmap_dropdown"])
	if palette_string == "Balance"
		params["palette"] = :balance
	elseif palette_string == "Viridis"
		params["palette"] = :viridis
	elseif palette_string == "Magma"
		params["palette"] = :magma
	else
		println("Not Implemented")
	end

	draw(params["canvas"])
end

function populate_ics(ics)	
	#Initial conditions
	params["ics"] = ics
	strlist = Gtk4.model(b["ic_dropdown"])
	if params["icid"] != nothing
		signal_handler_disconnect(b["ic_dropdown"], params["icid"])
	end

	empty!(strlist)
	for key in keys(ics)
		push!(strlist, key)
	end
	params["icid"] = signal_connect(change_ics, ic_dropdown, "notify::selected")
	selected_string!(b["ic_dropdown"], collect(keys(ics))[1])
end

function get_mesh()
	nx = parse(Int, b["nx_entry"].text)
	ny = parse(Int, b["ny_entry"].text)
	nh = 3
	Lx, Ly = (1,1)
		
	msk = zeros(nx, ny)
	msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

	simd = VectorizedCPU(16)

	return Arrays.Mesh(nx, ny, nh, simd, msk, Lx, Ly)
end

#Creating a model and setting it up TODO maybe only generate model once clicked on a validation button
function generate_model(widget)
	stop(widget)
	#Creating the model
	mesh = get_mesh()
	model, var, ics = params["model_func"](params["step_func"], params["interp_func"], mesh)
	params["model"] = model
	params["var"] = var

	populate_ics(ics)
	change_ics(widget, nothing)

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
main_container.child = canvas

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

signal_connect(change_model, model_dropdown, "notify::selected")
signal_connect(change_integrator, integrator_dropdown, "notify::selected")
signal_connect(change_interpolation, interp_dropdown, "notify::selected")

ic_dropdown = b["ic_dropdown"]
params["icid"] = signal_connect(change_ics, ic_dropdown, "notify::selected")

gen_model_button = b["gen_model_button"]
signal_connect(generate_model, gen_model_button, "clicked")

cmap_dropdown = b["cmap_dropdown"]
signal_connect(change_palette, cmap_dropdown, "notify::selected")

nx_entry = b["nx_entry"]
ny_entry = b["ny_entry"]
#signal_connect(change_mesh, nx_entry, "changed")
#signal_connect(change_mesh, ny_entry, "changed")

#Creating the initial model
change_model(nothing, nothing)
change_integrator(nothing, nothing)
change_interpolation(nothing, nothing)
generate_model(nothing)

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
