#TODO add vortex merging

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
dipole(x, y, x0,y0,d,sigma) = (gaussian(x, y, x0+d/2, y0, sigma) - gaussian(x, y, x0-d/2, y0, sigma))
tripole(x,y,x0,y0,r,sigma) = gaussian(x, y, x0+r*cos(0*2pi/3), y0+r*sin(0*2pi/3), sigma) + gaussian(x, y, x0+r*cos(2pi/3), y0+r*sin(2pi/3), sigma) + gaussian(x, y, x0+r*cos(2*2pi/3), y0+r*sin(2*2pi/3), sigma)

function set_dipole_euler(model)
	#Initial Conditions
	reset_state(model.state)

	for i in model.mesh.nh+1:model.mesh.nx-model.mesh.nh, j in model.mesh.nh+1:model.mesh.ny-model.mesh.nh
		x = model.mesh.xc[i,j]
		y = model.mesh.yc[i,j]
		model.state.omega[i,j] = dipole(x, y, 0.5,0.5,0.3,0.05) * model.mesh.msk2d[i,j]
	end
end

function set_vortex_merge_euler(model)
	#Initial Conditions
	reset_state(model.state)

	for i in model.mesh.nh+1:model.mesh.nx-model.mesh.nh, j in model.mesh.nh+1:model.mesh.ny-model.mesh.nh
		x = model.mesh.xc[i,j]
		y = model.mesh.yc[i,j]
		model.state.omega[i,j] = (gaussian(x, y, 0.5,0.5-0.15,0.05) + gaussian(x,y, 0.5, 0.5+0.15, 0.05)) * model.mesh.msk2d[i,j]
	end
end

function set_tripole_euler(model)
	#Initial Conditions
	reset_state(model.state)

	for i in model.mesh.nh+1:model.mesh.nx-model.mesh.nh, j in model.mesh.nh+1:model.mesh.ny-model.mesh.nh
		x = model.mesh.xc[i,j]
		y = model.mesh.yc[i,j]
		model.state.omega[i,j] = tripole(x, y, 0.5,0.5,0.3,0.05) * model.mesh.msk2d[i,j]
	end
end

function set_random_vortices(model)
	reset_state(model.state)
	n_vort = 10
	
	vorts = []
	for i in 1:n_vort
		push!(vorts,(2*rand(Float64)-1, rand(Float64), rand(Float64), 0.1 * rand(Float64)))
	end

	for i in model.mesh.nh+1:model.mesh.nx-model.mesh.nh, j in model.mesh.nh+1:model.mesh.ny-model.mesh.nh
		x = model.mesh.xc[i,j]
		y = model.mesh.yc[i,j]
		
		for vort in vorts
			model.state.omega[i,j] += vort[1] * gaussian(x, y, vort[2], vort[3], vort[4]) * model.mesh.msk2d[i,j]
		end
	end
end	

function get_eulerpsi(step_func, interp_func)
	#Defining our equation
	@Let omega = FormVariable{2, Dual}() #Vorticity

	@Let psi = InverseLaplacian(omega) #∇²ω = Ψ
	@Let u = Codifferential(psi) #u = δΨ
	@Let U = Sharp(u)

	#Time derivative
	@Let dtomega = - ExteriorDerivative(InteriorProduct(U, omega)) #dtω = L(U,ω)

	#Defining the parameters needed to explicit
	explparams = ExplicitParam(; interp = interp_func)

	#Generating the RHS
	euler_rhs! = to_kernel(dtomega; save = ["u_x", "u_y", "ι_U_omega"], explparams = explparams, verbose = false)

	#Testing the function

	#Defining the Mesh
	nx = 50
	ny = 50
	nh = 3

	msk = zeros(nx, ny)
	msk[nh+1:nx-nh, nh+1:ny-nh] .= 1
	#msk[nx÷2-nx÷5:nx÷2+nx÷5, 2*ny÷10:4*ny÷10] .= 0

	Lx, Ly = (1,1)
		
	#LoopManager
	scalar = PlainCPU()
	simd = VectorizedCPU(16)
	threads = MultiThread(scalar)
	
	mesh = Arrays.Mesh(nx, ny, nh, simd, msk, Lx, Ly)

	#Creating the State
	state = State(mesh)
	#Creating the Model
	model = Model(euler_rhs!, mesh, state, ["omega"]; cfl = 100., dtmax = 5., integratorstep! = step_func)
	
	println("Precompiling model...")
	Base.invokelatest(step!,model)
	println("Done !")
	reset_state(model.state)

	return model, :omega, Dict(("Dipole"=>set_dipole_euler, 
				    "Tripole"=>set_tripole_euler, 
				   "Soup"=>set_random_vortices,
				   "VortexMerge"=>set_vortex_merge_euler))
end
