
function set_vortex_rsw(model)
	reset_state(model.state)

	h0 = 0.05
	H = 1
	sigma = 0.05
	gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
	
	mesh = model.mesh
	state = model.state

	h = state.h

	for i in mesh.nh+1:mesh.nx-mesh.nh, j in mesh.nh+1:mesh.ny-mesh.nh
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]

		h[i,j] = (H + h0 * gaussian(x, y, 0.5, 0.5, sigma)) * mesh.A[5,5]
	end

	state.f .=  10 .* ones((mesh.nx,mesh.ny)) .* mesh.A .* mesh.msk2d
end

function set_dipole_rsw(model)
	reset_state(model.state)
	
	h0 = 0.05
	H = 1
	sigma = 0.05
	gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
	d=0.05

	mesh = model.mesh
	state = model.state

	h = state.h

	for i in mesh.nh+1:mesh.nx-mesh.nh, j in mesh.nh+1:mesh.ny-mesh.nh
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]

		h[i,j] = (H + h0 * (gaussian(x, y, 0.5+d/2, 0.5, sigma) - gaussian(x, y, 0.5-d/2, 0.5, sigma))) * mesh.A[5,5]
	end

	state.f .=  10 .* ones((mesh.nx,mesh.ny)) .* mesh.A .* mesh.msk2d
end

function get_rsw(step_func, interp_func, mesh)
	#Defining our equation
	@Let h = FormVariable{2, Primal}() #Height * A (h* technically)
	@Let u = FormVariable{1, Dual}() #Transported velocity

	@Let U = Sharp(u) # U = u#
	@Let k = 0.5 * Hodge(InnerProduct(u,u)) #k = 0.5 * hodge(innerproduct(u,u))
	@Let p = Hodge(h) # p = *(g(h*+b*))
	@Let zeta = ExteriorDerivative(u) # ζ* = du
	@Let f = FormVariable{2, Dual}() #Coriolis (f* so times A)
	@Let pv = (f + zeta) / h #TODO check what pv should be

	#Time derivative
	@Let dtu = -InteriorProduct(U, zeta + f) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
	@Let dth = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))

	#Defining the parameters needed to explicit
	explparams = ExplicitParam(; interp = interp_func)

	#Generating the RHS
	rsw_rhs! = to_kernel(dtu, dth, pv; save = ["zeta", "k"], explparams = explparams)

	#Initial Conditions
	state = State(mesh)

	#Creating the Model
	model = Model(rsw_rhs!, mesh, state, ["u_x", "u_y", "h"]; integratorstep! = step_func, cfl = 0.15, dtmax=0.15)
	
	println("Precompiling model...")
	Base.invokelatest(step!,model)
	println("Done !")
	reset_state(model.state)
	
	return model, :h, Dict(("Dipole"=>set_dipole_rsw, "Vortex"=>set_vortex_rsw))
end
