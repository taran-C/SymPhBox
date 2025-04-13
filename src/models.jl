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

	return model, :omega
end


