#TODO add vortex merging

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
dipole(x, y, x0,y0,d,sigma) = (gaussian(x, y, x0+d/2, y0, sigma) - gaussian(x, y, x0-d/2, y0, sigma))
tripole(x,y,x0,y0,r,sigma) = gaussian(x, y, x0+r*cos(0*2pi/3), y0+r*sin(0*2pi/3), sigma) + gaussian(x, y, x0+r*cos(2pi/3), y0+r*sin(2pi/3), sigma) + gaussian(x, y, x0+r*cos(2*2pi/3), y0+r*sin(2*2pi/3), sigma)

function set_dipole(var, model)
	#Initial Conditions
	reset_state(model.state)

	assvar = getproperty(model.state, var)
	for i in model.mesh.nh+1:model.mesh.nx-model.mesh.nh, j in model.mesh.nh+1:model.mesh.ny-model.mesh.nh
		x = model.mesh.xc[i,j]
		y = model.mesh.yc[i,j]
		assvar[i,j] = dipole(x, y, 0.5,0.5,0.3,0.05) * model.mesh.msk2d[i,j]
	end
end

function set_tripole(var, model)
	#Initial Conditions
	reset_state(model.state)

	assvar = getproperty(model.state, var)
	for i in model.mesh.nh+1:model.mesh.nx-model.mesh.nh, j in model.mesh.nh+1:model.mesh.ny-model.mesh.nh
		x = model.mesh.xc[i,j]
		y = model.mesh.yc[i,j]
		assvar[i,j] = tripole(x, y, 0.5,0.5,0.3,0.05) * model.mesh.msk2d[i,j]
	end
end
