# RDF Algorithm adapted from book: 
#       Understanding Molecular Simulation, From Algorithms to Applications 3rd edition

# Functions 

"""
Get nearest image of a particle
"""
function get_nearest_img(posi, posj, boxlength)
	dist = posi - posj
	return dist - boxlength*round(Int, dist/boxlength)
end

"""
Accumulate histogram of pair distances 
"""
function grsample(particles, g)
	# ngr += 1 # incerase counter

	# Loop over all particle pairs
	for i in range(1,length=nparts-1)	
		for j in range(i+1,length=nparts)
			xr = get_nearest_img(particles[i][2], particles[j][2], box[1])
			yr = get_nearest_img(particles[i][3], particles[j][3], box[2])
			zr = get_nearest_img(particles[i][4], particles[j][4], box[3])
			r = sqrt(xr^2 + yr^2 + zr^2)
			if r < box/2 # only consider dist less than box/2
				ig = Int(r/delg)
				g[ig] += 2 # increment histogram for both particles
			end
		end
	end
end

"""
Normalises the radial distribution at end of run
"""
function grnormalise()
	gfac = (4/3)*pi*(delg^3) # convert bins to 3d shells
	for i in range(1,length=nhis)
		vb = gfac*((i + 1)^3 - i^3) # volume of i-th bin
		nid = vb*rho # number of ideal gas particles in volume vb
		g(i) = g(i)/(ngr*npart*nid) # normalise g(r)
	end
end

# User-defined vars

nhis = 100 # number of bins
box = (50, 50, 50) # box size (x, y, z)
#particles = # vector of tuples (particleType, X, Y, Z) 

# Other vars

delg = box/(2*nhis) # bin size
ngr = 0 # count the number of times grsample() is called
g = zeros(Float64, nhis) # vector to store the g(r)
#nparts = 
rho = nparts/(box[1]*box[2]*box[3]) # number density of system

println(g)
