# RDF Algorithm adapted from book: 
#       Understanding Molecular Simulation, From Algorithms to Applications 3rd edition

# Imports

#using TarIterators

# Functions 

"""
Get nearest image of a particle
"""
function get_nearest_img(posi::Float64, posj::Float64, boxlength::Float64)
	dist = posi - posj
	return dist - boxlength*round(Int, dist/boxlength)
end

"""
Accumulate histogram of pair distances 
"""
function grsample(particles, g, nparts::Int64, delg::Float64, rmax::Float64, boxx::Float64, boxy::Float64, boxz::Float64)
	# ngr += 1 # incerase counter
	
	# Loop over all particle pairs
	Threads.@threads for i in range(1,stop=nparts-1)	
		Threads.@threads for j in range(i+1,stop=nparts)
			xr = get_nearest_img(parse(Float64,particles[i][2]), parse(Float64,particles[j][2]), boxx)
			yr = get_nearest_img(parse(Float64,particles[i][3]), parse(Float64,particles[j][3]), boxy)
			zr = get_nearest_img(parse(Float64,particles[i][4]), parse(Float64,particles[j][4]), boxz)
			r = sqrt(xr^2 + yr^2 + zr^2)
			if r < rmax #box/2 # only consider dist less than box/2
				ig = trunc(Int,r/delg)
				g[ig] += 2 # increment histogram for both particles
			end
		end
	end
end

"""
Normalises the radial distribution at end of run
"""
function grnormalise(g, delg::Float64, nhis::Int64, rho::Float64, ngr::Int64, npart::Int64)
	gfac = (4/3)*pi*(delg^3) # convert bins to 3d shells
	for i in range(1,stop=nhis)
		vb = gfac*((i + 1)^3 - i^3) # volume of i-th bin
		nid = vb*rho # number of ideal gas particles in volume vb
		g[i] = g[i]/(ngr*npart*nid) # normalise g(r)
	end
end

# User-defined vars

nhis = 100 # number of bins
box = (50, 50, 50) # box size (x, y, z)
#particles = # vector of tuples (particleType, X, Y, Z) 

# Other vars

#delg = box/(2*nhis) # bin size
ngr = 0 # count the number of times grsample() is called
g = zeros(Float64, nhis) # vector to store the g(r)
#nparts = 
#rho = nparts/(box[1]*box[2]*box[3]) # number density of system

# Read the contents of the LAMMPS traj file
# File must have columns in order: id element xu yu zu

open("test") do f
	nparts = 0
	boxx = 0
	boxy = 0
	boxz = 0
	delg = 0.2
	ngr = 0
	nhis = 100
	rho = 0

	rmax = nhis*delg
	g = zeros(Float64, nhis) # vector to store the g(r)

	while ! eof(f)
		line = readline(f) # read a line from the file
		
		# Extract number of atoms -> nparts
		if startswith(line, "ITEM: NUMBER OF ATOMS")
			nparts = parse(Int64,readline(f))
		# Extract the 	size of the box
		elseif startswith(line, "ITEM: BOX BOUNDS xy xz yz pp pp pp")	
			# xlo xhi xy
			# ylo yhi xz
			# zlo zhi yz
			line = readline(f)
			lo, hi, _ = split(line, " ")
			boxx = parse(Float64, hi) - parse(Float64, lo)
			line = readline(f)
			lo, hi, _ = split(line, " ")
			boxy = parse(Float64, hi) - parse(Float64, lo)
			line = readline(f)
			lo, hi, _ = split(line, " ")
			boxz = parse(Float64, hi) - parse(Float64, lo)
			
			rho = nparts/(boxx*boxy*boxz) # number density of system
		# Extract particle positions
		elseif startswith(line, "ITEM: ATOMS id element xu yu zu")
			particles = Vector{NTuple{4, String}}(undef, nparts)
			line = readline(f)
			count = 1
			while ! startswith(line,"ITEM:")
				_, elem, xu, yu, zu = split(line, " ")
				particles[count] = (elem, xu, yu, zu)

				count += 1
				line = readline(f)
			end

			# Run rdf algorithm after reading all particles of one timestep
			ngr += 1
			print("Analysing step $ngr... ")
			grsample(particles, g, nparts, delg, rmax, boxx, boxy, boxz)
			println("Done")
			break
		end
	end

	grnormalise(g, delg, nhis, rho, ngr, nparts) # normalise the rdf
	for i in range(1,stop=nhis)
		println("$(i*delg) $(g[i])")
	end
end
