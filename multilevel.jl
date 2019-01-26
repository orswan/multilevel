# multilevel.jl
# Evaluates the dynamics of a multilevel atom in the presence of lasers and magnetic fields.
# We use units with hbar=c=1
module multilevel
using LinearAlgebra#,DifferentialEquations

h = 1
c = 1
hbar = h/2pi

struct Level
	name::String
	S::Number				# Spin angular momentum
	L::Number				# Orbital...
	J::Number				# Total electronic...
	I::Number				# Nuclear...
	F::Number				# Nuclear + electronic...
	energy::Number			# Energy (can be w.r.t. any fixed zero point)
	E::Number				# Alias for energy
	Level(name,s,l,j,i,f,e) = new(name,s,l,j,i,f,e,e)
end

Level(name,s,l,j,i,f) = Level(name,s,l,j,i,f,nothing)		# Energy need not be specified
Level(name,s,l,j;energy) = Level(name,s,l,j,0,j,energy)	# Nuclear moment need not be specified

struct Atom
	levels::Array{Level,1}				# Levels
	# Linewidths between each pair of levels:
		# linewidths[i,j] is the ANGULAR linewidth from i to j.
		# So if i is lower energy than j, it is 0.
	linewidths::Array{Float64,2}
	lambda::Array{Float64,2}				# Vacuum wavelength connecting any pair of states
end

function Atom(levels::Array{Level,1},linewidths::Array{<:Number,2})
	# Generates an Atom object using the energy of the levels.  Fails if not all energies are specified. 
	energies = [L.energy for L in levels]
	if any(x->x==nothing,energies)			# Check that all energies are defined
		throw(ArgumentError("Not all Level energies are defined."))
	else
		return Atom(levels,
					convert(Array{Float64,2},linewidths),
					convert(Array{Float64,2},[h*c/abs(L1.energy-L2.energy) for L1 in levels, L2 in levels])
					)
	end
end

function Atom(A::Atom,levels::Array{Level,1})
	# Generates a new Atom from a subset of the levels of an existing atom.
	idxs = [findfirst(x->x==L,A.levels) for L in levels]
	return Atom(levels,A.linewidths[idxs,idxs],A.lambda[idxs,idxs])
end

function Base.show(io::IO, A::Atom)
	println("Levels: \n",["\t"*L.name for L in A.levels]...)
	println("Transitions: \n")
	for i in 1:length(A.levels)
		println("\t"*A.levels[i].name*" <--")
		for j in 1:length(A.levels)
			if A.linewidths[j,i]>0
				println("\t\t"*A.levels[j].name*": "*string(A.linewidths[j,i]))
			end
		end
	end
end

# We'll build the relevant parts of Sr88 now:
Sr88_1S0 = Level("1S0",0,0,0,energy=0.0)
Sr88_1P1 = Level("1P1",0,1,1,energy=21698.452e-7)
Sr88_3P0 = Level("3P0",1,1,0,energy=14317.507e-7)
Sr88_3P1 = Level("3P1",1,1,1,energy=14504.334e-7)
Sr88_3P2 = Level("3P2",1,1,2,energy=14898.545e-7)
Sr88_3S1 = Level("3S1",1,0,1,energy=29038.773e-7)

Sr88_linewidths = [0 0 0 0 0 0 ;				# Decay channels from the ground state (none)
				2.01e8 0 0 0 0 0;			# ... from 1P1
				1/157 0 0 0 0 0 ;			# ... from 3P0
				4.6911e4 0 0 0 0 0;			# ... from 3P1
				0 0 0 0 0 0;				# ... from 3P2
				0 0  8.98e6 2.7e7 4.2e7 0;]	# ... from 3S1

Sr88 = Atom([Sr88_1S0,Sr88_1P1,Sr88_3P0,Sr88_3P1,Sr88_3P0,Sr88_3S1],Sr88_linewidths)

struct Laser
	# For numerical stability, you should specify the k-vector only nominally, and specify the detuning separately.
	I::Number							# Intensity
	polarization::Array{Float64,1}			# Polarization vector (length gets scaled to 1 below)
	k::Array{Float64,1}					# k vector
	n::Number							# Index of refraction of ambient medium
	detuning::Number						# Detuning
	f::Number							# Frequency
	lambda::Number						# Wavelength
	Laser(I::Number,p::Array{Float64,1},k::Array{Float64,1},n::Number,detuning::Number) = 
		new(I,p/norm(p),k,n,detuning,norm(k)*c/(2pi*n),2pi/norm(k))
end

function Laser(intensity::Number,polarization::Array{<:Number,1};
					khat::Array{<:Number,1}=[1.0,0,0], f::Number, detuning::Number=0.0, n::Number=1.0)
	if f==nothing
		throw(ArgumentError("f or k must be specified."))
	else
		return Laser(intensity,convert(Array{Float64,1},polarization),khat*f*n*2pi,1.0,detuning)
	end
end

function Laser(A::Atom, ground::Level, excited::Level,
				saturation::Number, polarization::Array{<:Number,1}, detuning::Number, khat=[1.0,0,0], n=1.0)
	# Generates a laser using parameters of an atomic transition between "ground" and "excited" Levels.  
	# saturation specifies the laser intensity in units of the saturation intensity.
	# detuning specifies the laser frequency relative to the transition frequency, in units of linewidths.
	i = findfirst(x->x==ground,A.levels)
	j = findfirst(x->x==excited,A.levels)
	return Laser(saturation*Isat(A,ground,excited), polarization, khat=khat, f=c/A.lambda[j,i], n=n, detuning=detuning)
end

function Isat(A::Atom,ground::Level,excited::Level)
	# Returns the saturation intensity of the transition between "ground" and "excited" Levels.
	i = findfirst(x->x==ground,A.levels)
	j = findfirst(x->x==excited,A.levels)
	return (pi/3) * h*c*A.linewidths[j,i]/A.lambda[j,i]^3
end

#----------------- Bloch dynamics -------------------







end