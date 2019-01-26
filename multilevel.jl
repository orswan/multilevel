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
		# linewidths[i,j] is the linewidth from i to j.
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
	k::Array{Number,1}					# k vector
	polarization::Array{Number,1}			# Polarization vector
	n::Number							# Index of refraction of ambient medium
	f::Number							# Frequency
	lambda::Number						# Wavelength
	Laser(k,p,n,f,l) = new(k,p,n,norm(k)*c/(2pi*n),2pi/norm(k))
end

#struct Field							# The field is just a vector, so it probably doesn't need its own struct.
#end




end