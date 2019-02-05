# multilevel.jl
# Evaluates the dynamics of a multilevel atom in the presence of lasers and magnetic fields.
# We use units with hbar=c=1
module multilevel
using LinearAlgebra#,DifferentialEquations
using WignerSymbols
using Unitful

#------------------------ Units --------------------------

"""
The units I would like to use are such that h=c=1, 
a unit magnetic field generates a unit energy for a Bohr magneton.
"""

const h = 1
const c = 1
const hbar = h/2pi
const muB = 1

function unitConvert(s::String,baseUnits,targetUnits)
	# Converts a dimension given as s from baseUnits to targetUnits.
		# Symbols in the string should be delimited by spaces.
		# HACKY....
	SI_h = 6.626176e-34
	SI_hbar = SI_h/2pi
	SI_c = 299792458 
	pref = Dict(""=>1, "c"=>.01, "m"=>.001, "u"=>1e-6, "n"=>1e-9, "k"=>1000, "M"=>1e6, "G"=>1e9)
	if (baseUnits==:SI) & (targetUnits==:natural)
		ex = 1
		out = 1
		for i in split(s,' ')
			if i[end]=='s'
				out = out*pref[i[1:end-1]]^ex
			elseif i[end]=='m'
				out = out*(pref[i[1:end-1]]/SI_c)^ex
			elseif i[end]=='g'
				throw(ArgumentError,"mass not implemented in unit conversion yet")
			elseif i[end]=='/'
				ex = -1
			end
		end
		return out
	end
end

#--------------------- Structs -----------------------
struct State
	name::String
	S::Number				# Spin angular momentum
	L::Number				# Orbital...
	J::Number				# Total electronic...
	I::Number				# Nuclear...
	F::Number				# Nuclear + electronic...
	energy::Number			# Energy (can be w.r.t. any fixed zero point)
	E::Number				# Alias for energy
	State(name,s,l,j,i,f,e) = new(name,s,l,j,i,f,e,e)
end

State(name,s,l,j,i,f) = State(name,s,l,j,i,f,nothing)			# Energy need not be specified
State(name,s,l,j;energy) = State(name,s,l,j,0,j,energy)		# Nuclear moment need not be specified

struct Atom
	states::Array{State,1}				# States
	# Linewidths between each pair of states:
		# linewidths[i,j] is the ANGULAR linewidth from i to j.
		# So if i is lower energy than j, it is 0.
	linewidths::Array{Float64,2}
	lambda::Array{Float64,2}			# Vacuum wavelength connecting any pair of states
end

function Atom(states::Array{State,1},linewidths::Array{<:Number,2})
	# Generates an Atom object using the energy of the states.  Fails if not all energies are specified. 
	energies = [L.energy for L in states]
	if any(x->x==nothing,energies)			# Check that all energies are defined
		throw(ArgumentError("Not all State energies are defined."))
	else
		return Atom(states,
					convert(Array{Float64,2},linewidths),
					convert(Array{Float64,2},[h*c/abs(L1.energy-L2.energy) for L1 in states, L2 in states])
					)
	end
end

function Atom(A::Atom,states::Array{State,1})
	# Generates a new Atom from a subset of the states of an existing atom.
	idxs = [findfirst(x->x==L,A.states) for L in states]
	return Atom(states,A.linewidths[idxs,idxs],A.lambda[idxs,idxs])
end

function Base.show(io::IO, A::Atom)
	println("States: \n",["\t"*L.name for L in A.states]...)
	println("Transitions: \n")
	for i in 1:length(A.states)
		println("\t"*A.states[i].name*" <--")
		for j in 1:length(A.states)
			if A.linewidths[j,i]>0
				println("\t\t"*A.states[j].name*": "*string(A.linewidths[j,i]))
			end
		end
	end
end

mutable struct Laser
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

function Laser(A::Atom, ground::State, excited::State,
				saturation::Number, polarization::Array{<:Number,1}, detuning::Number, khat=[1.0,0,0], n=1.0)
	# Generates a laser using parameters of an atomic transition between "ground" and "excited" States.  
	# saturation specifies the laser intensity in units of the saturation intensity.
	# detuning specifies the laser frequency relative to the transition frequency, in units of linewidths.
	i = findfirst(x->x==ground,A.states)
	j = findfirst(x->x==excited,A.states)
	return Laser(saturation*Isat(A,ground,excited), polarization, khat=khat, f=c/A.lambda[j,i], n=n, detuning=detuning)
end

function Laser(A::Atom,ground::Int,excited::Int,
				saturation::Number,polarization::Array{<:Number,1}, detuning::Number,khat=[1.0,0,0],n=1.0)
	return Laser(saturation*Isat(A,A.states[ground],A.states[excited]), polarization, khat=khat, f=c/A.lambda[excited,ground], n=n, detuning=detuning)
end

function Base.show(io::IO, L::Laser)
	println("lambda: \t",L.lambda)
	println("Detuning: \t",L.detuning)
	println("Intensity: \t",L.I)
	println("Polarization: \t",L.polarization)
end

#----------------- Strontium -----------------------

Sr88_1S0 = State("1S0",0,0,0,energy=0.0)
Sr88_1P1 = State("1P1",0,1,1,energy=21698.452e-7)
Sr88_3P0 = State("3P0",1,1,0,energy=14317.507e-7)
Sr88_3P1 = State("3P1",1,1,1,energy=14504.334e-7)
Sr88_3P2 = State("3P2",1,1,2,energy=14898.545e-7)
Sr88_3S1 = State("3S1",1,0,1,energy=29038.773e-7)

Sr88_linewidths = [0 0 0 0 0 0 ;				# Decay channels from the ground state (none)
				2.01e8 0 0 0 0 0;				# ... from 1P1
				1/157 0 0 0 0 0 ;				# ... from 3P0
				4.6911e4 0 0 0 0 0;				# ... from 3P1
				0 0 0 0 0 0;					# ... from 3P2
				0 0  8.98e6 2.7e7 4.2e7 0;]*	# ... from 3S1
				unitConvert("s / nm",:SI,:natural)

Sr88 = Atom([Sr88_1S0,Sr88_1P1,Sr88_3P0,Sr88_3P1,Sr88_3P0,Sr88_3S1],Sr88_linewidths)
Sr88_689 = Atom(Sr88,[Sr88_1S0,Sr88_3P1])

#----------------- Utility -------------------------

function Isat(A::Atom,ground::State,excited::State)
	# Returns the saturation intensity of the transition between "ground" and "excited" States.
	i = findfirst(x->x==ground,A.states)
	j = findfirst(x->x==excited,A.states)
	if A.linewidths[j,i]==0
		println("WARNING: Linewidth zero encountered.  You may have mixed up your ground and excited states.")
	end
	return (pi/3) * h*c*A.linewidths[j,i]/A.lambda[j,i]^3
end

function Lande(S::State)
	# Lande g-factor for a state S.
	if S.J==0
		return 0
	else
		return 3/2 + (S.S*(S.S+1) - S.L*(S.L+1))/(2*S.J*(S.J+1))
	end
end

function magMoment(ground::State,gm::Number,excited::State,em::Number)
	# Computes magnetic moment between two states.
	return (em*Lande(excited) - gm*Lande(ground)) * muB
end

function sphericalBasis(axis::Array{<:Number,1})
	# Returns a set of spherical basis vectors along given axis.  
		# The particular basis vectors are generated by reflecting the canonical z-axis onto the given axis,
		# then swapping the x and y directions, and finally forming (x\pm iy)/sqrt(2) for the sigma_{\pm} vectors.
		# Return format is a matrix with rows [sigma-,z,sigma+]
	if size(axis)!=(3,)
		throw(ArgumentError,"axis must have shape (3,)")
	end
	if axis/norm(axis)==[0,0,1]
		return [1/sqrt(2) 0 1/sqrt(2);im/sqrt(2) 0 -im/sqrt(2); 0 1 0]
	else
		v = axis/norm(axis)-[0,0,1]			# This will be the normal to the plane of reflection
		v = v/norm(v)							# Normalize
		xyz = LinearAlgebra.I - 2v*v'			# xyz[:,i] is the reflected e_i unit vector
		# Note that xyz[:,1] is a reflected e_1, so to get a right handed coordinate system, we should use 
			# xyz[:,2] as x, xyz[:,1] as y, and xyz[:,3] as z.
		return [(xyz[:,2]-im*xyz[:,1])/sqrt(2) xyz[:,3] (xyz[:,2]+im*xyz[:,1])/sqrt(2)]
	end
end

function proj(p::Array{<:Number,1},Bhat::Array{<:Number,1},q::Number)
	# Projects p onto the spherical basis element q with axis B, where q=\pm 1 is sigma\pm and q=0 is pi.
	return dot(p,conj(sphericalBasis(Bhat)[:,2+q]))
end

function RabiFreq(A::Atom,ground::State,excited::State,L::Laser,Bhat::Array{<:Number,1};gm=nothing,em=nothing)
	# Computes Rabi frequency between ground and excited states from a given laser L.
		# Computes for given m states if provided.  Otherwise returns the Wigner-Eckart matrix element.
		# Currently doesn't work for half integers.
	i = findfirst(x->x==ground,A.states)
	j = findfirst(x->x==excited,A.states)
	Omega = sqrt(A.linewidths[j,i]^2 * L.I/(Isat(A,ground,excited)*2))		# Angular Rabi frequency
	if !(gm==nothing) & !(em==nothing)
		if (abs(gm-em)<=1)
			# TODO: Figure out which way the Clebsch-Gordan coefficients should go. ###############################
			Omega = Omega * clebschgordan(ground.J,gm,1,em-gm,excited.J,em) * proj(L.polarization,Bhat,em-gm)
		else			# If the difference in m states is not in {-1,0,1}, then the dipole matrix element is zero.
			Omega = 0
		end
	end
	return Omega
end

#----------------- Bloch dynamics -------------------

function RabiH(A::Atom,lasers::Array{Laser,1},B::Array{<:Number,1}=[0.0,0,0];Bhat=nothing)
	# Returns the matrix for the Hamiltonian given an Atom A, a list of lasers, and a magnetic field B.
		# Assumes that near detuned lasers are all that is relevant, so this captures a generalized Rabi system.
		# Currently only allows for one laser for any given transition. 
		# Currently does not account for nuclear spin.
	Jdims = [2*L.J+1 for L in A.states]
	Jfins = cumsum(Jdims)
	Jstarts = Jfins .- Jdims .+ 1
	N = Jfins[end]					# Total number of eigenstates, and dimension of the Hamiltonian matrix.
	H = zeros(N,N)					# Hamiltonian matrix
	b = norm(B)						# Magnitude of B
	if Bhat==nothing
		Bhat = (b==0 ? [0,0,1] : B/b)	# Direction of B
	end
	
	for L in lasers											# We add terms to the Hamiltonian laser by laser
		detuning = L.detuning									# Detuning of the laser from the transition
		idx = findmin(abs.(A.lambda .- L.lambda))[2]				# Index of transition closest to laser wavelength
		if A.linewidths[idx]==0								# This means idx[1] is the ground state
			gidx = idx[1]										# Ground state index in A
			eidx = idx[2]										# Excited state index in A
		else
			gidx = idx[2]
			eidx = idx[1]
		end
		ground = A.states[gidx]									# Ground state
		excited = A.states[eidx]									# Excited state
		
		for g = Jstarts[gidx]:Jfins[gidx]							# Run over all ground m levels
			gm = g .- Jstarts[gidx] .- A.states[gidx].J				# Ground state m level
			for e = Jstarts[eidx]:Jfins[eidx]						# Run over all excited m levels
				em = e .- Jstarts[eidx] .- A.states[eidx].J			# Excited state m level
				Omega = RabiFreq(A,ground,excited,L,Bhat,gm=gm,em=em)	# Rabi angular frequency
				H[g,e] = H[g,e] .+ hbar * Omega / 2				# Above diagonal
				H[e,g] = H[e,g] .+ hbar * conj(Omega) / 2			# Below diagonal
				H[e,e] = H[e,e] .- hbar * detuning .- magMoment(ground,gm,excited,em) * b	# Diagonal
			end
		end
	end
	return H
end

RabiH(A::Atom,L::Laser,B::Array{<:Number,1}=[0.0,0,0];Bhat=nothing) = RabiH(A,[L],B;Bhat=Bhat)

function RabiDynamics(ts::Array{<:Number,1},A::Atom,lasers::Union{Array{Laser,1},Laser},B::Array{<:Number,1}=[0.0,0,0];Bhat=nothing)
	# Returns population evolution for times ts, assuming initial state is the ground state.
	H = RabiH(A,lasers,B,Bhat=Bhat)			# Hamiltonian
	N = size(H)[1]						# Number of levels
	pops = zeros(length(ts),N)				# Populations
	for i in 1:length(ts)
		pops[i,:] = abs.(exp(im*ts[i]*H)[:,1]).^2
	end
	return pops
end

function RabiDynamics(ts::Union{StepRange,StepRangeLen},A::Atom,lasers::Union{Array{Laser,1},Laser},B::Array{<:Number,1}=[0.0,0,0];Bhat=nothing)
	# Returns population evolution for times ts, assuming initial state is the ground state.
	H = RabiH(A,lasers,B,Bhat=Bhat)			# Hamiltonian
	N = size(H)[1]						# Number of levels
	pops = zeros(length(ts),N)				# Populations
	pops[1,1] = 1.0						# Initial state
	expH = exp(im*H*(ts[2]-ts[1]))			# Time evolution operator for one time step
	v = zeros(N); v[1] = 1.0				# State vector
	for i in 2:length(ts)
		v = expH*v
		pops[i,:] = abs.(v).^2
	end
	return pops
end

end