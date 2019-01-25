# multilevel.jl
# Evaluates the dynamics of a multilevel atom in the presence of lasers and magnetic fields.


struct Level
	name::String
	S::Number				# Spin angular momentum
	L::Number				# Orbital...
	J::Number				# Total electronic...
	I::Number				# Nuclear...
	F::Number				# Nuclear + electronic...
	energy::Number			# Energy with respect to ionization
	E::Number				# Alias for energy
	Level(s,l,j,i,f,e) = new(s,l,j,i,f,e,e)
end

struct Atom
	levels::Array{Level,1}				# Levels
	# Linewidths between each pair of levels:
		# linewidths[i,j] is the linewidth from i to j.
		# So if i is lower energy than j, it is 0.
	linewidths::Array{Number,2}
	lambda::Array{Number,2}				# Vacuum wavelength connecting any pair of states
	function Atom(
end

struct Laser
	k::Array{Number,1}					# k vector
	polarization::Array{Number,1}		# Polarization vector
	n::Number							# Index of refraction of ambient medium
	f::Number							# Frequency
	lambda::Number						# Wavelength
end

#struct Field							# The field is just a vector, so it probably doesn't need its own struct.
#end
