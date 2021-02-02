import math
def generate_arc_vectors(s_positions,mu_positions,e_positions):
# takes the start, mid-arc, and end atoms to be averaged to produce parametric arc that intersects all points of a structure
# returns M, S and Mu
	i=0
	s_vector=[0,0,0] # note this is the start w.r.t the origin, not the S vector. i.e. is the <Al> positions
	for atom in s_positions:
		for j in range(3):
			s_vector[j]+=float(atom[j+1])
		i+=1
	for j in range(3):
		s_vector[j]/=i
	i=0
	mu_vector=[0,0,0]
	for atom in mu_positions:
		for j in range(3):
			mu_vector[j]+=float(atom[j+1])
		i+=1
	for j in range(3):
		mu_vector[j]/=i
	i=0
	e_vector=[0,0,0]
	for atom in e_positions:
		for j in range(3):
			e_vector[j]+=float(atom[j+1])
		i+=1
	for j in range(3):
		e_vector[j]/=i
	i=0
	M_vector,S_vector,Mu_vector=[0,0,0],[0,0,0],[0,0,0] # these are the vectors to: (i) The the midpoint of start and end of arc, (ii) the start of the arc from the midpoint, (iii) the top of the arc from the midpoint
	for j in range(3):
		M_vector[j]=(e_vector[j]+s_vector[j])/2
		S_vector[j]=s_vector[j]-M_vector[j]
		Mu_vector[j]=mu_vector[j]-M_vector[j]
	return M_vector,S_vector,Mu_vector
def project_onto_arc(M_vector,S_vector,Mu_vector,migrating_coords,t_steps=1800):
# this takes the M, S, and Mu vectors (later are relative to M), and the migrating coords/ vector of the metal which is being projected onto the curve
# returns vector of projection, length of the projection along curve C, the t value associated with the projection, the total length of the arc, and the magnitude of the error vector
	t=0
	error_array=[]
	error_vec=[0,0,0]
	unit_tangent=[0,0,0]
	curve_t=[0,0,0]
	norm=0
	dot_prod=0
	pi=math.pi
	total_arc_length=0
	length=0
	projection_vec=[0,0,0]
	while t < pi:
		for j in range(3):
			curve_t[j]=S_vector[j]*math.cos(t)+Mu_vector[j]*math.sin(t)+M_vector[j]
			error_vec[j]=migrating_coords[j]-curve_t[j]
			unit_tangent[j]=(-S_vector[j]*math.sin(t)+Mu_vector[j]*math.cos(t))
			norm+=unit_tangent[j]**2
		norm=math.sqrt(norm)
		total_arc_length+=norm*pi/t_steps
		for j in range(3):
			unit_tangent[j]/=norm
		norm=0
		for j in range(3):
			dot_prod+=error_vec[j]*unit_tangent[j]
			norm+=error_vec[j]**2
		norm=math.sqrt(norm)
		length=total_arc_length
		error_array.append([dot_prod,error_vec,norm,length,t])
	#	error_array.append(curve_t)
		t+=pi/t_steps
	error_norm=norm
	for residual in error_array:
		if residual[2] < error_norm:
			error_vec=residual[1]
			error_norm=residual[2]
			projection_length=residual[3]
			projection_t=residual[4]
	for j in range(3):
		projection_vec[j]=S_vector[j]*math.cos(projection_t)+Mu_vector[j]*math.sin(projection_t)+M_vector[j]
	return projection_vec, projection_length, projection_t, total_arc_length, error_vec, error_norm