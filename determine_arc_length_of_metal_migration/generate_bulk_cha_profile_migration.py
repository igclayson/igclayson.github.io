import generate_bulk_cha_profile_migration_mod as profile_mig
import math
import sys
s_atom_list={
	56:[0,0,0],
	58:[0,0,0],
	68:[0,0,0],
	70:[0,0,0],
	80:[0,0,0],
	82:[0,0,0]
}
mu_atom_list={
	13:[1,0,1],
	15:[1,0,1],
	24:[0,0,1],
	38:[1,0,1],
	39:[1,0,1],
	47:[1,0,0],
	64:[1,0,0],
	65:[1,0,0]
	
}
e_atom_list={
	52:[1,0,0],
	54:[1,0,0],
	64:[1,0,0],
	66:[1,0,0],
	76:[1,0,0],
	78:[1,0,0]
}
lattice_param=[[13.816308257179562,0.0,0.0],
[-6.8619717247947714,11.885283667445869,0.0],
[0.0,0.0,14.494423906405032]]
fileName="CuH2O6_fixed_He_moves_anim.xyz"
migration_atom='Cu'
zero_arc_length=0

fh=open(fileName,"r")
snapshot_array=[]
atoms_array=[]
current_atom_number=1
atom_number=0
for i, line in enumerate(fh):
	entry=line.split()
	if i==0:
		atom_number=int(line)
		continue
	if len(entry) != 4:
		continue
	if current_atom_number < atom_number:
		atoms_array.append(entry)
		current_atom_number+=1
	else:
		atoms_array.append(entry)
		snapshot_array.append(atoms_array)
		current_atom_number=1
		atoms_array=[]
		continue
fh.close()
projection_array=[]
#num_of_relAtoms=len(relative_atom_list)
for ff, frame in enumerate(snapshot_array):
	migrating_coords=[0,0,0]
	plane_av_coords=[0,0,0]
	s_positions, mu_positions, e_positions=[],[],[]
	for k, atom in enumerate(frame):
		# need to check for s, mu, and e positons
		for atom_num in s_atom_list:
			if (k+1)==atom_num:
# need to correct for crystal image
				if s_atom_list[atom_num] != [0,0,0]:
					for n in range(3):
						for m in range(3):
							atom[m+1]=float(atom[m+1])+float(s_atom_list[atom_num][n])*float(lattice_param[n][m])
				s_positions.append(atom)
		for atom_num in mu_atom_list:
			if (k+1)==atom_num:
# need to correct for crystal image
				if mu_atom_list[atom_num] != [0,0,0]:
					for n in range(3):
						for m in range(3):
							atom[m+1]=float(atom[m+1])+float(mu_atom_list[atom_num][n])*float(lattice_param[n][m])
				mu_positions.append(atom)
		for atom_num in e_atom_list:
			if (k+1)==atom_num:
# need to correct for crystal image
				if e_atom_list[atom_num] != [0,0,0]:
					for n in range(3):
						for m in range(3):
							atom[m+1]=float(atom[m+1])+float(e_atom_list[atom_num][n])*float(lattice_param[n][m])
				e_positions.append(atom)
		if atom[0]==migration_atom:
			for j in range(3):
				migrating_coords[j]=float(atom[j+1])
#				print(migrating_coords[j],end=" ")
#			print("")
	M_vector,S_vector,Mu_vector=profile_mig.generate_arc_vectors(s_positions,mu_positions,e_positions) # relative vectors to generate the arc which now the arc must be generated via the parametric curve
	#print(M_vector,S_vector,Mu_vector)
	#print(profile_mig.project_onto_arc(M_vector,S_vector,Mu_vector,migrating_coords))
	#projection_vec, projection_length, projection_t, total_arc_length, error_vec, error_norm
	#projection_vec, projection_length, projection_t, total_arc_length, error_vec, error_norm=profile_mig.project_onto_arc(M_vector,S_vector,Mu_vector,migrating_coords)
	if ff==0:
		projection_vec, projection_length, projection_t, total_arc_length, error_vec, error_norm=profile_mig.project_onto_arc(M_vector,S_vector,Mu_vector,migrating_coords)
		zero_arc_length=total_arc_length
	else:
		projection_vec, projection_length, projection_t, total_arc_length, error_vec, error_norm=profile_mig.project_onto_arc(M_vector,S_vector,Mu_vector,migrating_coords)
	#print("Cl", end=" ")
	print(projection_length)
#	for j in range(3):
#		print(migrating_coords[j],end=" ")
#	print("")
#print(total_arc_length)
#	r=0
#	for j in range(3):
#		plane_av_coords[j]=plane_av_coords[j]/num_of_relAtoms
#		r+=(migrating_coords[j]-plane_av_coords[j])**2
#	average_separation_array.append(math.sqrt(r))
