import copy
import itertools
def generate_crystal(lattice_paramaters,fileName): # this is an old one method that doesn't generate the cross-terms very well
	system=[]
	list_images=[]
	for seq in itertools.product("012", repeat=3):
		if seq==('1','1','1'):
			continue
		list_images.append([int(i)-1 for i in seq])
	# generating the 000 image
	O_coords=[]; T_coords=[]; base_coords=[]
	for i, line in enumerate(open(fileName, 'r')):
		entry=line.split()
		for i in range(1,4):
			entry[i]=float(entry[i])
		base_coords.append(entry)
# for atom in base_coords:
#	 for i in range(3):
#		 if atom[i+1] < 0:
#			 print(atom)
#			 for j in range(3):
#				 atom[j+1]+=lattice_paramaters[i][j]
#		 if atom[i+1] > lattice_paramaters[i][i]:
#			 for j in range(3):
#				 atom[j+1]-=lattice_paramaters[i][j]
	O_coords= copy.deepcopy(base_coords)
	for atom in O_coords:
		#if 'Si' in atom or 'Al' in atom:
		if 'Si' in atom or 'Al' in atom or 'Xe' in atom: # note that He acts as a T atom ghost
			T_coords.append(atom)
	for atom in T_coords:
		O_coords.remove(atom)
	# generating system array
	system.append([[0,0,0],base_coords])
	for cell_point_array in list_images:
		image_coords=[]
		image_coords=copy.deepcopy(base_coords)
		for atom in image_coords:
			TrsVec=[0.0,0.0,0.0] # this is the additive
			base_atom=copy.deepcopy(atom)
			for i in range(3):
				for j in range(3):
					TrsVec[i]+=lattice_paramaters[j][i]*cell_point_array[j]
				atom[i+1]+=TrsVec[i]
		system.append([cell_point_array,image_coords])
	return system # , T_coords, O_coords # can return the T and O coords in the base if desired
def old_generate_crystal(lattice_paramaters,periodic_conditions,fileName): # this is an old one method that doesn't generate the cross-terms very well
	seq_length=0
	system=[]
	for j in periodic_conditions:
		 if j==1:
			 seq_length+=1
	list_images=[]
	for seq in itertools.product("01", repeat=seq_length):
		check=0
		image_condition_p=[0,0,0]
		image_condition_n=[0,0,0]
		for j in range(3):
			if (check) == seq_length:
				break
			if periodic_conditions[j] == 1:
				if int(seq[check]) == 0:
					image_condition_p[j]=0
					image_condition_n[j]=0
				else:
					image_condition_p[j]=1
					image_condition_n[j]=-1
				check+=1
			else:
				image_condition_p[j]=0
				image_condition_n[j]=0
		if image_condition_p == image_condition_n:
			list_images.append(image_condition_p)
		else:
			list_images.append(image_condition_p) # image list is the array of all the images considered so 2nd order array
			list_images.append(image_condition_n)
	# generating the 000 image
	O_coords=[]; T_coords=[]; base_coords=[]
	for i, line in enumerate(open(fileName, 'r')):
		entry=line.split()
		for i in range(1,4):
			entry[i]=float(entry[i])
		base_coords.append(entry)
#	for atom in base_coords:
#		for i in range(3):
#			if atom[i+1] < 0:
#				print(atom)
#				for j in range(3):
#					atom[j+1]+=lattice_paramaters[i][j]
#			if atom[i+1] > lattice_paramaters[i][i]:
#				for j in range(3):
#					atom[j+1]-=lattice_paramaters[i][j]
	O_coords= copy.deepcopy(base_coords)
	for atom in O_coords:
		#if 'Si' in atom or 'Al' in atom:
		if 'Si' in atom or 'Al' in atom or 'Xe' in atom: # note that He acts as a T atom ghost
			T_coords.append(atom)
	for atom in T_coords:
		O_coords.remove(atom)
	# generating system array
	for cell_point_array in list_images:
		image_coords=[]
		image_coords=copy.deepcopy(base_coords)
		for atom in image_coords:
			TrsVec=[0.0,0.0,0.0] # this is the additive
			base_atom=copy.deepcopy(atom)
			for i in range(3):
				for j in range(3):
					TrsVec[i]+=lattice_paramaters[j][i]*cell_point_array[j]
				atom[i+1]+=TrsVec[i]
		system.append([cell_point_array,image_coords])
	return system # , T_coords, O_coords # can return the T and O coords in the base if desired
