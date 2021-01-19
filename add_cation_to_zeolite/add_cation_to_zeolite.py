import copy
import math
import plane_methods as pm
import itertools
import os
import crystal_methods as cm
import sys
#files={ # must be stripped pcoord file
#	'SSZ_13_T012_r3.28.pcoord' : "035",
#  'SSZ_13_T017_r3.24.pcoord' : "017",
#	'SSZ_13_T032_r3.07.pcoord' : "032",
#	'SSZ_13_T035_r3.07.pcoord' : "035"
#}
#for line in sys.stdin:
#  total_energies.append(float(line))
#	files.append(line)
#for fileName in files
directory=sys.argv[1] # note the ENTIRE path has to explicitly given
fileName=directory+"/"+sys.argv[2]
print("Performing analysis on file : {}".format(fileName))
#lattice_paramaters=[[14.919000,0.0,0.0],[0.0,14.919000,0.0],[0.0,0.0,14.919000]] # |A|, |B|, |C|
lattice_paramaters=[[13.6750001907,0.0000000000,0.0000000000],
[-6.8375000954,11.8428975619,0.0000000000],
[0.0000000000,0.0000000000,14.7670001984]]
periodic_conditions=[1,1,1] ## which periodic images to include, note this is a "yes" "no" set up and I recommend supercelling first if you need multiple unit cell interactions
# note coordiantes have to wrapped to the cell first which this script does automatically. This crystal_methods will automatically wrap atoms inside the cell to ensure that the correct images are evalauted and will print them wrapper

#               #
#               #
#								#
#################
#       				#
##             ##
### main file ###
##						 ##
#								#
#################
#               #
#               #
#               #


#system=cm.generate_crystal(lattice_paramaters,periodic_conditions,fileName)
system=cm.generate_crystal(lattice_paramaters,fileName)

base_T_coords=[]
base_O_coords=[]
image_T_coords=[]
image_O_coords=[]

relT_coords=[]
# relT_coords.append([i,T_arr,O_neighbours])
#                   ^ ^     ^
#                   0 1     3 (is a 4-element array)
relO_coords=[]
# relO_coords.append([j,O_arr,i,rji,rel_vectori,k,rik,rel_vectork])
#                     ^ ^     ^ ^   ^           ^ ^   ^    
#                     0 1     2 3   4           5 6   7
relO_image=[]
# relO_image.append([j,O_arr,i,rji,rel_vectori,k,rik,rel_vectork])
#                     ^ ^     ^ ^   ^           ^ ^   ^    
#                     0 1     2 3   4           5 6   7

base_structure=system.pop(0)
for atom in base_structure[1]:
	if atom[0]=='O' or atom[0]=='Ne':
		base_O_coords.append(atom)
	elif atom[0]=='Si' or atom[0]=='Al' or atom[0]=='Xe':
		base_T_coords.append(atom)
	else:
		raise Exception("Incorrect atoms in the base system")

for image in system:
  for atom in image[1]: # the base coordinates
    if atom[0]=='O' or atom[0]=='Ne':
      image_O_coords.append(atom)
    elif atom[0]=='Si' or atom[0]=='Al' or atom[0]=='Xe':
      image_T_coords.append(atom)
    else:
      raise Exception("Incorrect atoms in the images")
for i, T_atom in enumerate(base_T_coords):
	T_arr=[float(T_atom[n]) for n in range(1,4)]
	O_neighbours=[]
	base_O_coords_eval=copy.deepcopy(base_O_coords)
	base_T_coords_eval=copy.deepcopy(base_T_coords)
	
	for j, O_atom in enumerate(base_O_coords_eval): # relativise
		O_arr=[float(O_atom[m]) for m in range(1,4)]
		rji=pm.separation(T_arr,O_arr)
		replication=False
		if rji < 2:
			O_neighbours.append(j)
			replication=False
			rel_vectori=pm.obtain_relative_T_vectors(O_arr,T_arr)
			for prev_O_atom in relO_coords:
				if rji==prev_O_atom[4] or rel_vectori==prev_O_atom[7]:
					replication=True
			if replication == True:
				continue
			k=0; rjk=0; rel_vectork=[0,0,0]
			relO_coords.append([j,O_arr,i,rji,rel_vectori,k,rjk,rel_vectork])
#													^	^			^	^		^					 	^ ^		^		 
#													0	1			2	3		4						5 6		7
			for k, T2_atom in enumerate(base_T_coords_eval):
				if T2_atom == T_atom:
					continue
				T2_arr=[float(T2_atom[n]) for n in range(1,4)]
				rjk=pm.separation(O_arr,T2_arr)
				if rjk <2:
					rel_vectork=pm.obtain_relative_T_vectors(O_arr,T2_arr)
					relO_coords[-1][5]=k
					relO_coords[-1][6]=rjk
					relO_coords[-1][7]=rel_vectork
			for k, T2_atom in enumerate(image_T_coords):
				T2_arr=[float(T2_atom[n]) for n in range(1,4)]
				rjk=pm.separation(O_arr,T2_arr)
				if rjk <2:
					rel_vectork=pm.obtain_relative_T_vectors(O_arr,T2_arr)
					relO_coords[-1][5]="i{}".format(k)
					relO_coords[-1][6]=rjk
					relO_coords[-1][7]=rel_vectork
		else:
			for q, Ti_atom in enumerate(image_T_coords):
				Ti_arr=[float(Ti_atom[n]) for n in range(1,4)]
				rjq=pm.separation(O_arr,Ti_arr)
				replication=False
				if rjq < 2:
					rel_vectorq=pm.obtain_relative_T_vectors(O_arr,Ti_arr)
					for prev_O_atom in relO_coords:
						if rjq == prev_O_atom[6]:
							replication=True
					if replication == True:
						continue
					x=0; rjx=0; rel_vectorx=[0,0,0]
					rel_vectorq=pm.obtain_relative_T_vectors(O_arr,Ti_arr)
					relO_coords.append([j,O_arr,"i{}".format(q),rjq,rel_vectorq,x,rjx,rel_vectorx])
#	 													^ ^		 ^ ^	 ^					 ^ ^	 ^		
#	 													0 1		 2 3	 4					 5 6	 7
					for x, T2_atom in enumerate(base_T_coords_eval):
						if x == q:
							continue
						T2_arr=[float(T2_atom[n]) for n in range(1,4)]
						rjx=pm.separation(O_arr,T2_arr)
						if rjx <2:
							rel_vectorx=pm.obtain_relative_T_vectors(O_arr,T2_arr)
							relO_coords[-1][5]=x
							relO_coords[-1][6]=rjx
							relO_coords[-1][7]=rel_vectorx

	for j, O_atom in enumerate(image_O_coords): # relativise
		O_arr=[float(O_atom[m]) for m in range(1,4)]
		rji=pm.separation(T_arr,O_arr)
		if rji < 2:
			O_neighbours.append("i{}".format(j))
	relT_coords.append([i,T_arr,O_neighbours])
for j, Oim_atom in enumerate(image_O_coords):
	Oim_arr=[float(Oim_atom[i]) for i in range (1,4)]
	base_T_coords_eval=copy.deepcopy(base_T_coords)
	image_T_coords_eval=copy.deepcopy(image_T_coords)
	flag=0
	for i, Ti_atom in enumerate(base_T_coords_eval):
		T_arr=[float(Ti_atom[p]) for p in range(1,4)]
		rji=pm.separation(T_arr,Oim_arr)
		replication=False
		if rji < 2:
			rel_vectori=pm.obtain_relative_T_vectors(Oim_arr,T_arr)
			for prev_O_atom in relO_image: 
				if  rel_vectori == prev_O_atom[4] or rel_vectori == prev_O_atom[7]:
					replication=True
			if replication == True:
				continue
			k=0; rjk=0; rel_vectork=[0,0,0]
			relO_image.append([j,Oim_arr,i,rji,rel_vectori,k,rjk,rel_vectork])
#			 	 								 ^ ^		 	 ^ ^	 ^					 ^ ^	 ^		
#												 0 1			 2 3	 4					 5 6	 7
			for k, T2_atom in enumerate(base_T_coords_eval):
				if T2_atom == Ti_atom:
					continue
				T2_arr=[float(T2_atom[n]) for n in range(1,4)]
				rjk=pm.separation(Oim_arr,T2_arr)
				if rjk <2:
					rel_vectork=pm.obtain_relative_T_vectors(Oim_arr,T2_arr)
					relO_image[-1][5]=k
					relO_image[-1][6]=rjk
					relO_image[-1][7]=rel_vectork
					flag=1
			for k, T2_atom in enumerate(image_T_coords_eval):
				T2_arr=[float(T2_atom[n]) for n in range(1,4)]
				rjk=pm.separation(Oim_arr,T2_arr)
				if rjk <2:
					rel_vectork=pm.obtain_relative_T_vectors(Oim_arr,T2_arr)
					relO_image[-1][5]="i{}".format(k)
					relO_image[-1][6]=rjk
					relO_image[-1][7]=rel_vectork
					flag=1
		if flag ==1:
			break
	for q, Tq_atom in enumerate(image_T_coords_eval):
		Tq_arr=[float(Tq_atom[n]) for n in range(1,4)]
		rjq=pm.separation(Oim_arr,Tq_arr)
		replication=False
		if rjq < 2:
			rel_vectorq=pm.obtain_relative_T_vectors(Oim_arr,Tq_arr)
			for prev_O_atom in relO_image:
				if rel_vectorq == prev_O_atom[4] or rel_vectorq == prev_O_atom[7]:
					replication=True
			if replication == True:
				continue
			x=0; rjx=0; rel_vectorx=[0,0,0]
			relO_image.append([j,Oim_arr,"i{}".format(q),rjq,rel_vectorq,x,rjx,rel_vectorx])
#												 ^ ^				^							 ^	 ^					 ^ ^	 ^		
#												 0 1				2							 3	 4					 5 6	 7
			for x, T2_atom in enumerate(image_T_coords_eval):
				if T2_atom == Tq_atom:
					continue
				T2_arr=[float(T2_atom[n]) for n in range(1,4)]
				rjx=pm.separation(Oim_arr,T2_arr)
				if rjx <2:
					rel_vectorx=pm.obtain_relative_T_vectors(Oim_arr,T2_arr)
					relO_image[-1][5]=x
					relO_image[-1][6]=rjx
					relO_image[-1][7]=rel_vectorx
H_coords=[]
for relT_row in relT_coords:
	# relT_row[i, T_arr, [O_neighbours values]
	t=relT_row[0]
	entry=[]
	if base_T_coords[t][0] != 'Al':
		continue
	entry.append(t)
	O_neighbours=relT_row[2]
	for jVal in O_neighbours:
		if type(jVal) == int:	
			for relO_row in relO_coords:
				if relO_row[0] == jVal:
					Ti_vec=relO_row[4]
					Tk_vec=relO_row[7]
					O_pos=relO_row[1]
					H_coord=pm.obtain_unit_vec_trigPlanar2(Ti_vec,Tk_vec,[0,0,0])
					for j in range(3):
						H_coord[j]+=O_pos[j]
					break
		else:
			jVal=int(jVal.replace('i',''))
			for relO_row in relO_image:
				if relO_row[0] == jVal:
					Ti_vec=relO_row[4]
					Tk_vec=relO_row[7]
					O_pos=relO_row[1]
					H_coord=pm.obtain_unit_vec_trigPlanar2(Ti_vec,Tk_vec,[0,0,0])
					for j in range(3):
						H_coord[j]+=O_pos[j]
					break
		entry.append([H_coord,jVal]) # obtain relVec1 from relO_coords where the t is the i or k value
	H_coords.append(entry)
# for row in H_coords:
# 	for j in range(1,5):
# 		H_arrangement=row[j]
# 		H_coords=H_arrangement[0]
# 		print("H {} {} {}".format(H_coords[0],H_coords[1],H_coords[2]))
# for row in base_T_coords:
# 	print("Si {} {} {}".format(row[1],row[2],row[3]))
# for row in base_O_coords:
# 	print("O {} {} {}".format(row[1],row[2],row[3]))

N=len(H_coords)
permutations=[]
for seq in itertools.product("1234", repeat=N):
	string=""
	for x in seq:
		string+=x
	entry=[string]
	prev_O_site=-1
	for mu in range(N):
		seq_value=int(seq[mu])
		Bronsted_row=H_coords[mu]
		T_value=Bronsted_row[0]
		current_O_site=Bronsted_row[seq_value][1]
		if current_O_site== prev_O_site:
			for relO_row in relO_coords:
				if (relO_row[0]==current_O_site) and (relO_row[2]==T_value or relO_row[5]==T_value):
					O_coord=relO_row[1]
					plane=pm.find_plane(relO_row[4],relO_row[7],[0,0,0])
					H_position=pm.normalise_vec([plane[0],plane[1],plane[2]])
					for x in range(3):
						H_position[x]+=O_coord[x]
					
		else:
			H_position=H_coords[mu][seq_value][0]
		prev_O_site=current_O_site
		entry.append(H_position)
	permutations.append(entry)

for permutation in permutations:
	filePermut=sys.argv[1]+"/H"+permutation[0]+"_"+sys.argv[2]
#	os.chdir("T{}".format(files[fileName]))
	command="cp {} {}".format(fileName,filePermut)
	os.system(command)
	fh=open(filePermut,"a+")
	for q in range(N):
		H_atom=permutation[q+1]
		Bronsted_site="H {} {} {}\n".format(H_atom[0],H_atom[1],H_atom[2])
		fh.write(Bronsted_site)
	fh.close()
	os.chdir("..")
print("Finished analysis on : {}".format(fileName))
