import math
import numpy as np
# 1, [11.23606, 3.67248, 2.46166], 0, 1.6113680938568944, [-0.8868299999999998, 0.9912500000000004, 0.9096500000000003], 1, 1.6286227369467743, [0.9288900000000009, 0.5585600000000004, -1.21556]
def separation(arr1,arr2):
  r=0
  for j in range(3):
    r+=(arr1[j]-arr2[j])**2
  r=math.sqrt(r)
  return r
def obtain_relative_T_vectors(Oarr,Rarr): # give the relative vector w.r.t to the first argument
  retArr=[0,0,0]
  for j in range(3):
    retArr[j]=Rarr[j]-Oarr[j]
  return retArr
def dot_product(vec1, vec2):
	length=len(vec1)
	dotProd=0
	for j in range(length):
		dotProd+=vec1[j]*vec2[j]
	return dotProd
def cross_product(vec1, vec2): # returns the cross-product vector of a 3d vector
	i=(vec1[1]*vec2[2]-vec1[2]*vec2[1])
	j=-(vec1[0]*vec2[2]-vec1[2]*vec2[0])
	k=(vec1[0]*vec2[1]-vec1[1]*vec2[0])
	crossVec=[i,j,k]
	return crossVec
def modulus(vec):
	length=len(vec)
	mod=0
	for j in range(length):
		mod+=vec[j]**2
	mod=math.sqrt(mod)
	return mod
def normalise_vec(vec):
  normVec=[0,0,0]
  b=1/modulus(vec)
  for j in range(3):
    normVec[j]=vec[j]*b
  return normVec
def angle_between_vecs(vec1,vec2):
	if modulus(vec1) == 0 or modulus(vec2) == 0:
		print(vec1)	
		print(vec2)
	cTheta=(dot_product(vec1,vec2))/(modulus(vec1)*modulus(vec2))
	theta=math.acos(cTheta)
	cTheta_theta=[cTheta,theta]
	return cTheta_theta
def find_plane(relVec1, relVec2, point): # returns a 4-element list for a plane of a,b,c, and d for a.x+b.y+c.z=d
	plane=cross_product(relVec1, relVec2)
	d=0
	for j in range(3):
		d+=plane[j]*point[j]
	plane.append(d)
	return plane
def obtain_unit_vec_trigPlanar(relVec1,relVec2,point): 
######  NOTE  ######
### DOESN'T WORK ###
### DON'T USE IT ###
# this will return a unit vector forming a trigonal planar with an angle of separation minus theta away from vec1
# the relative vectors are from a point 
# thsi is predominantly used for adding Bronsted to O-(T)2 and thus the relVec will need the point added to them
# plane = a b c d
#					0 1 2 3
	plane=find_plane(relVec1,relVec2,point)
	ang_array=angle_between_vecs(relVec1,relVec2)
	f=modulus(relVec1)*math.cos(-ang_array[1]/2)
	Q=1-(plane[3]*relVec1[2])/(f*plane[2])
	n_x=(relVec1[0]/f)-(relVec1[2]*plane[0])/(f*plane[2])
	n_y=(relVec1[1]/f)-(relVec1[2]*plane[1])/(f*plane[2])
	Lambda=Q/n_y
	Omega=n_x/n_y
	beta=plane[3]/plane[2]-(plane[1]/plane[2])*(Q/n_y)
	alpha=plane[0]/plane[2]+(plane[1]/plane[2])*(n_x/n_y)
	s_alpha=1+Omega**2+alpha**2
	s_beta=-2*Lambda*Omega-2*beta*alpha
	s_gamma=Lambda**2+beta**2-1
	array=[f,Q,n_x,n_y,Lambda,Omega,beta,alpha,s_alpha,s_beta,s_gamma]
	if (s_beta**2-4*s_alpha*s_gamma) < 0:
#		print("WARNING!!! NO SOLUTION")
		return 1
	x_1=(-s_beta+math.sqrt(s_beta**2-4*s_alpha*s_gamma))/(2*s_alpha)
	x_2=(-s_beta-math.sqrt(s_beta**2-4*s_alpha*s_gamma))/(2*s_alpha)
	y_1=Lambda-Omega*x_1
	y_2=Lambda-Omega*x_2
	z_1=1/plane[2]*(plane[3]-plane[0]*x_1-plane[1]*y_1)
	z_2=1/plane[2]*(plane[3]-plane[0]*x_2-plane[1]*y_2)
	unit_1=[x_1,y_1,z_1]
	if separation(relVec1,unit_1) > 0.6 and separation(relVec2,unit_1) > 0.6:
		scale=2
		for j in range(3):
			unit_1[j]*=scale/math.sqrt(3)
		return unit_1
	else:
		unit_2=[x_2,y_2,z_2]
		for j in range(3):
			unit_2[j]*=scale/math.sqrt(3)
		return unit_2
def obtain_unit_vec_trigPlanar2(relVec1,relVec2,point):
  plane=find_plane(relVec1,relVec2,point)
  theta=angle_between_vecs(relVec1,relVec2)[1]
  J=math.pi-theta*1/2
  defVec=np.matrix([modulus(relVec1)*math.cos(J),modulus(relVec2)*math.cos(J),plane[3]], dtype=float)
 # J=2*math.pi-2*theta
 # defVec=np.matrix([modulus(relVec1)*math.cos(theta),modulus(relVec2)*math.cos(J),plane[3]], dtype=float)
  matrix=np.matrix([relVec1,relVec2,[plane[i] for i in range (3)]],dtype=float)
  matrix_inverse=np.linalg.inv(matrix)
  trigPlanVec=[0,0,0]
  for i in range (3):
    for j in range(3):
      trigPlanVec[i]+=matrix_inverse[i,j]*defVec[0,j]
  return trigPlanVec
	
