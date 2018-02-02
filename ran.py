import matplotlib.pyplot as plt
import numpy as np
from matplotlib import style
import math
from random import randrange

style.use('ggplot')

def ransac_fit_line(data,iterations=25,sigma=1,inliers=5):
	#data is 2*n list
	X=np.array(data)
	xx=X[:,0]
	yy=X[:,1] 
	plt.scatter(xx,yy,marker='+',color='red') #plotting the points
	sze=np.size(xx) # sixe of data
	Line1=0 # Best fitting line with largest number of inliers
	p1=0
	p2=0
	final_inliers=[]
	for i in range(iterations):
		k=np.random.permutation(sze)
		h=k[0:2]
		sample=X[h]
		linearray=sample[0]-sample[1]
		norm_linearray=np.linalg.norm(linearray)
		J=linearray/norm_linearray
		norm=np.matrix([[-J[1]],[J[0]]])
		distances=np.zeros(sze)
		selected=[]
		selected_count=0
		selected_inliers=[]
		
		for j in range(sze):
			V=np.matrix(X[j]-sample[0])
			#print(V)
			#print(norm)
			
			distances[j]=V*norm
			
			if distances[j]<=sigma:
				selected_count+=1
				selected_inliers.append(X[j])
				
				
			
		if selected_count>=inliers and selected_count>=Line1:
			Line1=selected_count
			final_inliers=np.array(selected_inliers)
			print(sample[0][0],sample[0][1],sample[1][1],sample[1][0])
			parameter1=(sample[1][1]-sample[0][1])/(sample[1][0]-sample[0][0])
			print(parameter1)
			parameter2=sample[0][1]-parameter1*sample[0][0]
			p1=parameter1
			p2=parameter2
		
	print(p1)
	print(p2)
	#print(final_inliers[:,1])
	#print(final_inliers[:,0])
	#plt.plot(final_inliers[:,0],final_inliers[:,1],'-y')
	
	#return(final_inliers)	
			
				
global NO_OF_TRAILS=200
global SEGMENT = 10
global GAP = 2/3
def ransac(raw_data):
	data=scan(raw_data)
	partitionLength=int((len(data)/360)*SEGMENT)
	lines=[]
	Visited=[0 for i in range(len(data))]
	for i in range(NO_OF_TRAILS):
		center = random.randInt(partitionLength/2,len(data)-partitionLength/2)
		if Visited[center]==0:
			begin=center-partitionLength/2
			end=center+partitionLength/2
			X,Y,TF=ransac_fit_line(data[begin:end])
			if (TF)==0:  # Condition if line is accepted
				for j in range(begin,end):
					if TF[j-begin]:
						Visited[j]=1
				lines.append([X,Y])
	return lines

data=[[1,1],[2,2],[1.25,1.25],[1.5,1.5],[1.75,1.75],[3,3],[4,1],[2.25,2.25],[2.5,2.5],[2.75,2.75],[3.25,3.25],[3.5,3.5],[5,3],[3,5]]	
ransac_fit_line(data)