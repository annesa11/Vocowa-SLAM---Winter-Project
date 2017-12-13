import matplotlib.pyplot as plt
import numpy as np
from matplotlib import style
import math
from random import randrange
from scipy import stats
from sklearn.linear_model import RANSACRegressor

def ransac_fit_line(data,sigma=1,iterations=100,inliers=10):# In this code i am specifying the minimum number of inliers that needs to be present
	#data is 2*n list                                                    
	X=np.array(data)
	xx=X[:,0]
	yy=X[:,1]
	plt.scatter(xx,yy,marker='+',color='red')
	xx=xx.reshape(-1,1)
	yy=yy.reshape(-1,1)
	
	ransac=RANSACRegressor(max_trials=iterations) # if you want to change the residual_threshold add that attribute
	ransac.fit(xx,yy)
	inlier_mask=ransac.inlier_mask_
	outlier_mask=np.logical_not(inlier_mask)
	X_selected=xx[inlier_mask]
	xx1=min(X_selected)
	xx2=max(X_selected)
	line_X=np.arange(X_selected.min(),X_selected.max()+1)[:,np.newaxis]
	line_y_ransac=ransac.predict(line_X)
	yy1=ransac.predict(xx1)[0]
	yy2=ransac.predict(xx2)[0]
	Z=[xx1[0],(xx1[0]+xx2[0])/2,xx2[0]]
	K=[yy1[0],(yy1[0]+yy2[0])/2,yy2[0]]
	
	plt.plot(line_X,line_y_ransac,color='cornflowerblue')
	plt.scatter(Z,K,marker='+',color='green')
	plt.show()
	uandme=np.size(np.where(inlier_mask==True))
	
	if uandme>=inliers:
		return(Z,K,inlier_mask)
	else:
		Z1=[0,0,0]
		inlier_mask_false=np.zeros((np.size(inlier_mask),1),dtype=bool)
		return(Z1,Z1,inlier_mask_false) # If the condition is based on no of inliers select this part of code
	
	


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
	
data=[[0,0],[20,100],[70,0],[100,80]]#,[1.75,1.75],[3,3],[4,1],[2.25,2.25],[2.5,2.5],[2.75,2.75],[3.25,3.25],[3.5,3.5],[5,3],[3,5]]	
X,Y,Z=ransac_fit_line(data)
print(X)
print(Y)
print(Z)
