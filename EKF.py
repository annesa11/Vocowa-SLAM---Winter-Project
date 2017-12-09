import numpy as np
import math
from random import randrange
         
def toFrame(F,p,nargout=1):
     # F is basically the reference frame,p is the point in global frame
	 # F=[F_x;F_y;F_alpha]
     np.matrix(F)
     np.matrix(p)
     t=F[0:2]
     a=F[2][0]
     R=[[math.cos(a),-math.sin(a)],[math.sin(a),math.cos(a)]]
     np.matrix(R)
     k=np.matrix(np.subtract(p,t))
     jf=np.transpose(R)*k
     pf=jf.tolist()
     	 
     if nargout>1:
       px=p[0][0]
       py=p[1][0]
       x=t[0][0]
       y=t[1][0]
       Pf_f=[[-math.cos(a),-math.sin(a),math.cos(a)*(py-y)-math.sin(a)*(px-x)],[math.sin(a),-math.cos(a),-math.cos(a)*(px-x)-math.sin(a)*(py-y)]]
       Pf_p=np.transpose(R)		
       Pf_p=Pf_p.tolist()
       return([pf,Pf_f,Pf_p])
     else:
       return([pf])	

	   

def fromFrame(F,pf,nargout=1):
	# F is basically the reference frame, pf is the point int the reference frame F_alpha
	# We have to basically pf into global frame
	# # F=[F_x;F_y;F_alpha]
	np.matrix(F)
	p1=np.matrix(pf) #Basically this part involves conversion of pf to matrix form multiplication
	t=F[0:2]
	a=F[2][0]
	R=[[math.cos(a),-math.sin(a)],[math.sin(a),math.cos(a)]]
	np.matrix(R)
	pw=np.add(R*p1,t) 
	pw=pw.tolist()
	if nargout>1:
	   px=pf[0][0]
	   py=pf[1][0]
	   
	   Pw_f=[[1,0,-py*math.cos(a)-px*math.sin(a)],[0,1,px*math.cos(a)-py*math.sin(a)]]
	   Pw_pf=R;
	   return([pw,Pw_f,Pw_pf])
	else:
	   return([pw])
	

def scan(p,nargout=1):
	# This function basically performs a range and bearing measurement of a 2D point.
	# p=[p_x;p_y] points in sensor frame
	
	px=p[0][0]
	py=p[1][0]
	j=px**2+py**2
	d=math.sqrt(j)
	a=math.atan(py/px)
	y=[[d],[a]]
	if nargout>1:
		Y_p=[[(px/d),(py/d)],[-(py/j),(px/j)]]
		return([y,Y_p])
	else:
		return([y])

def invscan(y,nargout=1):
	# This function is basically used to backproject a range and bearring measurement into a 2-d point
	# y=[range,bearing]
    d=y[0][0]
    a=y[1][0]
    px=d*math.cos(a)
    py=d*math.sin(a)
    p=[[px],[py]]
	
    if nargout>1:
        P_y=[[math.cos(a),-py],[math.sin(a),px]]
        return([p,P_y])
    else:
        return([p])
	

def move(r,u,n,nargout=1):
    #r=[x;y;alpha]->This is the robot pose 
    #u=[dx;dalpha]->This is the control input
    #n=[nx;nalpha]->This is the perturbation input
    #ro->updated robot pose	
    # Basically the control inputs are provided with respect to the robot frame which is then being converted to World Frame.
    a=r[2][0]
    
    dx=u[0][0]+n[0][0]
    da=u[1][0]+n[1][0]
    ao=a+da
    
    if ao>(math.pi):
        ao=ao-2*(math.pi)
	
    if ao<-(math.pi):
        ao=ao+2*(math.pi)
	
    dp=[[dx],[0]]
	
    if nargout==1:
        to=fromFrame(r,dp) # to is a list so ans is obtained using to[0]
        ro=[to[0][0],to[0][1],[ao]]		
        return([ro])
    else:
        po=fromFrame(r,dp,3)
        Ro_r=[po[1][0],po[1][1],[0,0,1]] # this is Jacobian wrt r		
        Ro_n=[[po[2][0][0],0],[po[2][1][0],0],[0,1]] # This is Jacobian wrt n		
        ro=[po[0][0],po[0][1],[ao]] # This is the point which contains the updated robot pose
        return([ro,Ro_r,Ro_n])


def observe(r,p,nargout=1):
    # This function converts the point p into global frame and then it uses to obtain Range and Bearing
    # r=[ro_x,ro_y,ro_alpha]-> Robot pose
    # p=[p_x,p_y]-> point in global frame
    # y=[range,bearing]	
    if nargout==1:
        y0=toFrame(r,p)		
        y=scan(y0[0])
        return(y)		
    else:
        J=toFrame(r,p,3)
        Pr_r=np.matrix(J[1])# Converted to matrix form for the purpose of multiplication		
        Pr_p=np.matrix(J[2]) # Converted to matrix for the purpose of multiplication
        K=scan(J[0],2)		
        Y_pr=np.matrix(K[1]) #Converted to matrix form for the purpose of multiplication
        Y_r=Y_pr*Pr_r # Jacobian wrt r
        Y_p=Y_pr*Pr_p # Jacobian wrt p				
        return([K[0],Y_r,Y_p])


def invobserve(r,y,nargout=1):
    # Backprojects a range and bearing measurement and transport to Map Frame 
    # y=[range;bearing]
    # r=[ro_x,ro_y,ro_alpha]-> Robot Pose
	
    if nargout==1:
        kim=invscan(y)	
        p=fromFrame(r,kim[0])
        return(p)		
    else:
        kim=invscan(y,2)
        p_r=kim[0]
        Pr_y=np.matrix(kim[1])
        p0=fromFrame(r,p_r,3)
        p=p0[0]
        P_r=p0[1]
        P_pr=np.matrix(p0[2])
        P_y=P_pr*Pr_y
        P_y=P_y.tolist() 
        return([p,P_r,P_y])		


#W=cloister(-4,4,-4,4,7) # Set of external landmarks of the form 2*N
W=[[1,2,3,4,5,6],[1,2,3,4,5,6]]
k=np.shape(W)
N=k[1]  # Total No of landmarks
print(N)
R=[[0],[-2],[0]] # Robot pose
print(R)
U=[[0.1],[0.05]] # Control parameter
print(U)
Y=np.zeros((2,N)) # Set of landmarks measurements
print(Y)

# Estimator
x=np.zeros(((3+2*N),1)) #State Vectors Mean
print(x)
P=np.zeros((3+2*N,3+2*N)) # Covariance Matrix 3+2*N by 3+2*N type of array
print(P)
q=[[0.01],[0.02]] # System noise ampli
print(q)
Q=[[0.01**2,0],[0,0.02**2]] # System noise covariance Matrix
print(Q)
s=[[0.1],[math.pi/180]] # Measurement noise ampli
print(s)
S=[[0.1**2,0],[0,(math.pi/180)**2]] # Measurement noise covariance Matrix
print(S)
mapspace=np.zeros((1,3+2*N),dtype=bool) 
landmarks=np.zeros((2,N))
r=np.where(mapspace==False)
r=r[1][0:3]

print(r)
mapspace[0][r]=True
print(mapspace)
for i in range(3):
	x[i]=R[i]
	for j in range(3):
		P[i][j]=0

print(x)
		
for t in range(1):
	q0=np.random.randn(2,1)
	n=np.multiply(q,q0)  # Perturbation Vector
	print(n)
	zeros=[[0],[0]]
	R=move(R,U,zeros) # We will perturb the estimator instead of the stimulator
	R=R[0] # THE ANS IS RETURN IN THE FORM OF LIST
	print(R)
	
	for i in range(N):
		s0=np.random.randn(2,1) 
		v=np.multiply(s,s0)  #measurement noise
		print(v)
		W1=[row[i] for row in W]
		W1=np.reshape(W1,(2,1)) #reshaping back in the form 2*1
		print(W1)
		Y0=observe(R,W1) 
		print(Y0[0])
		Y1=Y0[0]+v  #Adding the measurement noise to the measurement of landmark
		print(Y1)
		j=0
		for row in Y:
			row[i]=Y1[j][0]
			j=j+1
	print(Y)
	
	m0=np.where(np.transpose(landmarks)!=0) 
	o=np.transpose(landmarks)[m0]  #pointers to landmark in the form of 1*e
	m=np.reshape(o,(np.shape(o)[0],1)) #pointers to landmark in the form of e*1
	rm=np.union1d(r,o)  # all used states
	rm=rm.astype(int) # converting to int type for index usage
	delta=x[r]
	Estim_r=move(delta,U,n,3) #Estimator perturbated
	x[r]=Estim_r[0]
	
	
	A1=np.matrix(Estim_r[1])
	A2=np.matrix(P[np.ix_(r,o)])
	A3=np.matrix(P[np.ix_(r,r)])
	A4=np.matrix(Q)
	A5=np.matrix(Estim_r[2])
	P[np.ix_(r,o)]=np.array(A1*A2)
	P[np.ix_(o,r)]=np.transpose(P[np.ix_(r,o)])
	P[np.ix_(r,r)]=A1*A3*np.transpose(A1)+A5*A4*np.transpose(A5)  # The above part involves updation of Covariance Matrix
	
	G0=np.where(landmarks[1]!=0)
	G=G0[0] # Return landmarks pointer
	for i in G:
		l=landmarks[:,i]
		l=l.astype(int)
		E=observe(x[r],x[l],3) # This part is equivalent to h(x) in EKF slam
		rl=np.union1d(r,l) #Taking union of r and l
		rl=rl.astype(int) # converting to int type for index usage
		e=E[0]
		E_r=E[1]
		E_l=E[2]
		E_rl=[[E_r[0][0],E_r[0][1],E_r[0][2],E_l[0][0],E_l[0][1]],[E_r[1][0],E_r[1][1],E_r[1][2],E_l[1][0],E_l[1][1]]]
		A1=np.matrix(E_rl)
		A2=np.matrix(P[np.ix_(rl,rl)])
		Ealpha=A1*A2*np.transpose(A1)
		
		Yi=Y[:,1] # Measurement of landmark i
		Yi=np.reshape(Yi,(2,1))
		z=np.subtract(Yi-e) # innovation Gaussian(z,Z)
		if z[1][0]>(math.pi):  #We need values around zero for angles
			z[1][0]=z[1][0]-2*(math.pi)
		else:
			z[1][0]=z[1][0]+2*(math.pi)
		
		Z=S+Ealpha
		
		A3=np.matrix(Z)
		A4=np.matrix(z)
		A5=np.linalg.inv(A1)
		A6=np.matrix(P[np.ix_(rm,rl)])
		A8=np.matrix(P[np.ix_(rm,rm)])
		Qw=np.transpose(A4)*A5*A4
		
		if Qw<9:
			K=A6*np.transpose(A1)*A5 # Kalman Gain
			A7=np.matrix(K)
			# Map updation using rm
			X[rm]=X[rm]+A7*A4    
			P[np.ix_(rm,rm)]=A8-(A7*A3*np.transpose(A7))
	
	
	# Landmark initialization-one new landmark at each iteration	
	lids=np.where(landmarks[0]==0) # all non initialized landmarks
	if len(lids[0])!=0: # there are still landmarks to initialize
		print(len(lids[0]))
		rand_index=randrange(0,len(lids[0]))
		i=lids[0][rand_index] #pick one landmark randomly
		print(i)
		l=np.where(mapspace==False)
		l=l[1][0:2]
		if len(l)!=0:
			mapspace[0][l]=True
			landmarks[:,i]=l
			Yi=Y[:,i]
			Yi=np.reshape(Yi,(2,1))
			Landalpha=invobserve(x[r],Yi,3)
			x[l]=Landalpha[0]
			L_r=Landalpha[1]
			L_y=Landalpha[2]
			
			A1=np.matrix(L_r)
			A2=np.matrix(L_y)
			A3=np.matrix(P[np.ix_(r,rm)])
			A4=np.matrix(P[np.ix_(r,r)])
			A5=np.matrix(S)
			P[np.ix_(l,rm)]=A1*A3
			P[np.ix_(rm,l)]=np.transpose(P[np.ix_(l,rm)])
			P[np.ix_(l,l)]=A1*A4*np.transpose(A1)+A2*A5*np.transpose(A2)
			print(P)
			print(x)
			print(mapspace)
			print(landmarks)