import math
import numpy as np

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
        P_r=np.matrix(p0[1])
        P_pr=np.matrix(p0[2])
        P_y=P_pr*Pr_y
        P_y=P_y.tolist() 
        return([p,P_y])		
		
F=[[1],[2],[1.04]]
p=[[5],[6]]
n=[[0],[0]]
J=toFrame(F,p,3)
oi=J[0]
#print(J[0])
#print(J[1])
#print(J[2])
L=fromFrame(F,J[0],3)
#print(L[0])
#print(L[1])
#print(L[2])

k=[[1],[2]]
U=scan(k,2)
#print(U[0])
h=invscan(U[0],2)
#print(h[0])
#print(h[1])
io=move(F,p,n,3)
gg=observe(F,[],3)
print(gg[0])
#print(gg[1])
#print(gg[2])
pp=invobserve(F,gg[0],3)
#print(pp)
#print(pp[0])
#print(pp[1])