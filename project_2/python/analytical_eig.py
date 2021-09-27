import numpy as np 

def analytical_eigval(n):
    
    h = 1.0/n
    a = -1.0/h**2
    d = 2.0/h**2
    
    eigval_list = [] 
    for i in range(1,n+1):
        
        eigval = d + 2*a*np.cos((i*np.pi)/(n+1))
        eigval_list.append(eigval)
    eigval_list = (eigval_list)/np.linalg.norm(eigval_list)
    return eigval_list
    
def analytical_eigvec(n):

    h = 1.0/n
    a = -1.0/h**2
    d = 2.0/h**2
    
    eigvec_list = [] 
    for i in range(1,n+1):
        eigvec = []
        for j in range(1,n+1):
            x = np.sin((j*i*np.pi)/(n+1))
            eigvec.append(x)
        eigvec_list.append(eigvec)
    eigvec_list = eigvec_list / np.linalg.norm(eigvec_list,axis=1)
    return eigvec_list

