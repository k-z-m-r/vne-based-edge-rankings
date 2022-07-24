# def M_entropy(A):
#     L = np.diag(np.sum(A,axis=1)) - A
#     eigs, _ = np.linalg.eigh(L)
    
#     M = np.sum(A)/2
    
#     eigs = eigs[eigs > 0]
#     f = eigs/(2*M)
    
#     return -f*np.log2(f)


# def M_perturbation(A0,A1,e):
#     L0 = np.diag(np.sum(A0,axis=1)) - A0
#     L1 = np.diag(np.sum(A1,axis=1)) - A1
#     dL = L1 - L0
#     l,v = np.linalg.eigh(L0)
#     u = l[l > 0]
    
#     dh = np.log2(u) + 1/np.log(2)
    
#     M = np.sum(A0)/2
    
#     for i in range(len(l)):         
#         dh[i] *= (v[:,i].T @ dL @ v[:,i]) 
    
#     return -dh*e/(2*M)

# def Q_entropy(A0, A1, beta):

#     L0 = laplacian(A0); L1 = laplacian(A1)

#     eigs_0, _ = np.linalg.eigh(L0); eigs_1, _ = np.linalg.eigh(L1)

#     euler_0 = np.exp(-beta*eigs_0); euler_1 = np.exp(-beta*eigs_1)

#     q_i = euler_0 / np.sum(euler_0); p_i = euler_1 / np.sum(euler_1)

#     return p_i * np.log2(p_i / q_i)

# def Q_perturbation(A0, A1, epsilon, beta):

#     N = len(A0) #the number of nodes in the network

#     L0 = laplacian(A0); L1 = laplacian(A1) #initial and perturbed laplacians, respectively
#     dL = L1 - L0 #the change in laplacian

#     eigs, eigvecs = np.linalg.eigh(L0) #eigenvectors and eigenvalues according to the initial laplacian

#     f = np.exp(-beta*eigs)
#     f /= np.sum(f) #the function of e^(-\beta \lambda_i), f(\beta, \lambda).  We sometimes refer to it as u(\lambda)

#     dlambda = [eigvecs[:,j].T @ dL @ eigvecs[:,j] for j in range(N)] # d\lambda/d\epsilon, or how the eigenvectors change as \epsilon changes

#     dhdf = f/np.log(2) #dh/df, or dH/du, etc.  This is one part of the chain rule derivative
  
#     X = sum([-f[j] * dlambda[j] for j in range(N)])
#     dfdlambda = [dlambda[i] + X for i in range(N)] #this is the main results.  dhdf was already known basically from previous papers

#     return -(dhdf @ dfdlambda)*epsilon*beta