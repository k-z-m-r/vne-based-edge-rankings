############################################################
#  Jeremy Kazimer and Dane Taylor
#  Basic Graph Utility
#  5/24/21
#  Version 4
#############################################################

import numpy as np 

#Returns the max edges of a network for a given N.  Assuming unweighted.
def max_edges(n,self_edges = False):
    if (self_edges == False):
        return n*(n - 1)/2
    else:
        return n*(n-1)/2 + n

#Returns the Unnormalized Laplacian of a given Adjacency Matrix.
def laplacian(A):
    return np.diag(A.sum(axis = 1)) - A

#Forms an Adjacency Matrix out of an Edge List.
def form_A(edge_list, N):
    A = np.zeros((N, N))
    for edge in edge_list:
        A[edge[0], edge[1]] = 1

    A = A + A.T

    return A

def get_edge_type(in_hot, edges):
   
    in_mask = in_hot @ in_hot.T
    
    locs = np.array([])
    
    for edge in edges:
        i, j = edge
        if in_mask[i, j] == 1:
            locs  = np.append(locs , 1)
        else:
            locs  = np.append(locs, 0)
    
    return locs

#Forms Adjacency Matrix given the number of nodes (n) and edge probability (p).
def Erdos_Renyi_Gnp(n, p, self_edges = False):
    
    B = np.array(np.random.rand(n,n)<p,dtype=int)
    B = np.triu(B)

    if self_edges == True:
        A = B + B.T - np.diag(np.diag(B))
    else:
        A = B + B.T - 2*np.diag(np.diag(B))
        
    return A

def Erdos_Renyi_Gnm(N, M):
    A = np.zeros((N, N))

    X = np.ones((N, N))
    X = np.triu(X)
    np.fill_diagonal(X, 0)
    
    edges = np.random.permutation(np.argwhere(X == 1))[:M]
    
    for edge in edges:
        i, j = edge
        
        A[i, j] = 1
        A[j, i] = 1
    
    return A

def SBM(n, p_in, p_out, k, self_edges = False):
    N = n//k

    one_hot = np.zeros((n, k))

    for dim in range(k):
        one_hot[(N*dim):(N*(dim + 1)), dim] = 1

    in_mask = one_hot @ one_hot.T
    out_mask = 1 - in_mask

    A_in = Erdos_Renyi_Gnp(n, p_in) * in_mask
    A_out = Erdos_Renyi_Gnp(n, p_out) * out_mask

    A = A_in + A_out

    return A, one_hot



def Multilayer(N, M, k, Dx):
    
    n = N*k
    A = Erdos_Renyi_Gnm(n, M)
    
    one_hot = np.zeros((k, n)).T
    for dim in range(k):
        one_hot[(N*dim):(N*(dim + 1)), dim] = 1
        
    inmask = one_hot @ one_hot.T
    
    A = A*inmask
 
    N_ = (k - 1)*N
       
    A = A + (Dx*(np.diag(np.ones(N_), N) + np.diag(np.ones(N_), -N))).astype(float)
    
    return A, one_hot
             
def get_edges(A):
    return np.argwhere(np.triu(A) != 0)


############################################################
# Rewiring functions
############################################################

def remove_edge(A, edge):
    i, j = edge
    
    A_ij, A_ji = A[i, j], A[j, i]
    
    A[i, j] = 0
    A[j, i] = 0
    
    return A_ij, A_ji

def add_edge(A, edge, weights = (1, 1)):
    i, j = edge
    
    A_ij, A_ji = weights
    
    A[i, j] = A_ij
    A[j, i] = A_ji
    
    return None

def rewire_graph(A, r_edge, a_edge):
    
    weights = remove_edge(A, r_edge)
    add_edge(A, a_edge, weights)
    
    return A
    
    
    























