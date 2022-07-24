import numpy as np 
from matplotlib import pyplot as plt
import networkx as nx

from graphs_util import *


def B_entropy(A, beta):
    
    L = laplacian(A)
    evals, _ = np.linalg.eigh(L)
    beta_evals = -beta*evals
    
    z = np.exp(beta_evals)

    Z = np.sum(z)
    f = z/Z

    H = -f*np.log2(f)

    return H

def B_perturbation(dL, eigs, eigvecs, f, dhdf, epsilon, beta):

    dlambda = np.diag(eigvecs.T @ dL @ eigvecs) # d\lambda/d\epsilon, or how the eigenvectors change as \epsilon changes
  
    X = -f * dlambda
    X = X.sum()

    dfdlambda = dlambda + X

    return -(dhdf @ dfdlambda)*epsilon*beta

def B_perturbation_1_edge(dL, i, j, f, dhdf, eigvecs, epsilon, beta):
    
    dlambda = -(eigvecs[i] - eigvecs[j])**2
    
    dfdlambda = (-f * dlambda).sum() + dlambda
    
    return -dL*np.sum(dhdf * dfdlambda)*epsilon*beta


def edge_rankings(A, beta, epsilon):
    
    A0 = A.copy()
    L0 = laplacian(A0)
    h0 = B_entropy(A0, beta).sum()
    eigs, vecs = np.linalg.eigh(L0)
    
    f = np.exp(-beta*eigs)
    f /= f.sum()
    dhdf = f*(-np.log2(f) + 1/np.log(2))
    
    edges = get_edges(A0)
    
    M = edges.shape[0]
    Hs = np.zeros((2, M))
    for idx, edge in enumerate(edges):
        weights = remove_edge(A0, edge)
        i, j = edge
        Hs[0, idx] = B_entropy(epsilon*A0 + (1 - epsilon)*A, beta).sum() - h0
        Hs[1, idx] = B_perturbation_1_edge(weights[0], i, j, f, dhdf, vecs, epsilon, beta).sum()
        add_edge(A0, edge, weights)
        
        print('{:.2f}% done!'.format((idx + 1)/M * 100), '\r', end = '')
    
    edge_sort = np.argsort(-Hs, axis = 1)
    edge_ranks = np.argsort(edge_sort, axis = 1)
    
    return Hs, edge_sort, edge_ranks

def modified_edge_rankings(A, beta, epsilon, prints = False):
    
    A0 = A.copy()
    L0 = laplacian(A0)
    h0 = B_entropy(A0, beta).sum()
    eigs, vecs = np.linalg.eigh(L0)
    
    f = np.exp(-beta*eigs)
    f /= f.sum()
    dhdf = f*(-np.log2(f) + 1/np.log(2))
    
    edges = get_edges(A)
    
    M = edges.shape[0]
    Hs = np.zeros((M))
    for idx, edge in enumerate(edges):
        weights = remove_edge(A, edge)
        i, j = edge
        Hs[idx] = B_perturbation_1_edge(weights[0], i, j, f, dhdf, vecs, epsilon, beta).sum()
        add_edge(A, edge, weights)
        if prints == True:
            print('{:.2f}% done!'.format((idx + 1)/M * 100), '\r', end = '')
    
    edge_sort = np.argsort(-Hs)
    edge_ranks = np.argsort(edge_sort)
    
    return Hs, edge_sort, edge_ranks

def transit_edge_rankings(A, eigs, vecs, beta, epsilon):
    
    A0 = A.copy()
    
    f = np.exp(-beta*eigs)
    f /= f.sum()
    dhdf = f*(-np.log2(f) + 1/np.log(2))
    
    edges = get_edges(A0)
    
    M = edges.shape[0]
    Hs = np.zeros((M))
    for idx, edge in enumerate(edges):
        weights = remove_edge(A0, edge)
        i, j = edge
        Hs[idx] = B_perturbation_1_edge(weights[0], i, j, f, dhdf, vecs, epsilon, beta).sum()
        add_edge(A0, edge, weights)
        
        print('{:.2f}% done!'.format((idx + 1)/M * 100), '\r', end = '')
    
    edge_sort = np.argsort(-Hs)
    edge_ranks = np.argsort(edge_sort)
    
    return Hs, edge_sort, edge_ranks











