import numpy as np

def get_voting_array(df, drop = True):
    
    '''
        df -> the congressional dataframe for a particular session
        drop -> remove people who don't vote, a flag of Yes or No
    '''
    
    ids = df['icpsr']
    unq_ids = np.unique(ids)
    
    bills = df['rollnumber']
    unq_bills = np.unique(bills)
    
    V = np.zeros((unq_ids.shape[0], unq_bills.shape[0]))
    
    to_drop = []
    
    for idx, c_id in enumerate(unq_ids):
        
        c_df = df[df['icpsr'] == c_id]
        c_bills = np.array((c_df['rollnumber'])) - 1
        c_cast = c_df['cast_code']
        
        V[idx, c_bills] = c_cast
        
        if np.sum(c_cast) == 0:
            to_drop.append(idx)

    if drop == True:
        V = np.delete(V, to_drop, axis = 0)
        
    return V, to_drop

def get_one_hot(votes_df, ideology_df, to_drop = None, drop = True):
    
    '''
        votes_df -> the dataframe encoding votes that is used to get V.
        ideology_df -> this stores party information.
    '''

    ids = ideology_df['icpsr']
    
    p_ids = np.zeros((ids.shape[0]))
    names = np.array([])
    
    for idx, c_id in enumerate(np.unique(votes_df['icpsr'])):
        X = ideology_df[c_id == ids]
        
        p_ids[idx] = X['party_code']
        names = np.append(names, X['bioname'].astype(str))

        
    if drop == True and to_drop is not None:
        p_ids = np.delete(p_ids, to_drop, axis = 0)
        names = np.delete(names, to_drop, axis = 0)
        
    unq_parties = np.unique(p_ids)
    
    one_hot = np.zeros((p_ids.shape[0], unq_parties.shape[0]))
    
    for idx, p in enumerate(unq_parties):
        
        one_hot[:, idx] = (p == p_ids).astype(int)
    
    return one_hot, p_ids, names

def get_adjacency(V, p_ids, one_hot, names):
    
    '''
        V -> the array containing voting trends.
    '''
    
    N = V.shape[0]
    A = np.zeros((N, N))
    
    for idx, v_0 in enumerate(V):
        for jdx, v_1 in enumerate(V):

            if idx != jdx:
                locs = v_0 == v_1
                sim = np.sum(locs)/locs.shape[0]
                A[idx, jdx] = sim
                
    locs = np.sum(A, axis = 1) > 10
    
    A = A[locs].T[locs]
    p_ids = p_ids[locs]
    names = names[locs]
    one_hot = one_hot[locs]
    
    return A, p_ids, one_hot, names