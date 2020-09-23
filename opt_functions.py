import numpy as np
import scipy.sparse as spr
from scipy.spatial import KDTree


def rmin_elnodes(ne,rmin,cent,me,mn,nt):
    yt = np.array([],dtype=np.uint32)
    yn = np.array([],dtype=np.uint32)
    rij = np.array([],dtype=np.float64)
    rijt = np.zeros(4,dtype=np.float64)
    for i in range(0,ne):
        x = cent[i][0]
        y = cent[i][1]
        z = cent[i][2]
        n1 = me[i][0]-1
        n2 = me[i][1]-1
        n3 = me[i][2]-1
        n4 = me[i][3]-1
        rijt[0] = np.sqrt((mn[n1][0]-x)**2+(mn[n1][1]-y)**2+(mn[n1][2]-z)**2)
        rijt[1] = np.sqrt((mn[n2][0]-x)**2+(mn[n2][1]-y)**2+(mn[n2][2]-z)**2)
        rijt[2] = np.sqrt((mn[n3][0]-x)**2+(mn[n3][1]-y)**2+(mn[n3][2]-z)**2)
        rijt[3] = np.sqrt((mn[n4][0]-x)**2+(mn[n4][1]-y)**2+(mn[n4][2]-z)**2)
        yt = np.concatenate((yt, n1, n2, n3, n4), axis=None)
        yn = np.concatenate((yn, np.repeat(i, nt)), axis=None)
        rij = np.concatenate((rij, rijt), axis=None)  
    yr = spr.csr_matrix((rij, (yt, yn))).toarray()
    return yr

def allocH1H2(ne,nauxIx,nauxIy,nauxIz,rmin,vifac,dim,t,cent):
    ii1 = np.array([],dtype=np.uint32)
    ii2 = np.array([],dtype=np.uint32)
    jj = np.array([],dtype=np.uint32)
    kk1 = np.array([],dtype=np.float64)
    kk2 = np.array([],dtype=np.float64)
    elemid = np.array([],dtype=np.uint32)
    rij = np.array([],dtype=np.float64)
    tree = KDTree(cent)
    if dim == 2:
        vmin = np.pi*rmin*rmin*t
    else:
        vmin = (4/3)*np.pi*rmin*rmin*rmin
    for i in range(0,ne):
        rij, elemid = tree.query(tree.data[i], 44, eps=0, p=2, distance_upper_bound=rmin)
        n = len(np.where(~np.isinf(rij))[0])
        rij = rij[:n]
        elemid = elemid[:n]
        sumf = 0
        for j in range(0,n):
            fac = rmin-rij[j]
            facn = vmin-vifac[elemid[j]][0]
            sumf = sumf+np.maximum(0,fac)*np.maximum(0,facn) # Update
            # sumf = sumf+np.maximum(0,fac) # 99line
            ii1 = np.concatenate((ii1, i), axis=None)
            jj = np.concatenate((jj, elemid[j]), axis=None)     
            kk1 = np.concatenate((kk1, np.maximum(0,fac)*np.maximum(0,facn) ), axis=None)  # Update
            # kk1 = np.concatenate((kk1, np.maximum(0,fac)), axis=None)  # 99line
        ii2 = np.concatenate((ii2, i), axis=None)
        kk2 = np.concatenate((kk2, sumf), axis=None)
    H1 = spr.csr_matrix((kk1, (ii1, jj))).toarray()
    H2 = spr.csr_matrix((kk2, (ii2, ii2))).toarray()
    return H1, H2

def allocH1H2H3(mn,ne,nn,nauxIx,nauxIy,nauxIz,rmin,vifac,dim,nt,cent,nod_inrad):
    ii = np.array([],dtype=np.uint32)
    jj = np.array([],dtype=np.uint32)
    kk = np.array([],dtype=np.float64)
    for i in range(0,nn):
        elemid = np.where(nod_inrad[i]!=0)[0]
        nn = len(elemid)
        rij = nod_inrad[i][elemid]
        for j in range(0,nn):
            fac2 = (1/((nn-1))*(1-(rij[j]/np.sum(rij))))
            ii = np.concatenate((ii, i), axis=None)
            jj = np.concatenate((jj, elemid[j]), axis=None)
            kk = np.concatenate((kk, fac2), axis=None)
    H1 = spr.csr_matrix((kk, (ii, jj))).toarray()
    ii1 = np.array([],dtype=np.uint32)
    ii2 = np.array([],dtype=np.uint32)
    jj1 = np.array([],dtype=np.uint32)
    kk1 = np.array([],dtype=np.float64)
    kk2 = np.array([],dtype=np.float64)
    elemid = np.array([],dtype=np.uint32)
    rij = np.array([],dtype=np.float64)
    tree = KDTree(mn)
    for i in range(0,ne):
        rij, elemid = tree.query(cent[i], 44, eps=0, p=2, distance_upper_bound=rmin)
        n = len(np.where(~np.isinf(rij))[0])
        rij = rij[:n]
        elemid = elemid[:n]
        sumf = 0
        for j in range(0,n):
            fac = rmin-rij[j]
            sumf = sumf+np.maximum(0,fac)
            ii1 = np.concatenate((ii1, i), axis=None)
            jj1 = np.concatenate((jj1, elemid[j]), axis=None)     
            kk1 = np.concatenate((kk1, fac), axis=None)
        ii2 = np.concatenate((ii2, i), axis=None)
        kk2 = np.concatenate((kk2, sumf), axis=None)
    H2 = spr.csr_matrix((kk1, (ii1, jj1))).toarray()
    H3 = spr.csr_matrix((kk2, (ii2, ii2))).toarray()
    return H1, H2, H3