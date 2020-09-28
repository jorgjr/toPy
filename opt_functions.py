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
        yt = np.append(yt, np.array([n1, n2, n3, n4]))
        yn = np.append(yn, np.repeat(i, nt))
        rij = np.append(rij, rijt)  
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
    vr = np.uint32(vmin/np.min(vifac))
    for i in range(0,ne):
        rij, elemid = tree.query(tree.data[i], vr, eps=0, p=2, distance_upper_bound=rmin)
        n = len(np.where(~np.isinf(rij))[0])
        rij = rij[1:n]
        elemid = elemid[1:n]
        fac = np.reshape(rmin-rij,(n-1,1))
        sumf = np.sum(fac)
        ii1 = np.append(ii1, np.repeat(i,n-1))
        jj = np.append(jj, elemid)
        kk1 = np.append(kk1, fac)
        ii2 = np.append(ii2, i)
        kk2 = np.append(kk2, sumf)
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
        fac2 = (1/((nn-1))*(1-(rij/np.sum(rij))))
        ii = np.append(ii, np.repeat(i,nn))
        jj = np.append(jj, elemid)
        kk = np.append(kk, fac2)
    H1 = spr.csr_matrix((kk, (ii, jj))).toarray()
    ii1 = np.array([],dtype=np.uint32)
    ii2 = np.array([],dtype=np.uint32)
    jj1 = np.array([],dtype=np.uint32)
    kk1 = np.array([],dtype=np.float64)
    kk2 = np.array([],dtype=np.float64)
    elemid = np.array([],dtype=np.uint32)
    rij = np.array([],dtype=np.float64)
    tree = KDTree(mn)
    vmin = (4/3)*np.pi*rmin*rmin*rmin
    vr = np.uint32(vmin/np.min(vifac))
    for i in range(0,ne):
        rij, elemid = tree.query(cent[i], vr, eps=0, p=2, distance_upper_bound=rmin)
        n = len(np.where(~np.isinf(rij))[0])
        rij = rij[1:n]
        elemid = elemid[1:n]
        fac = rmin-rij
        sumf = np.sum(fac)
        ii1 = np.append(ii1, np.repeat(i,n-1))
        jj1 = np.append(jj1, elemid)     
        kk1 = np.append(kk1, fac)
        ii2 = np.append(ii2, i)
        kk2 = np.append(kk2, sumf)
    H2 = spr.csr_matrix((kk1, (ii1, jj1))).toarray()
    H3 = spr.csr_matrix((kk2, (ii2, ii2))).toarray()
    return H1, H2, H3