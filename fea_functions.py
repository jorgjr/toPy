import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

    
# Common Functions

def delete_from_csr(mat, row_indices=[], col_indices=[]):
    if not isinstance(mat, sp.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    rows = []
    cols = []
    if row_indices:
        rows = list(row_indices)
    if col_indices:
        cols = list(col_indices)
    if len(rows) > 0 and len(cols) > 0:
        row_mask = np.ones(mat.shape[0], dtype=bool)
        row_mask[rows] = False
        col_mask = np.ones(mat.shape[1], dtype=bool)
        col_mask[cols] = False
        return mat[row_mask][:,col_mask]
    elif len(rows) > 0:
        mask = np.ones(mat.shape[0], dtype=bool)
        mask[rows] = False
        return mat[mask]
    elif len(cols) > 0:
        mask = np.ones(mat.shape[1], dtype=bool)
        mask[cols] = False
        return mat[:,mask]
    else:
        return mat

def SolverFEM(Ut,sK,sF,U):
    ff = sF
    uu = Ut
    auxu = np.where(Ut==0)
    auxi = np.where(np.isnan(Ut))
    ff = np.delete(ff,auxu[0],0)
    uu = np.delete(uu,auxu[0],0)
    kk = delete_from_csr(sK, auxu[0].tolist(), auxu[0].tolist())
    # uu, _ = spl.minres(kk,ff,x0=None, shift=0.0, tol=1e-05, maxiter=1000, M=None, callback=None, show=False, check=False)
    uu = spl.spsolve(kk,ff)
    for i in range(0,len(auxi[0])):
        U[auxi[0][i]] = uu[i]
    F = np.asmatrix(sK.toarray())*np.asmatrix(U)
    return np.asarray(F), U

def SolverFEM_opt(Ut,sK,sF,U):
    ff = sF
    uu = Ut
    auxu = np.where(Ut==0)
    auxi = np.where(np.isnan(Ut))
    ff = np.delete(ff,auxu[0],0)
    uu = np.delete(uu,auxu[0],0)
    kk = delete_from_csr(sK, auxu[0].tolist(), auxu[0].tolist())
    # uu, _ = spl.minres(kk,ff,x0=None, shift=0.0, tol=1e-05, maxiter=1000, M=None, callback=None, show=False, check=False)
    uu = spl.spsolve(kk,ff)
    for i in range(0,len(auxi[0])):
        U[auxi[0][i]] = uu[i]
    return U

def TetrahedronElementPStresses(sigma):
    s1 = sigma[0][0]+sigma[1][0]+sigma[2][0]
    s2 = sigma[0][0]*sigma[1][0]+sigma[0][0]*sigma[2][0]+sigma[1][0]*sigma[2][0]- \
         sigma[3][0]*sigma[3][0]-sigma[4][0]*sigma[4][0]-sigma[5][0]*sigma[5][0]
    ms3 = np.array([[sigma[0][0],sigma[3][0],sigma[5][0]], \
                    [sigma[3][0],sigma[1][0],sigma[4][0]], \
                    [sigma[5][0],sigma[4][0],sigma[2][0]]],dtype=np.float64)
    s3 = np.linalg.det(ms3)
    y = np.array([[s1],[s2],[s3]],dtype=np.float64)
    return y


# 2D Functions

def LinearTriangleElementStiffness(A,D,t,xi,yi,xj,yj,xm,ym):
    betai = yj-ym
    betaj = ym-yi
    betam = yi-yj
    gammai = xm-xj
    gammaj = xi-xm
    gammam = xj-xi
    aux = np.asarray([[betai,0,betaj,0,betam,0],[0,gammai,0,gammaj,0,gammam], \
                      [gammai,betai,gammaj,betaj,gammam,betam]],dtype=np.float64)
    aux = aux/(2*A)
    B = np.asmatrix(aux)
    y = (t*A[0])*B.T*D*B
    return y, B

def LinearTriangleAssemble(k,me,ne):
    row = np.reshape(np.hstack((np.repeat(2*np.reshape(me,(-1,1))-2,6,axis=1), \
                                np.repeat(2*np.reshape(me,(-1,1))-1,6,axis=1))),(1,-1))
    row = np.ravel(row)
    col = np.reshape(np.hstack((2*np.reshape(me,(-1,1))-2, \
                                2*np.reshape(me,(-1,1))-1)),(ne,6))
    col = np.hstack((col,col,col,col,col,col))
    col = np.ravel(col)
    data = np.ravel(k)
    y = sp.csr_matrix((data, (row, col)))
    return y


# 3D Functions

def TetrahedronElementStiffness(V,D,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4):
    mbeta1 = np.asarray([[1,y2,z2],[1,y3,z3],[1,y4,z4]],dtype=np.float64)
    mbeta2 = np.asarray([[1,y1,z1],[1,y3,z3],[1,y4,z4]],dtype=np.float64)
    mbeta3 = np.asarray([[1,y1,z1],[1,y2,z2],[1,y4,z4]],dtype=np.float64)
    mbeta4 = np.asarray([[1,y1,z1],[1,y2,z2],[1,y3,z3]],dtype=np.float64)
    mgamma1 = np.asarray([[1,x2,z2],[1,x3,z3],[1,x4,z4]],dtype=np.float64)
    mgamma2 = np.asarray([[1,x1,z1],[1,x3,z3],[1,x4,z4]],dtype=np.float64)
    mgamma3 = np.asarray([[1,x1,z1],[1,x2,z2],[1,x4,z4]],dtype=np.float64)
    mgamma4 = np.asarray([[1,x1,z1],[1,x2,z2],[1,x3,z3]],dtype=np.float64)
    mdelta1 = np.asarray([[1,x2,y2],[1,x3,y3],[1,x4,y4]],dtype=np.float64)
    mdelta2 = np.asarray([[1,x1,y1],[1,x3,y3],[1,x4,y4]],dtype=np.float64)
    mdelta3 = np.asarray([[1,x1,y1],[1,x2,y2],[1,x4,y4]],dtype=np.float64)
    mdelta4 = np.asarray([[1,x1,y1],[1,x2,y2],[1,x3,y3]],dtype=np.float64)
    beta1 = -1*np.linalg.det(mbeta1)
    beta2 = np.linalg.det(mbeta2)
    beta3 = -1*np.linalg.det(mbeta3)
    beta4 = np.linalg.det(mbeta4)
    gamma1 = np.linalg.det(mgamma1)
    gamma2 = -1*np.linalg.det(mgamma2)
    gamma3 = np.linalg.det(mgamma3)
    gamma4 = -1*np.linalg.det(mgamma4)
    delta1 = -1*np.linalg.det(mdelta1)
    delta2 = np.linalg.det(mdelta2)
    delta3 = -1*np.linalg.det(mdelta3)
    delta4 = np.linalg.det(mdelta4)
    B1 = np.array([[beta1,0,0],[0,gamma1,0],[0,0,delta1], \
         [gamma1,beta1,0],[0,delta1,gamma1],[delta1,0,beta1]],dtype=np.float64)
    B2 = np.array([[beta2,0,0],[0,gamma2,0],[0,0,delta2], \
         [gamma2,beta2,0],[0,delta2,gamma2],[delta2,0,beta2]],dtype=np.float64)
    B3 = np.array([[beta3,0,0],[0,gamma3,0],[0,0,delta3], \
         [gamma3,beta3,0],[0,delta3,gamma3],[delta3,0,beta3]],dtype=np.float64)
    B4 = np.array([[beta4,0,0],[0,gamma4,0],[0,0,delta4], \
         [gamma4,beta4,0],[0,delta4,gamma4],[delta4,0,beta4]],dtype=np.float64)
    aux = np.hstack((B1,B2,B3,B4))/(6*V)
    B = np.asmatrix(aux)
    y = V[0]*B.T*D*B
    return y, B

def TetrahedronAssemble(k,me,ne):
    row = np.reshape(np.hstack((np.repeat(3*np.reshape(me,(-1,1))-3,12,axis=1), \
                                np.repeat(3*np.reshape(me,(-1,1))-2,12,axis=1), \
                                np.repeat(3*np.reshape(me,(-1,1))-1,12,axis=1))),(1,-1))
    row = np.ravel(row)
    col = np.reshape(np.hstack((3*np.reshape(me,(-1,1))-3, \
                                3*np.reshape(me,(-1,1))-2, \
                                3*np.reshape(me,(-1,1))-1)),(ne,12))
    col = np.hstack((col,col,col,col,col,col,col,col,col,col,col,col))
    col = np.ravel(col)
    data = np.ravel(k)
    y = sp.csr_matrix((data, (row, col)))
    return y
