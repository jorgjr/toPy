import numpy as np
from stl import mesh
import gmsh
import sys
import fea_functions
import opt_functions
import time


def TicTocGenerator():
    ti = 0
    tf = time.time()
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti

def toc(tempBool=True):
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    toc(False)


gmsh.initialize()

def setoptions():
    global boundnod, TicToc
    boundnod = False
    TicToc = TicTocGenerator()
    gmsh.option.setNumber("General.ShowModuleMenu", 1)
    gmsh.option.setNumber("General.MenuWidth", 340)
    gmsh.option.setNumber("General.ColorScheme", 0)
    gmsh.option.setNumber("General.FltkColorScheme", 1)
    gmsh.option.setNumber('Geometry.Tolerance', 1e-4)
    gmsh.option.setNumber('Geometry.Points', 1)
    gmsh.option.setNumber('Geometry.Lines', 1)
    gmsh.option.setNumber('Geometry.Surfaces', 1)
    gmsh.option.setNumber('Mesh.Points', 0)
    gmsh.option.setNumber('Mesh.PointType',1)
    gmsh.option.setNumber('Mesh.SurfaceEdges', 1)
    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
    gmsh.option.setNumber('Mesh.VolumeEdges', 1)
    gmsh.option.setNumber('Mesh.VolumeFaces', 1)
    gmsh.option.setNumber('Mesh.ColorCarousel', 0)
    gmsh.option.setColor('Mesh.Points', 0,0,210,255)
    gmsh.option.setColor('Mesh.Triangles', 210,210,0,255)
    gmsh.option.setColor('Mesh.Quadrangles', 210,210,0,255)
    gmsh.option.setColor('Mesh.Tetrahedra', 210,210,0,255)
    gmsh.option.setColor('Mesh.Hexahedra', 210,210,0,255)
    gmsh.option.setColor('Mesh.Prisms', 210,210,0,255)
    gmsh.option.setColor('Mesh.Pyramids', 210,210,0,255)
    gmsh.option.setColor('Mesh.Trihedra', 210,210,0,255)
setoptions()

parametermat = """
[
  { "type":"number", "name":"Material Data/Young Modulus [MPa]", "values":[1], "min":1e-2,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Material Data/Poisson's Ratio", "values":[0.3], "min":1e-2,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Material Data/[2D] Thickness [mm]", "values":[1], "min":1e-2,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"string", "name":"Material Data/[2D] Analysis Plane", "values":["Stress"],
   "choices":["Stress", "Strain"] },
  { "type":"string", "name":"FEA/Button", "values":["fea", ""],
    "visible":false }
]"""
gmsh.onelab.set(parametermat)

parameteropt = """
[
  { "type":"number", "name":"Optimization Data/Penalty", "values":[3.5], "min":1e-2,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Optimization Data/Volume Fraction [%]", "values":[50], "min":1,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Optimization Data/Minimum Radius [mm]", "values":[2.5], "min":1e-2,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Optimization Data/Convergency", "values":[0.1], "min":1e-2,
    "max":1e9, "step":1,"attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Optimization Data/Evolutionary Rate", "values":[2], "min":1,
    "max":1e9, "step":1,"attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Optimization Data/Max. Iterations", "values":[300], "min":11,
    "max":1e9, "step":1,"attributes":{"Highlight":"AliceBlue"} },
  { "type":"string", "name":"Optimization Data/Method", "values":["BESO-MP"],
   "choices":["BESO", "BESO-MP", "SIMP", "BESO-3D-Stress"] },
  { "type":"string", "name":"OPT/Button", "values":["opt", ""],
    "visible":false }
]"""
gmsh.onelab.set(parameteropt)

parameternod = """
[
  { "type":"number", "name":"Set Node/Tag", "values":[0], "min":0,
    "max":1e9, "step":1, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Set Node/X", "values":[0], "min":-1e9,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Set Node/Y", "values":[0], "min":-1e9,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Set Node/Z", "values":[0], "min":-1e9,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Set Node/Size at GeoPoint", "values":[0.05], "min":1e-4,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"string", "name":"ADD/Button", "values":["add", ""],
    "visible":false }
]"""
gmsh.onelab.set(parameternod)

parameterbc = """
[
  { "type":"number", "name":"Boundary at Node/Tag", "values":[1], "min":1,
    "max":1e9, "step":1, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"string", "name":"Boundary at Node/Type", "values":["Displacement"],
   "choices":["Displacement", "Force"] },
  { "type":"number", "name":"Boundary at Node/X", "values":[0], "min":-1e9,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Boundary at Node/Y", "values":[0], "min":-1e9,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"number", "name":"Boundary at Node/Z", "values":[0], "min":-1e9,
    "max":1e9, "attributes":{"Highlight":"AliceBlue"} },
  { "type":"string", "name":"BND/Button", "values":["bnd", ""],
    "visible":false }
]"""
gmsh.onelab.set(parameterbc)

parametermenu = """
[
  { "type":"string", "name":"SubMenu/Feature", "values":["F.E.A."],
   "choices":["F.E.A.", "Topology Optimization", "Reset Physical Groups", 
              "Set Node", "Set Size at GeoPoint", "Insert B.C. at Node",
              "Reset B.C. at Nodes"] },
  { "type":"string", "name":"ONELAB/Button", "values":["F.E.A.", "fea"],
    "visible":false }
]"""
gmsh.onelab.set(parametermenu)


def preset():
    global matnod, matel, nodtags, nodcoords, numnod, matnodpr, numel, cent, gl, nt, dim, aov, \
    ke, ket, ue, uet, K, Kt, F, Ft, U, Ut, B, D, nauxIx, nauxIy, nauxIz, fixel, freeel,maxit
    # fe, fet
    # Data about the model
    maxit = np.uint64(gmsh.onelab.getNumber("Optimization Data/Max. Iterations"))
    dim = gmsh.model.getDimension()
    if dim == 2:
        eltags, elnodtags = gmsh.model.mesh.getElementsByType(2,-1)
        gl = 2
        nt = 3
    elif dim == 3:
        eltags, elnodtags = gmsh.model.mesh.getElementsByType(4,-1)
        gl = 3
        nt = 4
    nodtags, nodcoords, _ = gmsh.model.mesh.getNodes(-1,-1)
    numnod = np.uint32(len(nodtags))
    matnod = np.zeros((numnod,3),dtype=np.float64)
    for i in range(0,numnod):
        matnod[i]=gmsh.model.mesh.getNode(i+1)[0]
    matnodpr = np.zeros((numnod,3),dtype=np.float64)
    print(numnod)
    numel = np.uint32(len(eltags))
    matel = np.reshape(np.asarray(elnodtags),(numel,nt))
    matel = np.uint32(matel)
    print(numel)
    # Pre-allocating memory
    aov = np.zeros((numel,1),dtype=np.float64)
    if dim == 2:
        B = np.zeros((numel,3,gl*nt),dtype=np.float64)
        D = np.zeros((3,3),dtype=np.float64)
        cent = gmsh.model.mesh.getBarycenters(2,-1, False, 'primary')
        for i in range(0,numel):
            v0 = matel[i][0]-1
            v1 = matel[i][1]-1
            v2 = matel[i][2]-1
            d01 = np.sqrt((matnod[v1][0]-matnod[v0][0])**2+(matnod[v1][1]-matnod[v0][1])**2+(matnod[v1][2]-matnod[v0][2])**2)
            d12 = np.sqrt((matnod[v2][0]-matnod[v1][0])**2+(matnod[v2][1]-matnod[v1][1])**2+(matnod[v2][2]-matnod[v1][2])**2)
            d20 = np.sqrt((matnod[v0][0]-matnod[v2][0])**2+(matnod[v0][1]-matnod[v2][1])**2+(matnod[v0][2]-matnod[v2][2])**2)
            s = (d01+d12+d20)/2
            aov[i] = np.sqrt((s*(s-d01)*(s-d12)*(s-d20)))
    elif dim == 3:
        B = np.zeros((numel,6,gl*nt),dtype=np.float64)
        D = np.zeros(1,dtype=np.float64)
        cent = gmsh.model.mesh.getBarycenters(4,-1, False, 'primary')
        for i in range(0,numel):           
            v0 = matel[i][0]-1
            v1 = matel[i][1]-1
            v2 = matel[i][2]-1
            v3 = matel[i][3]-1
            xyz = np.array([[1,1,1,1], \
                            [matnod[v0][0],matnod[v1][0],matnod[v2][0],matnod[v3][0]], \
                            [matnod[v0][1],matnod[v1][1],matnod[v2][1],matnod[v3][1]], \
                            [matnod[v0][2],matnod[v1][2],matnod[v2][2],matnod[v3][2]]])
            aov[i] = np.sqrt((np.linalg.det(xyz)/6)**2)
    cent = np.reshape(cent,(numel,3))
    if len(np.where(aov==0)[0])>0 :
        aov2 = aov
        for i in range(0,len(np.where(aov==0)[0])):
            matel = np.delete(matel,np.where(aov==0)[0][i]-i,0)
            aov2 = np.delete(aov2,np.where(aov==0)[0][i]-i,0)
            cent = np.delete(cent,np.where(aov==0)[0][i]-i,0)
            numel = numel-1
        aov = aov2
    if len(np.where(np.isnan(aov))[0])>0:
        aov2 = aov
        for i in range(0,len(np.where(np.isnan(aov))[0])):
            matel = np.delete(matel,np.where(np.isnan(aov))[0][i]-i,0)
            aov2 = np.delete(aov2,np.where(np.isnan(aov))[0][i]-i,0)
            cent = np.delete(cent,np.where(aov==0)[0][i]-i,0)
            numel = numel-1
        aov = aov2
    nauxIx = cent[:,0]
    nauxIy = cent[:,1]
    nauxIz = cent[:,2]
    ke = np.zeros((numel,gl*nt,gl*nt),dtype=np.float64)
    ket = np.zeros((numel,gl*nt,gl*nt),dtype=np.float64)
    ue = np.zeros((numel,gl*nt,1),dtype=np.float64)
    uet = np.zeros((numel,gl*nt,1),dtype=np.float64)
    # fe = np.zeros((numel,gl*nt,1),dtype=np.float64)
    # fet = np.zeros((numel,gl*nt,1),dtype=np.float64)
    K = np.zeros((gl*numnod,gl*numnod),dtype=np.float64)
    Kt = np.zeros((gl*numnod,gl*numnod),dtype=np.float64)
    U = np.zeros((gl*numnod,1),dtype=np.float64)  
    Ut = np.zeros((gl*numnod,1),dtype=np.float64)
    Ut[:] = np.nan
    F = np.zeros((gl*numnod,1),dtype=np.float64)   
    Ft = np.zeros((gl*numnod,1),dtype=np.float64)   
    # Boundary conditions
    if boundnod == True:
        for i in range(0,numbc):
            if bctype[i] == "Displacement":
                if dim == 2:
                    Ut[2*bctag[i]-2] = bcx
                    Ut[2*bctag[i]-1] = bcy
                elif dim == 3:
                    Ut[3*bctag[i]-3] = bcx
                    Ut[3*bctag[i]-2] = bcy
                    Ut[3*bctag[i]-1] = bcz
            elif bctype[i] == "Force":
                if dim == 2:
                    Ft[2*bctag[i]-2] = bcx
                    Ft[2*bctag[i]-1] = bcy
                elif dim == 3: 
                    Ft[3*bctag[i]-3] = bcx
                    Ft[3*bctag[i]-2] = bcy
                    Ft[3*bctag[i]-1] = bcz
    bounds = gmsh.model.getPhysicalGroups()
    fixel = np.array([],dtype=np.uint64)
    for b in bounds:
        bc = gmsh.model.getPhysicalName(b[0],b[1])
        print(bc)
        nbc = (bc.split("_"))
        if nbc[0] == "fixel":
            ntags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(b[0], b[1])
            for i in range(0,len(ntags)):
                auxfix = np.where(matel==ntags[i])[0]
                fixel = np.concatenate((fixel,auxfix),axis=None)
            # fixel = np.append(fixel,np.where(matel==ntags)[0])
            # fixel = np.uint64(np.unique(np.array(fixel)))
        else:
            if nbc[2] == "nan":
                x = np.nan
            elif nbc[2] == "fix":
                x = 0
            else:
                x = np.float64(nbc[2])
            if nbc[4] == "nan":
                y = np.nan
            elif nbc[4] == "fix":
                y = 0
            else:
                y = np.float64(nbc[4]) 
            if dim == 3:
                if nbc[6] == "nan":
                    z = np.nan
                elif nbc[6] == "fix":
                    z = 0
                else:
                    z = np.float64(nbc[6]) 
            if nbc[0] == "displacement":
                ntags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(b[0], b[1])
                ntags = np.uint32(ntags)
                if dim == 2:
                    for n in range(0,len(ntags)):
                        Ut[2*ntags[n]-2] = x
                        Ut[2*ntags[n]-1] = y
                elif dim == 3:
                    for n in range(0,len(ntags)):
                        Ut[3*ntags[n]-3] = x
                        Ut[3*ntags[n]-2] = y
                        Ut[3*ntags[n]-1] = z
            elif nbc[0] == "force":
                ntags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(b[0], b[1])
                ntags = np.uint32(ntags)
                xn = x/len(ntags)
                yn = y/len(ntags)
                if dim == 2:
                    for n in range(0,len(ntags)):
                        Ft[2*ntags[n]-2] = xn
                        Ft[2*ntags[n]-1] = yn
                elif dim == 3:
                    zn = z/len(ntags)
                    for n in range(0,len(ntags)):   
                        Ft[3*ntags[n]-3] = xn
                        Ft[3*ntags[n]-2] = yn
                        Ft[3*ntags[n]-1] = zn
    fixel = np.uint64(np.unique(fixel))
    freeel = np.arange(numel)
    jj = np.uint64(0)
    hum = np.uint64(1)
    for i in range(0,len(fixel)):
        freeel = np.delete(freeel,[fixel[i]-jj])
        jj = jj+hum
     
def insertbc():
    global boundnod, bctag, bctype, bcx, bcy, bcz, numbc
    if boundnod == False:
        boundnod = True
        bctag = []
        bctype = []
        bcx = []
        bcy = []
        bcz = []
        numbc = 0
    bctag.append(np.uint32(gmsh.onelab.getNumber("Boundary at Node/Tag"))[0])
    bctype.append(gmsh.onelab.getString("Boundary at Node/Type")[0])
    bcx.append(np.float64(gmsh.onelab.getNumber("Boundary at Node/X")))
    bcy.append(np.float64(gmsh.onelab.getNumber("Boundary at Node/Y")))
    bcz.append(np.float64(gmsh.onelab.getNumber("Boundary at Node/Z")))
    numbc = numbc+1
    print(numbc)

def feanalysis_2D():
    global E, NU, t, D, B, ke, K, F, U, ue, matnodpr
    preset()
    E = np.float64(gmsh.onelab.getNumber("Material Data/Young Modulus [MPa]"))
    NU = np.float64(gmsh.onelab.getNumber("Material Data/Poisson's Ratio"))
    t = np.float64(gmsh.onelab.getNumber("Material Data/[2D] Thickness [mm]"))
    if gmsh.onelab.getString("Material Data/[2D] Analysis Plane")[0] == "Stress":
        aux = np.asarray([[1,NU,0],[NU,1,0],[0,0,(1-NU)/2]],dtype=np.float64)
        D = (E/(1-NU*NU))*aux
        D = np.asmatrix(D)
    else:
        aux = np.asarray([[1-NU,NU,0],[NU,1-NU,0],[0,0,(1-2*NU)/2]],dtype=np.float64)
        D = (E/(1+NU)/(1-2*NU))*aux
        D = np.asmatrix(D)
    for i in range(0,numel):
        ke[i], B[i] = fea_functions.LinearTriangleElementStiffness(aov[i],D,t,\
                            matnod[matel[i][0]-1][0],matnod[matel[i][0]-1][1], \
                            matnod[matel[i][1]-1][0],matnod[matel[i][1]-1][1], \
                            matnod[matel[i][2]-1][0],matnod[matel[i][2]-1][1])
    K = fea_functions.LinearTriangleAssemble(ke,matel,numel)
    F, U = fea_functions.SolverFEM(Ut,K,Ft,U)
    # Elemental displacement
    # for i in range(0,numel):
    #     ue[i] = np.array(([[U[2*matel[i][0]-2][0]], [U[2*matel[i][0]-1][0]], \
    #                        [U[2*matel[i][1]-2][0]], [U[2*matel[i][1]-1][0]], \
    #                        [U[2*matel[i][2]-2][0]], [U[2*matel[i][2]-1][0]]]),dtype=np.float64)
    # Here goes stress and pstresses
    for i in range(0,numnod):
        matnodpr[i][0] = matnod[i][0]+U[2*i][0]
        matnodpr[i][1] = matnod[i][1]+U[2*i+1][0]
    for i in range(0,numnod):
        gmsh.model.mesh.setNode(i+1,(matnodpr[i][0],matnodpr[i][1],0),[])
    gmsh.model.mesh.generate()     
    gmsh.model.occ.synchronize()

def feanalysis_3D():
    global E, NU, t, D, B, ke, K, F, U, ue, matnodpr
    preset()
    E = np.float64(gmsh.onelab.getNumber("Material Data/Young Modulus [MPa]"))
    NU = np.float64(gmsh.onelab.getNumber("Material Data/Poisson's Ratio"))
    aux = np.asarray([[1-NU,NU,NU,0,0,0], \
                      [NU,1-NU,NU,0,0,0], \
                      [NU,NU,1-NU,0,0,0], \
                      [0,0,0,(1-2*NU)/2,0,0], \
                      [0,0,0,0,(1-2*NU)/2,0], \
                      [0,0,0,0,0,(1-2*NU)/2]],dtype=np.float64)
    D = (E/((1+NU)*(1-2*NU)))*aux
    D = np.asmatrix(D) 
    for i in range(0,numel):
        ke[i], B[i] = fea_functions.TetrahedronElementStiffness(aov[i],D,\
                            matnod[matel[i][0]-1][0],matnod[matel[i][0]-1][1],matnod[matel[i][0]-1][2], \
                            matnod[matel[i][1]-1][0],matnod[matel[i][1]-1][1],matnod[matel[i][1]-1][2], \
                            matnod[matel[i][2]-1][0],matnod[matel[i][2]-1][1],matnod[matel[i][2]-1][2], \
                            matnod[matel[i][3]-1][0],matnod[matel[i][3]-1][1],matnod[matel[i][3]-1][2])
    K = fea_functions.TetrahedronAssemble(ke,matel,numel)
    F, U = fea_functions.SolverFEM(Ut,K,Ft,U)
    # Elemental displacement
    # for i in range(0,numel):
    #     ue[i] = np.array(([[U[3*matel[i][0]-3][0]], [U[3*matel[i][0]-2][0]], [U[3*matel[i][0]-1][0]], \
    #                        [U[3*matel[i][1]-3][0]], [U[3*matel[i][1]-2][0]], [U[3*matel[i][1]-1][0]], \
    #                        [U[3*matel[i][2]-3][0]], [U[3*matel[i][2]-2][0]], [U[3*matel[i][2]-1][0]], \
    #                        [U[3*matel[i][3]-3][0]], [U[3*matel[i][3]-2][0]], [U[3*matel[i][3]-1][0]]]),dtype=np.float64)
    # Here goes stress and pstresses
    for i in range(0,numnod):
        matnodpr[i][0] = matnod[i][0]+U[3*i][0]
        matnodpr[i][1] = matnod[i][1]+U[3*i+1][0]
        matnodpr[i][2] = matnod[i][2]+U[3*i+2][0]
    for i in range(0,numnod):
        gmsh.model.mesh.setNode(i+1,(matnodpr[i][0],matnodpr[i][1],matnodpr[i][2]),[])
    gmsh.model.mesh.generate()     
    gmsh.model.occ.synchronize()

def topy_opt_2D():
    global E, NU, t, D, ke, ket, uet, B, K, Kt, F, U, x, alpha
    preset()
    E = np.float64(gmsh.onelab.getNumber("Material Data/Young Modulus [MPa]"))
    NU = np.float64(gmsh.onelab.getNumber("Material Data/Poisson's Ratio"))
    t = np.float64(gmsh.onelab.getNumber("Material Data/[2D] Thickness [mm]"))
    pe = np.float64(gmsh.onelab.getNumber("Optimization Data/Penalty"))
    rmin = np.float64(gmsh.onelab.getNumber("Optimization Data/Minimum Radius [mm]"))
    conv = np.float64(gmsh.onelab.getNumber("Optimization Data/Convergency"))/100
    er = np.float64(gmsh.onelab.getNumber("Optimization Data/Evolutionary Rate"))/100
    vfrac = np.float64(gmsh.onelab.getNumber("Optimization Data/Volume Fraction [%]"))/100
    idxf = np.where(Ft!=0)[0]
    Ff = np.zeros((len(idxf),gl*numnod,1),dtype=np.float64) 
    for i in range(0,len(idxf)):
        Ff[i][idxf[i]] = Ft[idxf[i]]
    x = np.zeros((numel,1),dtype=np.float64) 
    alpha = np.zeros((numel,1),dtype=np.float64)
    alphan = np.zeros((numel,1),dtype=np.float64) 
    oldalpha = np.zeros((numel,1),dtype=np.float64) 
    ce = np.zeros((numel,1),dtype=np.float64) 
    kk = np.arange(numel)
    vifac = t*aov
    x[:] = 1
    cc = []
    vol = np.nanmax((1*(1-er),vfrac))
    xmin = np.float64(0.001)
    change = 1
    ii = 0
    if gmsh.onelab.getString("Material Data/[2D] Analysis Plane")[0] == "Stress":
        aux = np.asarray([[1,NU,0],[NU,1,0],[0,0,(1-NU)/2]],dtype=np.float64)
        D = (E/(1-NU*NU))*aux
        D = np.asmatrix(D)
    else:
        aux = np.asarray([[1-NU,NU,0],[NU,1-NU,0],[0,0,(1-2*NU)/2]],dtype=np.float64)
        D = (E/(1+NU)/(1-2*NU))*aux
        D = np.asmatrix(D)
    for i in range(0,numel):
        ke[i], B[i] = fea_functions.LinearTriangleElementStiffness(aov[i],D,t, \
                            matnod[matel[i][0]-1][0],matnod[matel[i][0]-1][1], \
                            matnod[matel[i][1]-1][0],matnod[matel[i][1]-1][1], \
                            matnod[matel[i][2]-1][0],matnod[matel[i][2]-1][1])
    H1, H2 = opt_functions.allocH1H2(numel,nauxIx,nauxIy,nauxIz,rmin,vifac,dim,t,cent)
    if gmsh.onelab.getString("Optimization Data/Method")[0] == "SIMP":
       xold = np.zeros((numel,1),dtype=np.float64) 
       x[:] = vfrac
       xnew = np.zeros((numel,1),dtype=np.float64)
       while change > conv:
            tic()
            cc.append(0)
            alpha.fill(0)
            np.copyto(xold,x)
            # F.E.A.
            for i in range(0,numel):
                ket[i]=ke[i]*(x[i]**pe)
            Kt = fea_functions.LinearTriangleAssemble(ket,matel,numel)
            for j in range(0,len(idxf)):
                U = fea_functions.SolverFEM_opt(Ut,Kt,Ff[j],U)
                for i in range(0,numel):
                    uet[i] = np.array(([[U[2*matel[i][0]-2][0]], [U[2*matel[i][0]-1][0]], \
                                        [U[2*matel[i][1]-2][0]], [U[2*matel[i][1]-1][0]], \
                                        [U[2*matel[i][2]-2][0]], [U[2*matel[i][2]-1][0]]]),dtype=np.float64)
                # Objective and Sensitivity
                for i in range(0,numel):
                    ce[i] = np.asmatrix(uet[i]).T*np.dot(ket[i],uet[i])
                    cc[ii] = cc[ii]+np.asarray((x[i][0]**pe)*ce[i])[0]
                    alpha[i] = alpha[i]+np.asarray(((-pe*x[i][0]**(pe-1))*ce[i])/vifac[i][0])[0]
            # Filtering Sensitivity
            xp = alpha*x
            for i in range(0,numel):
                alphan[i] = np.dot(H1[i,kk],xp)/(H2[i,i]*x[i])
            # Design Update and Plot
            l1 = 0
            l2 = 1e9
            move = 0.2
            n = 0.5
            xnew[fixel] = 1
            while (l2-l1 > 1.0e-4):
               th = 0.5*(l2+l1)
               xnew[freeel] = np.maximum(0.001,np.maximum(x[freeel]-move,np.minimum(1,np.minimum(x[freeel]+move,x*(-alphan[freeel]/th)**n))))
               if np.sum(np.sum(xnew))-np.sum(np.sum(numel*vfrac)) > 0:
                  l1 = th
               else:
                  l2 = th
            np.copyto(x,xnew)
            change = np.nanmax(np.absolute(x-xold))
            print(change)
            print(cc[ii])
            ii = ii+1
            print(ii)
            if ii > maxit:
                break
            toc()
    else:     
        while change > conv:
            tic()
            cc.append(0)
            alpha.fill(0)
            vol = np.maximum(vol*(1-er),vfrac)
            if ii > 0:
                np.copyto(oldalpha,alphan)
            # F.E.A.
            for i in range(0,numel):
                ket[i]=ke[i]*(x[i]**pe)
            Kt = fea_functions.LinearTriangleAssemble(ket,matel,numel)
            for j in range(0,len(idxf)):
                U = fea_functions.SolverFEM_opt(Ut,Kt,Ff[j],U)
                for i in range(0,numel):
                    uet[i] = np.array(([[U[2*matel[i][0]-2][0]], [U[2*matel[i][0]-1][0]], \
                                        [U[2*matel[i][1]-2][0]], [U[2*matel[i][1]-1][0]], \
                                        [U[2*matel[i][2]-2][0]], [U[2*matel[i][2]-1][0]]]),dtype=np.float64)
                # Objective and Sensitivity
                for i in range(0,numel):
                    ce[i] = np.asmatrix(uet[i]).T*np.dot(ket[i],uet[i])
                    cc[ii] = cc[ii]+np.asarray((0.5*x[i][0]**pe)*ce[i])[0]
                    alpha[i] = alpha[i]+np.asarray(((0.5*x[i][0]**(pe-1))*ce[i])/vifac[i][0])[0]
            # Filtering Sensitivity
            for i in range(0,numel):
                alphan[i] = np.dot(H1[i,kk],alpha)/H2[i,i]
            # Stabilization  
            if ii>0:
                alphan = (oldalpha+alphan)/2    
            # Design Update and Plot
            l1 = np.nanmin(np.nanmin(alphan))
            l2 = np.nanmax(np.nanmax(alphan))
            x[fixel] = 1
            while ((l2-l1)/l2>1.0e-5):
               th = (l1+l2)/2.0
               x[freeel] = np.maximum(xmin,np.sign(alphan[freeel]-th))
               if np.sum(np.sum(aov*x))-np.sum(np.sum(aov*vol))>0:
                  l1 = th
               else:
                  l2 = th
            if ii>9:
                change = np.abs(np.sum(cc[ii-9:ii-4])-np.sum(cc[ii-4:ii+1]))/np.sum(cc[ii-4:ii+1])
                print(change)
                print(cc[ii])
            ii = ii+1
            print(ii)
            if ii > maxit:
                break
            toc()

def topy_opt_3D():
    global E, NU, D, ke, ket, uet, B, K, Kt, F, U, x, alpha, poc
    preset()
    E = np.float64(gmsh.onelab.getNumber("Material Data/Young Modulus [MPa]"))
    NU = np.float64(gmsh.onelab.getNumber("Material Data/Poisson's Ratio"))
    pe = np.float64(gmsh.onelab.getNumber("Optimization Data/Penalty"))
    rmin = np.float64(gmsh.onelab.getNumber("Optimization Data/Minimum Radius [mm]"))
    conv = np.float64(gmsh.onelab.getNumber("Optimization Data/Convergency"))/100
    er = np.float64(gmsh.onelab.getNumber("Optimization Data/Evolutionary Rate"))/100
    vfrac = np.float64(gmsh.onelab.getNumber("Optimization Data/Volume Fraction [%]"))/100
    idxf = np.where(Ft!=0)[0]
    Ff = np.zeros((len(idxf),gl*numnod,1),dtype=np.float64) 
    for i in range(0,len(idxf)):
        Ff[i][idxf[i]] = Ft[idxf[i]]
    x = np.zeros((numel,1),dtype=np.float64)
    alpha = np.zeros((numel,1),dtype=np.float64)
    alphan = np.zeros((numel,1),dtype=np.float64) 
    oldalpha = np.zeros((numel,1),dtype=np.float64)
    ce = np.zeros((numel,1),dtype=np.float64)
    kk = np.arange(numel)
    x[:] = 1
    cc = []
    vol = np.nanmax((1*(1-er),vfrac))
    xmin = np.float64(0.001)
    change = 1
    ii = 0
    aux = np.asarray([[1-NU,NU,NU,0,0,0], \
                      [NU,1-NU,NU,0,0,0], \
                      [NU,NU,1-NU,0,0,0], \
                      [0,0,0,(1-2*NU)/2,0,0], \
                      [0,0,0,0,(1-2*NU)/2,0], \
                      [0,0,0,0,0,(1-2*NU)/2]],dtype=np.float64)
    D = (E/((1+NU)*(1-2*NU)))*aux
    D = np.asmatrix(D) 
    for i in range(0,numel):
        ke[i], B[i] = fea_functions.TetrahedronElementStiffness(aov[i],D, \
                            matnod[matel[i][0]-1][0],matnod[matel[i][0]-1][1],matnod[matel[i][0]-1][2], \
                            matnod[matel[i][1]-1][0],matnod[matel[i][1]-1][1],matnod[matel[i][1]-1][2], \
                            matnod[matel[i][2]-1][0],matnod[matel[i][2]-1][1],matnod[matel[i][2]-1][2], \
                            matnod[matel[i][3]-1][0],matnod[matel[i][3]-1][1],matnod[matel[i][3]-1][2])
    if gmsh.onelab.getString("Optimization Data/Method")[0] == "SIMP":
       H1, H2 = opt_functions.allocH1H2(numel,nauxIx,nauxIy,nauxIz,rmin,aov,dim,1,cent)
       xold = np.zeros((numel,1),dtype=np.float64) 
       x[:] = vfrac
       xnew = np.zeros((numel,1),dtype=np.float64)
       while change > conv:
            tic()
            cc.append(0)
            alpha.fill(0)
            np.copyto(xold,x)
            # F.E.A.
            for i in range(0,numel):
                ket[i]=ke[i]*(x[i]**pe)
            Kt = fea_functions.TetrahedronAssemble(ket,matel,numel)
            for j in range(0,len(idxf)):
                U = fea_functions.SolverFEM_opt(Ut,Kt,Ff[j],U)
                for i in range(0,numel):
                    uet[i] = np.array(([[U[3*matel[i][0]-3][0]], [U[3*matel[i][0]-2][0]], [U[3*matel[i][0]-1][0]], \
                                        [U[3*matel[i][1]-3][0]], [U[3*matel[i][1]-2][0]], [U[3*matel[i][1]-1][0]], \
                                        [U[3*matel[i][2]-3][0]], [U[3*matel[i][2]-2][0]], [U[3*matel[i][2]-1][0]], \
                                        [U[3*matel[i][3]-3][0]], [U[3*matel[i][3]-2][0]], [U[3*matel[i][3]-1][0]]]),dtype=np.float64)
                # Objective and Sensitivity
                for i in range(0,numel):
                    ce[i] = np.asmatrix(uet[i]).T*np.dot(ket[i],uet[i])
                    cc[ii] = cc[ii]+np.asarray((x[i][0]**pe)*ce[i])[0]
                    alpha[i] = alpha[i]+np.asarray(((-pe*x[i][0]**(pe-1))*ce[i])/aov[i][0])[0]
            # Filtering Sensitivity
            xp = alpha*x
            for i in range(0,numel):
                alphan[i] = np.dot(H1[i,kk],xp)/(H2[i,i]*x[i])
            # Design Update and Plot
            l1 = 0
            l2 = 1e9
            move = 0.2
            n = 0.5
            xnew[fixel] = 1
            while (l2-l1 > 1.0e-4):
               th = 0.5*(l2+l1)
               xnew[freeel] = np.maximum(0.001,np.maximum(x[freeel]-move,np.minimum(1,np.minimum(x[freeel]+move,x*(-alphan[freeel]/th)**n))))
               if np.sum(np.sum(xnew))-np.sum(np.sum(numel*vfrac)) > 0:
                  l1 = th
               else:
                  l2 = th
            np.copyto(x,xnew)
            change = np.nanmax(np.absolute(x-xold))
            print(change)
            print(cc[ii])
            ii = ii+1
            print(ii)
            if ii > maxit:
                break
            toc()
    else:     
        if gmsh.onelab.getString("Optimization Data/Method")[0] == "BESO":
            nod_inrad = opt_functions.rmin_elnodes(numel,rmin,cent,matel,matnod,nt)
            H1, H2, H3 = opt_functions.allocH1H2H3(matnod,numel,numnod,nauxIx,nauxIy,nauxIz,rmin,aov,dim,nt,cent,nod_inrad)
            alphann = np.zeros((numnod,1),dtype=np.float64)
            kkk = np.arange(numnod)
        else:
            H1, H2 = opt_functions.allocH1H2(numel,nauxIx,nauxIy,nauxIz,rmin,aov,dim,1,cent)
        while change > conv:
            tic()
            cc.append(0)
            alpha.fill(0)
            vol = np.maximum(vol*(1-er),vfrac)
            if ii > 0:
                np.copyto(oldalpha,alphan)
            # F.E.A.
            for i in range(0,numel):
                ket[i]=ke[i]*(x[i]**pe)
            Kt = fea_functions.TetrahedronAssemble(ket,matel,numel)
            for j in range(0,len(idxf)):
                U = fea_functions.SolverFEM_opt(Ut,Kt,Ff[j],U)
                # Elemental displacement
                for i in range(0,numel):
                    uet[i] = np.array(([[U[3*matel[i][0]-3][0]], [U[3*matel[i][0]-2][0]], [U[3*matel[i][0]-1][0]], \
                                        [U[3*matel[i][1]-3][0]], [U[3*matel[i][1]-2][0]], [U[3*matel[i][1]-1][0]], \
                                        [U[3*matel[i][2]-3][0]], [U[3*matel[i][2]-2][0]], [U[3*matel[i][2]-1][0]], \
                                        [U[3*matel[i][3]-3][0]], [U[3*matel[i][3]-2][0]], [U[3*matel[i][3]-1][0]]]),dtype=np.float64)
                # Objective and Sensitivity
                if gmsh.onelab.getString("Optimization Data/Method")[0] == "BESO":
                    for i in range(0,numel):
                        ce[i] = np.asmatrix(uet[i]).T*np.dot(ket[i],uet[i])
                        cc[ii] = cc[ii]+np.asarray(0.5*ce[i])[0]
                        alpha[i] = alpha[i]+np.asarray((0.5*ce[i])/aov[i][0])[0]
                else:
                    for i in range(0,numel):
                        ce[i] = np.asmatrix(uet[i]).T*np.dot(ket[i],uet[i])
                        cc[ii] = cc[ii]+np.asarray((0.5*x[i][0]**pe)*ce[i])[0]
                        alpha[i] = alpha[i]+np.asarray(((0.5*x[i][0]**(pe-1))*ce[i])/aov[i][0])[0]
            # Filtering Sensitivity
            if gmsh.onelab.getString("Optimization Data/Method")[0] == "BESO":
                for i in range(0,numnod):
                    alphann[i] = np.dot(H1[i,kk],alpha)
                for i in range(0,numel):
                    alphan[i] = np.dot(H2[i,kkk],alphann)/H3[i,i]
            else:
                for i in range(0,numel):
                    alphan[i] = np.dot(H1[i,kk],alpha)/H2[i,i]
            # Stabilization
            if ii>0:
                alphan = (oldalpha+alphan)/2
            # Design Update and Plot
            l1 = np.nanmin(np.nanmin(alphan))
            l2 = np.nanmax(np.nanmax(alphan))
            x[fixel] = 1
            while ((l2-l1)/l2>1.0e-5):
               th = (l1+l2)/2.0
               x[freeel] = np.maximum(xmin,np.sign(alphan[freeel]-th))
               if np.sum(np.sum(aov*x))-np.sum(np.sum(aov*vol))>0:
                  l1 = th
               else:
                  l2 = th
            if ii>9:
                change = np.abs(np.sum(cc[ii-9:ii-4])-np.sum(cc[ii-4:ii+1]))/np.sum(cc[ii-4:ii+1])
                print(change)
                print(cc[ii])
            ii = ii+1
            print(ii)
            if ii > maxit:
                break
            toc()
            
def topy_opt_3D_stress():
    global E, NU, D, ke, ket, uet, B, K, Kt, F, U, x, alpha, poc
    preset()
    E = np.float64(gmsh.onelab.getNumber("Material Data/Young Modulus [MPa]"))
    NU = np.float64(gmsh.onelab.getNumber("Material Data/Poisson's Ratio"))
    pe = np.float64(gmsh.onelab.getNumber("Optimization Data/Penalty"))
    rmin = np.float64(gmsh.onelab.getNumber("Optimization Data/Minimum Radius [mm]"))
    conv = np.float64(gmsh.onelab.getNumber("Optimization Data/Convergency"))/100
    er = np.float64(gmsh.onelab.getNumber("Optimization Data/Evolutionary Rate"))/100
    vfrac = np.float64(gmsh.onelab.getNumber("Optimization Data/Volume Fraction [%]"))/100
    idxf = np.where(Ft!=0)[0]
    Ff = np.zeros((len(idxf),gl*numnod,1),dtype=np.float64) 
    for i in range(0,len(idxf)):
        Ff[i][idxf[i]] = Ft[idxf[i]]
    stresses = np.zeros((numel,6,1),dtype=np.float64)
    pstresses = np.zeros((numel,3,1),dtype=np.float64)
    vmstresses = np.zeros((numel,1),dtype=np.float64) 
    x = np.zeros((numel,1),dtype=np.float64) 
    alpha = np.zeros((numel,1),dtype=np.float64) 
    alphan = np.zeros((numel,1),dtype=np.float64) 
    oldalpha = np.zeros((numel,1),dtype=np.float64) 
    ce = np.zeros((numel,1),dtype=np.float64) 
    kk = np.arange(numel)
    x[:] = 1
    cc = []
    vol = np.nanmax((1*(1-er),vfrac))
    xmin = np.float64(0.001)
    change = 1
    ii = 0
    aux = np.asarray([[1-NU,NU,NU,0,0,0], \
                      [NU,1-NU,NU,0,0,0], \
                      [NU,NU,1-NU,0,0,0], \
                      [0,0,0,(1-2*NU)/2,0,0], \
                      [0,0,0,0,(1-2*NU)/2,0], \
                      [0,0,0,0,0,(1-2*NU)/2]],dtype=np.float64)
    D = (E/((1+NU)*(1-2*NU)))*aux
    D = np.asmatrix(D) 
    for i in range(0,numel):
        ke[i], B[i] = fea_functions.TetrahedronElementStiffness(aov[i],D, \
                            matnod[matel[i][0]-1][0],matnod[matel[i][0]-1][1],matnod[matel[i][0]-1][2], \
                            matnod[matel[i][1]-1][0],matnod[matel[i][1]-1][1],matnod[matel[i][1]-1][2], \
                            matnod[matel[i][2]-1][0],matnod[matel[i][2]-1][1],matnod[matel[i][2]-1][2], \
                            matnod[matel[i][3]-1][0],matnod[matel[i][3]-1][1],matnod[matel[i][3]-1][2])
    # H1, H2 = opt_functions.allocH1H2(numel,nauxIx,nauxIy,nauxIz,rmin,aov,dim,1,cent)
    nod_inrad = opt_functions.rmin_elnodes(numel,rmin,cent,matel,matnod,nt)
    H1, H2, H3 = opt_functions.allocH1H2H3(matnod,numel,numnod,nauxIx,nauxIy,nauxIz,rmin,aov,dim,nt,cent,nod_inrad)
    alphann = np.zeros((numnod,1),dtype=np.float64)
    kkk = np.arange(numnod)
    while change > conv:
        tic()
        cc.append(0)
        alpha.fill(0)
        vol = np.maximum(vol*(1-er),vfrac)
        if ii > 0:
            np.copyto(oldalpha,alphan)
        # F.E.A.
        for i in range(0,numel):
            ket[i]=ke[i]*(x[i]**pe)
        Kt = fea_functions.TetrahedronAssemble(ket,matel,numel)
        # Elemental displacement
        for j in range(0,len(idxf)):
            U = fea_functions.SolverFEM_opt(Ut,Kt,Ff[j],U)
            # Elemental displacement
            for i in range(0,numel):
                uet[i] = np.array(([[U[3*matel[i][0]-3][0]], [U[3*matel[i][0]-2][0]], [U[3*matel[i][0]-1][0]], \
                                    [U[3*matel[i][1]-3][0]], [U[3*matel[i][1]-2][0]], [U[3*matel[i][1]-1][0]], \
                                    [U[3*matel[i][2]-3][0]], [U[3*matel[i][2]-2][0]], [U[3*matel[i][2]-1][0]], \
                                    [U[3*matel[i][3]-3][0]], [U[3*matel[i][3]-2][0]], [U[3*matel[i][3]-1][0]]]),dtype=np.float64)
                # stresses[i] = x[i][0]*D*B[i]*uet[i]
                stresses[i] = D*B[i]*uet[i]
                pstresses[i] = fea_functions.TetrahedronElementPStresses(stresses[i])
                vmstresses[i] = (0.5*((pstresses[i][0]-pstresses[i][1])**2+ \
                                      (pstresses[i][1]-pstresses[i][2])**2+ \
                                      (pstresses[i][2]-pstresses[i][0])**2))
            # Objective and Sensitivity
            for i in range(0,numel):
                ce[i] = np.asmatrix(uet[i]).T*np.dot(ket[i],uet[i])
                cc[ii] = cc[ii]+np.asarray((0.5*x[i][0]**pe)*ce[i])[0]
                alpha[i] = alpha[i]+np.asarray(((x[i][0]**(pe-1))*vmstresses[i]))[0]
        # Filtering Sensitivity
        # for i in range(0,numel):
        #     alphan[i] = np.dot(H1[i,kk],alpha)/H2[i,i]
        for i in range(0,numnod):
            alphann[i] = np.dot(H1[i,kk],alpha)
        for i in range(0,numel):
            alphan[i] = np.dot(H2[i,kkk],alphann)/H3[i,i]
        # Stabilization
        if ii>0:
            alphan = (oldalpha+alphan)/2
        # Design Update and Plot
        l1 = np.nanmin(np.nanmin(alphan))
        l2 = np.nanmax(np.nanmax(alphan))
        x[fixel] = 1
        while ((l2-l1)/l2>1.0e-5):
           th = (l1+l2)/2.0
           x[freeel] = np.maximum(xmin,np.sign(alphan[freeel]-th))
           if np.sum(np.sum(aov*x))-np.sum(np.sum(aov*vol))>0:
              l1 = th
           else:
              l2 = th
        if ii>9:
            change = np.abs(np.sum(cc[ii-9:ii-4])-np.sum(cc[ii-4:ii+1]))/np.sum(cc[ii-4:ii+1])
            print(change)
            print(cc[ii])
        ii = ii+1
        print(ii)
        if ii > maxit:
            break
        toc()
        
def post_topy():
    nodeTags = {}
    nodeCoords = {}
    elementTypes = {}
    elementTags = {}
    elementNodeTags = {}
    entities = gmsh.model.getEntities()
    for e in entities:
        nodeTags[e], nodeCoords[e], _ = gmsh.model.mesh.getNodes(e[0], e[1])
        elementTypes[e], elementTags[e], elementNodeTags[e] = gmsh.model.mesh.getElements(e[0], e[1])
    nmatel = np.reshape(elementNodeTags[dim,1],(numel,nt))[np.where(x==1)[0]]
    if dim == 2:
        geom = mesh.Mesh(np.zeros(nmatel.shape[0], dtype=mesh.Mesh.dtype))
        for i, f in enumerate((nmatel-1)):
            for j in range(3):
                geom.vectors[i][j] = matnod[f[j],:]
    else:
        geom = mesh.Mesh(np.zeros(nmatel.shape[0]*4, dtype=mesh.Mesh.dtype))
        for i, f in enumerate((nmatel-1)):
            geom.vectors[4*i][0] = matnod[f[0],:]
            geom.vectors[4*i][1] = matnod[f[1],:]
            geom.vectors[4*i][2] = matnod[f[2],:]
            geom.vectors[4*i+1][0] = matnod[f[0],:]
            geom.vectors[4*i+1][1] = matnod[f[1],:]
            geom.vectors[4*i+1][2] = matnod[f[3],:]
            geom.vectors[4*i+2][0] = matnod[f[0],:]
            geom.vectors[4*i+2][1] = matnod[f[2],:]
            geom.vectors[4*i+2][2] = matnod[f[3],:]
            geom.vectors[4*i+3][0] = matnod[f[1],:]
            geom.vectors[4*i+3][1] = matnod[f[2],:]
            geom.vectors[4*i+3][2] = matnod[f[3],:]
    geom.save('opt_geom.stl')
    nmatel = np.reshape(nmatel,-1)    
    nmin = np.uint64(len(nmatel)/nt)
    neltags = elementTags[dim,1]
    elementTags[dim,1] = [neltags[0][:nmin]]
    elementNodeTags[dim,1] = [nmatel]
    gmsh.model.mesh.clear()
    for e in entities:
        gmsh.model.mesh.addNodes(e[0], e[1], nodeTags[e], nodeCoords[e])
        if e[0] == dim:
            gmsh.model.mesh.addElements(e[0], e[1], elementTypes[e], elementTags[e], elementNodeTags[e])
    if dim == 2:
        gmsh.option.setNumber('Mesh.SurfaceEdges', 0)
        gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
        gmsh.option.setNumber('Mesh.VolumeEdges', 0)
        gmsh.option.setNumber('Mesh.VolumeFaces', 0)
    else:
        gmsh.option.setNumber('Mesh.SurfaceEdges', 0)
        gmsh.option.setNumber('Mesh.SurfaceFaces', 0)
        gmsh.option.setNumber('Mesh.VolumeEdges', 0)
        gmsh.option.setNumber('Mesh.VolumeFaces', 1)
    

if '-nopopup' not in sys.argv:
    gmsh.fltk.initialize()
gmsh.fltk.update()

# GUI
while 1:
    if gmsh.fltk.isAvailable() == 0:
        break
    gmsh.fltk.wait()
    while 1:
        btn = gmsh.onelab.getString("SubMenu/Feature")
        if btn[0] == gmsh.onelab.getString("ONELAB/Button")[0]:
            break
        else:
            if btn[0] == "F.E.A.":
                bt = "fea"
            elif btn[0] == "Topology Optimization":
                bt = "topy"
            elif btn[0] == "Set Node":
                bt = "set"
            elif btn[0] == "Set Size at GeoPoint":
                bt = "sizing"
            elif btn[0] == "Insert B.C. at Node":
                bt = "insertbc"
            elif btn[0] == "Reset B.C. at Nodes":
                bt = "resetbc"
            elif btn[0] == "Reset Physical Groups":
                bt = "resetpg"
            gmsh.onelab.setString("ONELAB/Action", [""])
            gmsh.onelab.setString("ONELAB/Button", [btn[0], bt])
            gmsh.fltk.update()
            break
    action = gmsh.onelab.getString("ONELAB/Action")
    if len(action) < 1:
        continue
    elif action[0] == "fea":
        if gmsh.model.getDimension() == 2:
            feanalysis_2D()
        elif gmsh.model.getDimension() == 3:
            feanalysis_3D()
        gmsh.onelab.setString("ONELAB/Action", [""])
        gmsh.onelab.setString("ONELAB/Button", ["F.E.A.", "fea"])
        gmsh.fltk.update()
    elif action[0] == "topy":
        if gmsh.model.getDimension() == 2:
            topy_opt_2D()
        elif gmsh.model.getDimension() == 3:
            if gmsh.onelab.getString("Optimization Data/Method")[0] == "BESO-3D-Stress": 
                topy_opt_3D_stress()
            else:
                topy_opt_3D()
        post_topy()
        gmsh.onelab.setString("ONELAB/Action", [""])
        gmsh.onelab.setString("ONELAB/Button", ["Topology Optimization", "topy"])
        gmsh.fltk.update()
    elif action[0] == "set":
        tc, _ = gmsh.model.mesh.getNode(np.uint32(gmsh.onelab.getNumber("Set Node/Tag"))[0])
        print(tc)
        gmsh.model.mesh.setNode(np.uint32(gmsh.onelab.getNumber("Set Node/Tag"))[0],
                                (gmsh.onelab.getNumber("Set Node/X"),
                                 gmsh.onelab.getNumber("Set Node/Y"),
                                 gmsh.onelab.getNumber("Set Node/Z")),[])
        gmsh.onelab.setString("ONELAB/Action", [""])
        gmsh.onelab.setString("ONELAB/Button", ["Set Node", "set"])
        tc, _ = gmsh.model.mesh.getNode(np.uint32(gmsh.onelab.getNumber("Set Node/Tag"))[0])
        gmsh.model.mesh.generate()
        print(tc)
        gmsh.fltk.update()
    elif action[0] == "sizing":
        gmsh.model.geo.mesh.setSize((0,np.uint32(gmsh.onelab.getNumber("Set Node/Tag"))[0]),
                                gmsh.onelab.getNumber("Set Node/Size at GeoPoint"))
        gmsh.onelab.setString("ONELAB/Action", [""])
        gmsh.onelab.setString("ONELAB/Button", ["Set Size at GeoPoint", "sizing"])
        gmsh.model.mesh.generate()
        gmsh.fltk.update()
    elif action[0] == "insertbc":
        gmsh.onelab.setString("ONELAB/Action", [""])
        gmsh.onelab.setString("ONELAB/Button", ["Insert B.C. at Node", "insertbc"])
        insertbc()
        gmsh.fltk.update()
    elif action[0] == "resetbc":
        global boundnod
        boundnod = False
        gmsh.onelab.setString("ONELAB/Action", [""])
        gmsh.onelab.setString("ONELAB/Button", ["Reset B.C. at Nodes", "resetbc"])
    elif action[0] == "resetpg":
        bounds = gmsh.model.getPhysicalGroups()
        for b in bounds:
            gmsh.model.removePhysicalName(gmsh.model.getPhysicalName(b[0],b[1]))
        gmsh.model.removePhysicalGroups()
        gmsh.onelab.setString("ONELAB/Action", [""])
        gmsh.onelab.setString("ONELAB/Button", ["Reset Physical Groups", "resetpg"])
        gmsh.fltk.update()
    
gmsh.finalize()