import numpy as np
import strassen

def conjugate_grad(A, b, x=None):
    n = b.shape[0]
    if not x:
        x = np.ones((n,1))
    r = np.dot(A, x) - b
    p = - r
    r_k_norm = strassen.FastMultiply(np.transpose(r), r)
    for _ in range(n*n):
        Ap = np.dot(A, p)
        alpha = r_k_norm / strassen.FastMultiply(np.transpose(p), Ap)
        x += alpha * p
        r += alpha * Ap
        r_kplus1_norm = strassen.FastMultiply(np.transpose(r), r)
        if np.sqrt(r_kplus1_norm) < 1e-6:
            break
        beta = r_kplus1_norm / r_k_norm
        r_k_norm = r_kplus1_norm
        p = beta * p - r
    return x

def DistNodes(f,s):
    if(s[0] == f[0]):
        aux = np.pi/2
        if(s[1] < f[1]):
             aux *= -1
        return (np.sqrt((s[0]-f[0])**2+(s[1]-f[1])**2),aux)
    else:
        aux = np.arctan((s[1]-f[1])/(s[0]-f[0]))
    if(aux < 0 and s[1] > f[1]):
        aux += np.pi
    if(s[1] < f[1]):
        aux += np.pi
        if(s[0] > f[0]):
            aux += np.pi
    return (np.sqrt((s[0]-f[0])**2+(s[1]-f[1])**2),aux)
 
 
def UBoundaryCondition(nU, u, i, NodesCondition):
    nU[i][0] = u
    NodesCondition.append(i)
 
def FBoundaryCondition(nF, f, i, ForcesCondition):
    nF[i][0] += f
    ForcesCondition.append(i)

def ElementStiffness(l, angle, A, I, E):
    w = [A*np.cos(angle)**2+12*I*np.sin(angle)**2/(l*l),
         A*np.sin(angle)**2+12*I*np.cos(angle)**2/(l*l),
         (A-12*I/(l*l))*np.cos(angle)*np.sin(angle),
         6*I*np.sin(angle)/l,
         6*I*np.cos(angle)/l]
    aux = [[w[0],w[2],-w[3],-w[0],-w[2],-w[3]],
           [w[2],w[1],w[4],-w[2],-w[1],w[4]],
           [-w[3],w[4],4*I,w[3],-w[4],2*I],
           [-w[0],-w[2],w[3],w[0],w[2],w[3]],
           [-w[2],-w[1],-w[4],w[2],w[1],-w[4]],
           [-w[3],w[4],2*I,w[3],-w[4],4*I]]
    aux = np.array(aux)
    aux = aux*E/l
    return aux

def AssemblyStiffness(nStiffnessMatrix, k, i, j):
    for p in range(0,3):
        for m in range(0,3):
            nStiffnessMatrix[3*i+p][3*i+m] += k[p][m]
            nStiffnessMatrix[3*i+p][3*j+m] += k[p][3+m]
            nStiffnessMatrix[3*j+p][3*i+m] += k[p+3][m]
            nStiffnessMatrix[3*j+p][3*j+m] += k[p+3][3+m]
 
def Initialize(nStiffnessMatrix, Elements, K, nU, nF, Nodes, NumberOfElement):
    for i in range(0,Nodes):
        nU[i][0] = 0
        nF[i][0] = 0
    for i in range(0,NumberOfElement):
        AssemblyStiffness(nStiffnessMatrix, K[i], int(Elements[i][0]), int(Elements[i][1]))

def PreSolvingStiffness(nStiffnessMatrix, Nodes, NodesCondition):
    nsize = Nodes-len(NodesCondition)
    newStiffness = np.zeros((nsize,nsize))
    contr = -1
    for i in range(0,Nodes):
        contc = -1
        flagr = False
        for k in range(0,len(NodesCondition)):
            if(i == NodesCondition[k]):
                flagr = True
                break
        if(flagr):
            continue
        contr += 1
        for j in range(0,Nodes):
            flagc = False
            for k in range(0,len(NodesCondition)):
                if(j == NodesCondition[k]):
                    flagc = True
                    break
            if(flagc):
                continue
            contc += 1
            newStiffness[contr][contc] = nStiffnessMatrix[i][j]
    return newStiffness
 
 
def PreSolvingF(nF, nS, nU, Nodes, NodesCondition):
    nsize = Nodes-len(NodesCondition)
    newF = np.zeros(nsize).reshape(nsize,1)
    contr = -1
    for i in range(0,Nodes):
        flagr = False
        for k in range(0,len(NodesCondition)):
            if(i == NodesCondition[k]):
                flagr = True
                break
        if(flagr):
            for k in range(0,Nodes):
                nF[k][0] = nF[k][0]-nS[k][i]*nU[i][0]
            continue
 
            
    for i in range(0,Nodes):
        flagr = False
        for k in range(0,len(NodesCondition)):
            if(i == NodesCondition[k]):
                flagr = True
                break
        if(flagr):
            continue
        contr += 1
        newF[contr][0] = nF[i][0]
    
    return newF
                      
 
def Solve(nStiffnessMatrix, nU, nF, StiffnessMatrix, Nodes, NodesCondition):
    newStiffness = PreSolvingStiffness(nStiffnessMatrix, Nodes, NodesCondition)
    newF = PreSolvingF(nF,nStiffnessMatrix,nU, Nodes, NodesCondition)
    u = conjugate_grad(newStiffness,newF)    
    #u = np.linalg.solve(newStiffness,newF)
    contr = -1
    for i in range(0,Nodes):
        flagr = False
        for k in range(0,len(NodesCondition)):
            if(i == NodesCondition[k]):
                flagr = True
                break
        if(flagr):
            continue
        contr += 1
        nU[i][0] = u[contr][0]
    nnF = strassen.FastMultiply(StiffnessMatrix,nU)
    return nU,nnF