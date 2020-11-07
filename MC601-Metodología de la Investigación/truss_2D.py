import numpy as np
import matplotlib.pyplot as plt
import strassen
import FEA_functions as fea

def new_model(PosNodes, Elements, NodesConditions, ForcesConditions):
    x = set()
    
    for i in range(len(Elements)):
        x.add(Elements[i][0])
        x.add(Elements[i][1])
    
    sz = len(set)

    nPosNodes = [0]*sz
    nElements = [0]*len(Elements)
    nNodesConditions = [0]*len(NodesConditions)
    nForcesConditions = [0]*len(ForcesConditions)
    
    sz = 0
    map_nodes = {}
    for i in x:
        map_nodes[i] = sz
        nPosNodes[sz] = PosNodes[i]
        sz += 1

    for i in range(0, len(Elements)):
        nElements[i] = (map_nodes[Elements[i][0]], map_nodes[Elements[i][1]])
    
    for i in range(0, len(NodesConditions)):
        nNodesConditions = (map_nodes[NodesConditions[i][0]], NodesConditions[i][1], NodesConditions[i][2], NodesConditions[i][3])

    for i in range(0, len(ForcesConditions)):
        nForcesConditions = (map_nodes[ForcesConditions[i][0]], ForcesConditions[i][1], ForcesConditions[i][2], ForcesConditions[i][3])

    return nPosNodes, nElements, nNodesConditions, nForcesConditions


def expand(ex_P, ex_E, newPosNodes, cont_nodes, cont_elements, s, e, n):
    ex_E[cont_elements] = (s, cont_nodes)
    cont_elements += 1
    for i in range(n-1):
        ex_E[cont_elements] = (cont_nodes + i, cont_nodes + i + 1)
        cont_elements += 1
    ex_E[cont_elements] = (cont_nodes+n-1, e)
    cont_elements += 1

    x1 = newPosNodes[s][0]
    x2 = newPosNodes[e][0]
    y1 = newPosNodes[s][1]
    y2 = newPosNodes[e][1]
    
    for i in range(n):
        ex_P[cont_nodes] = (x1+(i+1)*(x2-x1)/(n+1), y1+i*(y2-y1)/(n+1))
        cont_nodes += 1
        

def __main__(PosNodes, Elements, NodesConditions, ForcesConditions, Properties):
    newPosNodes, newElements, newNodesConditions, newForcesConditions = new_model(PosNodes, Elements, NodesConditions, ForcesConditions)
    
    n = 8
    Nodes = 3*(len(newPosNodes)+len(newElements)*n)
    NumberOfElements = len(newElements)*(n+1)
    
    expand_PosNodes = []*(Nodes//3)
    expand_Elements = []*(NumberOfElements)


    for i in range(len(newPosNodes)):
        expand_PosNodes[i] = newPosNodes[i]

    cont_nodes = len(newPosNodes)
    cont_elements = 0

    for i in range(len(Elements)):
        start_node = newElements[i][0]
        end_node = newElements[i][1]
        expand(expand_PosNodes, expand_Elements, newPosNodes, cont_nodes, cont_elements, start_node, end_node, n)
    
    expand_PosNodes = np.array(expand_PosNodes)
    expand_Elements = np.array(expand_Elements)

    L = []*(NumberOfElements)
    K = []*(NumberOfElements)
    
    for i in range(0, NumberOfElements):
        L[i] = fea.DistNodes(expand_PosNodes[expand_Elements[i][0]], expand_PosNodes[expand_Elements[i][1]])
    L = np.array(L)

    for i in range(0, NumberOfElements):
        K[i] = fea.ElementStiffness(L[i][0], L[i][1], Properties[0], Properties[1], Properties[0])
    K = np.array(K)

    StiffnessMatrix = np.zeros((Nodes, Nodes))

    U = np.zeros(Nodes).reshape(Nodes,1)
    F = np.zeros(Nodes).reshape(Nodes,1)

    fea.Initialize(StiffnessMatrix, expand_Elements, K, U, F, Nodes, NumberOfElements)
    boundaries_nodes = []
    loads_nodes = []

    #Node in UBoundary = Node*3+(x=0,y=1,theta=2)
    for i in range(len(newNodesConditions)):
        for j in range(1,4):
            if newNodesConditions[i][j][1]:
                fea.UBoundaryCondition(U, newNodesConditions[i][j][0], 3*newNodesConditions[i][0] + j-1, boundaries_nodes)
    
    for i in range(len(newForcesConditions)):
        for j in range(1,4):
            if newForcesConditions[i][j][1]:
                fea.FBoundaryCondition(F, newForcesConditions[i][j][0], 3*newForcesConditions[i][0] + j-1, loads_nodes)


    U,F = fea.Solve(StiffnessMatrix, U, F, StiffnessMatrix, Nodes, boundaries_nodes)
    print("Stiffness Matrix:\n",StiffnessMatrix,'\n')
    print("Displacements:")
    for i in range((U.shape[0])//3):
        print('[',U[3*i][0],',',U[3*i+1][0],',',U[3*i+2][0],']')
    print("\nForces:")
    for i in range((F.shape[0])//3):
        print('[',F[3*i][0],',',F[3*i+1][0],',',F[3*i+2][0],']')

    return 0, 0