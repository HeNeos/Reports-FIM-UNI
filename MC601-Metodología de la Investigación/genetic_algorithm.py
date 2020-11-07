import numpy as np
from geneticalgorithm import geneticalgorithm as ga
import truss_2D as truss

INF = 1e20

Nodes = 10
Dimension = (Nodes*(Nodes-1))/2
PosNodes = [] # POSITIONS
NodesConditions = [] # NODES WITH RESTRICTIONS (BOUNDS)
ForcesConditions = [] # NODES WITH RESTRICTION (LOADS)
Properties = [] # PROPERTIES OF THE MODEL (A, I, E)

def make_adjacency(X, Elements):
    sz = len(X)
    N = int((np.sqrt(8*sz+1)+1)/2)
    position = 0
    adjacency_matrix = np.eye(N, dtype=int)
    for j in range(0,N-1):
        for i in range(0,N-1-j):
            if(X[position]):
                Elements.append((j,j+i+1))
                adjacency_matrix[j][j+i+1] = 1
                adjacency_matrix[j+i+1][j] = 1
                position += 1
    return adjacency_matrix

def count_nodes(N, F):
    nodes = set()
    for i in range(len(N)):
        nodes.add(N[i][0])
    for i in range(len(F)):
        nodes.add(F[i][0])
    return nodes


def dfs(g, s):
    return 0

def f(X):
    # PASS X AS A (MATRIX -> VECTOR) ADJACENCY MATRIX?
    Elements = []
    adjacency = make_adjacency(X, Elements)

    # VERIFY IF X IS A CONNECTED GRAPH WITH THE BOUNDARIES
    # REMEMBER, IN A FEM MODEL YOU HAVE SOME FIXED ELEMENTS: LOADS AND BOUNDARIES
    # THE GENERAL NODES ARE JUST A REPRESENTATION OF THE GEOMETRY
    
    fixed_nodes = count_nodes(NodesConditions, ForcesConditions)

    for i in fixed_nodes:
        touched = 1 + dfs(adjacency, i)
        if touched != len(fixed_nodes):
            return INF
    
    # RUN FINITE ELEMENT MODEL
    MaxStress, Mass = truss.__main__(PosNodes, Elements, NodesConditions, ForcesConditions, Properties)

    # RETURN cost(MAX STRESS VALUE, MASS)
    return (MaxStress+Mass) # :ppp

model = ga(function = f,
        dimension = Dimension,
        variable_type = 'bool',
        function_timeout = 10)



