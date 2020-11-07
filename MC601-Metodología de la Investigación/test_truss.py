import numpy as np
import matplotlib.pyplot as plt
import truss_2D as truss
PosNodes = [(0,0), (0,1500), (1500,0), (1500,3000)]
Elements = [(0,1), (0,2), (0,3), (1,3), (2,3)]
Properties = [(0.25*np.pi*(50)**2),(np.pi*50**4)/64,3.2e5]
NodesConditions = [(0,(0,1),(0,1),(0,1)),(1,(0,1),(0,1),(0,1))]
LoadsConditions = [(2,(5000,1),(1000,1),(5000,1))]

truss.__main__(PosNodes, Elements, NodesConditions, LoadsConditions, Properties)

