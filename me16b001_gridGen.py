"""
Reference(s):
http://caefn.com/cfd/hyperbolic-tangent-stretching-grid
"""
from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from array import *
from scipy.sparse import *

# TODO: fix the grid function for odd value of grid elements
def grid(nx,ny, gama):
    # use ty tx to change number of elements
    tx=(2)/((nx+1)+1)
    ty=(2)/((ny+1)+1)
    x=[]
    y=[]
    nx=[]
    ny=[]

    # x elements on left half
    for i in np.arange(0., 1., tx):
        nx.append(i)
    for j in nx:
        x.append(1-(np.tanh(gama*(1-(2*j)/len(nx))))/(np.tanh(gama)))

    # y elements on right half
    for i in np.arange(0., 1., ty):
        ny.append(i)
    for i in ny:
        y.append(1-(np.tanh(gama*(1-(2*i)/len(ny))))/(np.tanh(gama)))

    # mirroring x and y elements for the right half
    for i in range(len(nx)-1):
        x.append(x[len(nx)+i-1]+x[len(nx)-(i+1)]-x[len(nx)-(i+2)])
    for i in range(len(ny)-1):
        y.append(y[len(ny)+i-1]+y[len(ny)-(i+1)]-y[len(ny)-(i+2)])

    xd=[]
    yh=[]
    for i in x:
        xd.append(i/(x[len(x)-1]))
    for i in y:
        yh.append(i/(y[len(y)-1]))
    return(xd,yh)

def gridPlot(c,d):
    for i in c:
        for k in d:
            plt.scatter(i,k,color='black',marker='+')
    plt.show()

def spacing(test_list):
    res = [test_list[i + 1] - test_list[i] for i in range(len(test_list)-1)]
    return(res)

def uGrid(nx,ny):
    l1=[]
    l2=[]
    for i in range(nx+1):
        l1.append(1/(nx+1))
    for i in range(ny+1):
        l2.append(1/(ny+1))
    return(l1,l2)

# # print(uGrid(10,10))
# print(grid(10,10,3))
# # NOTE: call from any program using
# # from gridGen import grid,gridPlot
# x,y=grid(10,10,3)   # accepts (nx, ny, game)
# gridPlot(grid(10,10,3)[0],grid(10,10,3)[1])
