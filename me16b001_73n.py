from numpy import *
import scipy.linalg
from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt
from gridGen import grid,gridPlot,spacing,uGrid

def plotContour(matrix):
    cs = plt.contourf(matrix,  extend='both')
    cs.cmap.set_over('red')
    cs.cmap.set_under('blue')
    cs.changed()
    plt.show()

def initialize():
    return(zeros((nx+1,ny+1)),zeros((nx+1,ny+1)),zeros((nx+1,ny+1)),zeros((nx+1,ny+1)))

def bcsApply(wMat,psiMat,u,v):
    # applying bc on omega
    # wMat[0,:]=-2*psiMat[1,:]*(r**2)/(dx**2)
    # wMat[nx,:]=-2*psiMat[nx,:]*(r**2)/(dx**2)
    # wMat[:,1]=-2*psiMat[:,1]/(dy**2)
    # wMat[:,nx]=-2*psiMat[:,nx]/(dy**2)-2*u[:,nx]/dy

    # bcs as per the book
    for i in range(ny+1):
            wMat[0,i] = -((2*U)/dy)
    wMat[0,:]=-(3/(dy**2))*psiMat[1,:]-0.5*wMat[1,:]-3*U/dy    # top wall
    wMat[ny,:]=-(3/(dy**2))*psiMat[ny-1,:]-0.5*wMat[ny-1,:]     # bottom wall
    wMat[:,nx]=-(3*r*r/(dx**2))*psiMat[:,nx-1]-0.5*wMat[:,nx-1]    # right wall
    wMat[:,0]=-(3*r*r/(dx**2))*psiMat[:,1]-0.5*wMat[:,1]     # left wall

    # applying BC, psi is zero at all the boundaries
    psiMat[0,:]=0
    psiMat[nx,:]=0
    psiMat[:,0]=0
    psiMat[:,ny]=0
    return(wMat,psiMat)

def bcwMat(wMat,psiMat):
    wMat[0,:]=-(3/(dy**2))*psiMat[1,:]-0.5*wMat[1,:]-3*U/dy    # top wall
    wMat[ny,:]=-(3/(dy**2))*psiMat[ny-1,:]-0.5*wMat[ny-1,:]     # bottom wall
    wMat[:,nx]=-(3*r*r/(dx**2))*psiMat[:,nx-1]-0.5*wMat[:,nx-1]    # right wall
    wMat[:,0]=-(3*r*r/(dx**2))*psiMat[:,1]-0.5*wMat[:,1]     # left wall
    return(wMat)

def calculation():
    a,b,u,v = initialize()
    u[0,:]=1  #  top layer with vel = U
    wMat,psiMat = bcsApply(a,b,u,v)
    print("\nBC applied\n wMat",wMat,"\n psiMat\n",psiMat)

    pItMax=1000
    pIt=0
    while pIt<pItMax:
        # vorticity stream function relation
        for i in range(1,ny):
            for j in range(1,nx):
                wMatOld = wMat
                RHS=(1/Re)*(((wMat[i+1,j]+wMat[i-1,j])/(dx*dx)) + (1/(r*r))*((wMat[i,j+1]+wMat[i,j-1])/(dy*dy)))
                # 1st order upwind over one inner layer
                if i==1 or i==ny-1 or j==1 or j==nx-1:
                    if u[i,j]<=0 and v[i,j]>=0:
                        wMat[i,j]=((-u[i,j]/dx + v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS-(u[i,j]/dx)*wMat[i+1,j]+(v[i,j]/dy)*wMat[i,j-1])
                    elif u[i,j]<=0 and v[i,j]<=0:
                        wMat[i,j]=((-u[i,j]/dx - v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS-(u[i,j]/dx)*wMat[i+1,j]-(v[i,j]/dy)*wMat[i,j+1])
                    elif u[i,j]>=0 and v[i,j]<=0:
                        wMat[i,j]=((u[i,j]/dx - v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS+(u[i,j]/dx)*wMat[i-1,j]-(v[i,j]/dy)*wMat[i,j+1])
                    elif u[i,j]>=0 and v[i,j]>=0:
                        wMat[i,j]=((u[i,j]/dx + v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS+(u[i,j]/dx)*wMat[i-1,j]+(v[i,j]/dy)*wMat[i,j-1])
                # 2nd oder upwind for interor nodes
                else:
                    if u[i,j]<=0 and v[i,j]>=0:
                        wMat[i,j]=((-1.5*u[i,j]/dx + 1.5*v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS-(2*u[i,j]/dx)*wMat[i+1,j]+(0.5*u[i,j]/dx)*wMat[i+2,j]+(2*v[i,j]/dy)*wMat[i,j-1]-(0.5*v[i,j]/dy)*wMat[i,j-2])
                    elif u[i,j]<=0 and v[i,j]<=0:
                        wMat[i,j]=((-1.5*u[i,j]/dx - 1.5*v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS-(2*u[i,j]/dx)*wMat[i+1,j]+(0.5*u[i,j]/dx)*wMat[i+2,j]-(2*v[i,j]/dy)*wMat[i,j+1]+(0.5*v[i,j]/dy)*wMat[i,j+2])
                    elif u[i,j]>=0 and v[i,j]<=0:
                        wMat[i,j]=((1.5*u[i,j]/dx - 1.5*v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS+(2*u[i,j]/dx)*wMat[i-1,j]-(0.5*u[i,j]/dx)*wMat[i-2,j]-(2*v[i,j]/dy)*wMat[i,j+1]+(0.5*v[i,j]/dy)*wMat[i,j+2])
                    elif u[i,j]>=0 and v[i,j]>=0:
                        wMat[i,j]=((1.5*u[i,j]/dx + 1.5*v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS+(2*u[i,j]/dx)*wMat[i-1,j]-(0.5*u[i,j]/dx)*wMat[i-2,j]+(2*v[i,j]/dy)*wMat[i,j-1]-(0.5*v[i,j]/dy)*wMat[i,j-2])
                # print("old",wMatOld)
                # print("new",wMat)
                rs1 =np.amax(abs(wMatOld-wMat))
                print("rs1 ",rs1)
                print("___________________________________")

        # stream function equation
        for i in range(1,nx):
            for j in range(1,ny):
                psiMatOld=psiMat
                psiMat[i,j] = 0.5*((1/(dx**2) + (1/(r**2))*(1/dy*dy))**(-1))*(wMat[i,j]+(psiMat[i+1,j]+psiMat[i-1,j])/(dx**2)+(1/(r**2))*(psiMat[i,j+1]+psiMat[i,j-1])/(dx**2))
                rs2 = np.amax(abs(psiMatOld-psiMat))
                print("###################################")
                print("rs2",rs2)
        # # break if resiudal is close to zero
        # if abs(rs1)<1 and abs(rs2)<1:
        #     break
        print("pIt",pIt)

        # find in the velocities u,v from psiMat
        for i in range(1,ny):
            for j in range(1,nx):
                u[i,j] = (psiMat[i,j+1]-psiMat[i,j-1])/(2*dy)
                v[i,j] = -(psiMat[i+1,j]-psiMat[i-1,j])/(2*dx)
        print("\nu",u,"\nv",v)

        # bcs on omega
        wMat=bcwMat(wMat,psiMat)

        pIt=pIt+1 # counter for psiMat Iterations
    print("\npsiMat at ",pIt," Iteration ",psiMat)
    plotContour(flipud(psiMat))

nx=10  # elements in x dir
ny=10  # elements in y dir
U = 1  # mormalized plate velocity
Re = 100 # reynolds number
r = 1 # aspect ratio
H=1
D=1
dx=H/nx
dy=D/ny
calculation()
