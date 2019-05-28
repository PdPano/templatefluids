import numpy as np
import flagsAndShifts
import constants
def printImmersedBoundary(imagePointsPair,gridMask,xgrid,ygrid,fname="./results/immersedInterface.dat"):
    nPointsI,nPointsJ=gridMask.shape
    X,Y = np.meshgrid(xgrid,ygrid)
    dx = xgrid[1]-xgrid[0]
    try:
        dy = ygrid[1]-ygrid[0]
    except:
        dy=0.1
    xmin = xgrid.min()
    ymin = ygrid.min()
    with open(fname,"w") as fout:
        fout.write("GHIAS\n")
        fout.write("%d\n"%(len(list(imagePointsPair.keys())),))
        for ind in sorted(imagePointsPair.keys()):
            i = int(np.floor((imagePointsPair[ind][2][1]-ymin)/dy))
            j = int(np.floor((imagePointsPair[ind][2][0]-xmin)/dx))
            neighb = [(i,j), (i+1,j),(i+1,j+1),(i,j+1)]
            if((0<=np.array(neighb)[:,0]).all() and (0<=np.array(neighb)[:,1]).all() and (np.array(neighb)[:,0]<nPointsI).all() and (np.array(neighb)[:,1]<nPointsJ).all()):
                fout.write("%d %d\n"%(ind[0],ind[1]))
                for n in neighb:
                    if n in imagePointsPair:
                        bound = imagePointsPair[n][1]
                        norm = imagePointsPair[n][3]
                        fout.write("b %d %d %1.15e %1.15e %1.15e %1.15e\n"%(n[0],n[1],bound[0],bound[1],norm[0],norm[1]))
                    else:
                        fout.write("f %d %d %1.15e %1.15e %1.15e %1.15e\n"%(n[0],n[1],X[n],Y[n],0.0,0.0))
                imageCoord = imagePointsPair[ind][2]
                fout.write("%1.15e %1.15e\n\n"%(imageCoord[0],imageCoord[1]))
            else:
                gridMask[ind]=1

def printImmersedDiscontinuities(discontinuities,fname="./results/immersedInterface.dat"):
    with open(fname,"w") as fout:
        fout.write("KARAGIOZIS\n")
        fout.write("%d\n"%(len(discontinuities),))
        for point in discontinuities:
            fout.write("%c %d %d %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n"%(point))
def showImmersedBoundary(imagePointsPair,gridMask,xgrid,ygrid):
    nPointsI,nPointsJ=gridMask.shape
    X,Y = np.meshgrid(xgrid,ygrid)
    dx = xgrid[1]-xgrid[0]
    dy = ygrid[1]-ygrid[0]
    xmin = xgrid.min()
    ymin = ygrid.min()
    for ind in sorted(imagePointsPair.keys()):
        print("Ind:   ",ind)
        print("Coord: ",imagePointsPair[ind][0])
        print("Bound: ",imagePointsPair[ind][1])
        print("Image: ",imagePointsPair[ind][2])
        print("Norml: ",imagePointsPair[ind][3])
        i = int(np.floor((imagePointsPair[ind][2][1]-ymin)/dy))
        j = int(np.floor((imagePointsPair[ind][2][0]-xmin)/dx))
        print("Neigh: ", (i,j), (i+1,j),(i+1,j+1),(i,j+1))
        neighb = [(i,j), (i+1,j),(i+1,j+1),(i,j+1)]

        if((0<=np.array(neighb)[:,0]).all() and (0<=np.array(neighb)[:,1]).all() and (np.array(neighb)[:,0]<nPointsI).all() and (np.array(neighb)[:,1]<nPointsJ).all()):
            for n in neighb:
                if n in imagePointsPair:
                    print('b',imagePointsPair[n][1],imagePointsPair[n][3])
                else:
                    print('f',(X[n],Y[n]))
        else:
            gridMask[ind]=1
        print(" ")


def printFlagMap(flagMap,xgrid,ygrid,fname="./results/gridInfo.dat"):
    xmin = xgrid.min()
    ymin = ygrid.min()
    dx=xgrid[1]-xgrid[0]
    try:
        dy=ygrid[1]-ygrid[0]
    except:
        dy=0.1
    ni,nj = flagMap.shape
    with open(fname,"w") as fout:
        fout.write("%d %d\n"%(ni,nj))
        fout.write("%f %f\n"%(dx,dy))
        fout.write("%f %f\n"%(xmin,ymin))
        for i in range(ni):
            for j in range(nj):
                fout.write("%d "%(flagMap[i,j],))
            fout.write("\n")




def printBoundaryConfiguration(boundaryDictionary,fname="./results/boundary.dat"):
    fout=open(fname,'w')
    fout.write(str(len(boundaryDictionary))+"\n")

    for boundary in boundaryDictionary:
        isbx=boundaryDictionary[boundary]['isbx']
        isby=boundaryDictionary[boundary]['isby']
        bxtype=boundaryDictionary[boundary]['bxtype']
        bytype=boundaryDictionary[boundary]['bytype']
        rho=boundaryDictionary[boundary]['rho']
        u=boundaryDictionary[boundary]['u']
        v=boundaryDictionary[boundary]['v']
        T=boundaryDictionary[boundary]['T']
        p=boundaryDictionary[boundary]['p']
        timeFunctionType=boundaryDictionary[boundary]['timeFunctionType']
        timeParameterA=boundaryDictionary[boundary]['timeParameterA']
        timeParameterB=boundaryDictionary[boundary]['timeParameterB']
        fout.write("{:d}\n".format(boundary))
        fout.write("{:d} {:d}\n".format(isbx,isby))
        fout.write("{:d} {:d}\n".format(bxtype,bytype))
        fout.write("{:1.15e} {:1.15e} {:1.15e} {:1.15e} {:1.15e}\n".format(rho,u,v,T,p))
        fout.write("{:d} {:1.15e} {:1.15e}\n\n".format(timeFunctionType,timeParameterA,timeParameterB))

        #fout.write("{:d} {:1.15e} {:1.15e} {:1.15e} {:1.15e} {:1.15e} {:d} {:1.15e} {:1.15e}\n".format(boundary,rho,u,v,T,p,timeFunctionType,timeParameterA,timeParameterB))

def printInitialCondition(gridMask,params,fname="./results/initialCondition.dat"):
    gam = params['gamma']
    mach = params['mach']
    T_inf = params['T_inf']
    E_inf = T_inf*1.0/(gam*(gam-1)*mach**2)
    nPointsI,nPointsJ = gridMask.shape
    fout=open(fname,'w')
    fout.write("{:d} {:d}\n".format(nPointsI,nPointsJ))
    if('init' in params):
        if(params['init']=='explosion'):
            centerI,centerJ=int(nPointsI/2),int(nPointsJ/2)
            if('explosion_pos' in params):
                centerI,centerJ=params['explosion_pos']
            radius=15
            for i in range(nPointsI):
                for j in range(nPointsJ):
                    r2 = (i-centerI)**2 + (j-centerJ)**2
                    if(r2<radius**2):
                        rho=10
                        e=16*E_inf
                    else:
                        rho=1
                        e=E_inf
                    fout.write("{:1.15e} {:1.15e} {:1.15e} {:1.15e}\n".format(rho,0.,0.,e))
            return
    for i in range(nPointsI):
        for j in range(nPointsJ):
            fout.write("{:1.15e} {:1.15e} {:1.15e} {:1.15e}\n".format(1.0,0.0,0.0,E_inf))

def printSplinesToFile(pointsFromSplines,fname="./results/splinesPoints.csv"):
    with open(fname,"w") as outfile:
        outfile.write("x,y\n")
        for x,y in pointsFromSplines:
            if(x is not None):
                outfile.write("{:1.15e},{:1.15e}\n".format(x,y))

