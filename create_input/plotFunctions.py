import matplotlib.pyplot as plt
import numpy as np

def plotEverything(gridMask,imagePointsPair,pointsFromSplines,xgrid,ygrid,discontinuityPos):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    xpairs=[]
    ypairs=[]
    imagePoints=[]
    boundaryPoints=[]
    for i in imagePointsPair:
        xpairs.append([imagePointsPair[i][0][0],imagePointsPair[i][2][0]])
        ypairs.append([imagePointsPair[i][0][1],imagePointsPair[i][2][1]])
        imagePoints.append((imagePointsPair[i][2][0],imagePointsPair[i][2][1]))
        boundaryPoints.append((imagePointsPair[i][1][0],imagePointsPair[i][1][1]))

    imagePoints=np.array(imagePoints)
    boundaryPoints=np.array(boundaryPoints)

    xlist = []
    ylist = []
    for xends,yends in zip(xpairs,ypairs):
        xlist.extend(xends)
        xlist.append(None)
        ylist.extend(yends)
        ylist.append(None)


    for x1,y1 in zip(xpairs,ypairs):
        ax.plot(x1,y1,'k-',alpha=0.3)

    X,Y = np.meshgrid(xgrid,ygrid)
    meshPoints = np.dstack([X.ravel(),Y.ravel()])[0]
    ax.scatter(meshPoints[:,0],meshPoints[:,1],marker="+",s=30,color="k", label="Mesh points")

    if(len(pointsFromSplines)>0):
        ax.plot(pointsFromSplines[:,0],pointsFromSplines[:,1],lw=2,color='b', label="Interpolating cubic spline")
    ax.scatter(X[gridMask%4==3],Y[gridMask%4==3],marker="o",s=50,lw=2,facecolors="none",edgecolor="b",label="Boundary points")
    ax.scatter(X[gridMask%4==2],Y[gridMask%4==2],marker="o",s=50,lw=2,facecolors="none",edgecolor="r",label="Ghost points")
    ax.scatter(X[gridMask%4==1],Y[gridMask%4==1],marker="o",s=50,lw=2,facecolors="none",edgecolor="y",label="Solid points")
    if(len(imagePoints)>0):
        ax.scatter(imagePoints[:,0],imagePoints[:,1],marker="o",s=50,lw=2,facecolors="none",edgecolor="c",label="Image points")
    if(len(boundaryPoints)>0):
        ax.scatter(boundaryPoints[:,0],boundaryPoints[:,1],marker="o",s=30,lw=1,facecolors="k",edgecolor="k",label="Interface points")
    ax.set_aspect('equal')
    ax.set_xlabel("x",fontsize=15)
    ax.set_ylabel("y",fontsize=15)

    discontinuityPos = np.array(discontinuityPos)
    ax.scatter(discontinuityPos[:,0],discontinuityPos[:,1],marker="x",color='r')
    #ax.legend(loc="lower right")

    plt.show()
