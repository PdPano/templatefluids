import numpy as np
import math
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import constants


def splinesFromPoints(curves):
    splines=[]
    for points in curves:
        points = np.array(points)
        nPoints = points.shape[0]
        matrix = np.eye(nPoints,k=-1) + 4*np.eye(nPoints,k=0)+np.eye(nPoints,k=1)
        matrix[0,-1]=1
        matrix[-1,0]=1

        difs = 3*(np.roll(points,-1,axis=0)-np.roll(points,1,axis=0))

        Ds = np.array(np.linalg.solve(matrix,difs))

        coefs = np.zeros((nPoints,2,4))

        #Polynomial coefficients: a*t^0 + b*t^1 + c*t^2 + d*t^3
        coefs[:,:,0] = points[:,:] #a
        coefs[:,:,1] = Ds[:,:]     #b
        coefs[:,:,2] = 3*(np.roll(points,-1,axis=0)-points) - 2*Ds - np.roll(Ds,-1,axis=0) #c
        coefs[:,:,3] = 2*(points-np.roll(points,-1,axis=0))+Ds+np.roll(Ds,-1,axis=0) #d

        splines.append([(poly.Polynomial(coefs[i,0,:]),poly.Polynomial(coefs[i,1,:])) for i in range(nPoints)])

    return splines

def getPointsFromSplines(splines,amplification):
    if amplification<1:
        amplification=1

    ts = np.arange(0,1+0.5/amplification,1.0/amplification)
    pointsFromSplines = []
    for body in splines:
        for s in body:
            for t in ts:
                pointsFromSplines.append([s[0](t),s[1](t)])
        pointsFromSplines.append([None,None])
    return np.array(pointsFromSplines)



def getPointImagePair(gridMask,curves,xgrid,ygrid):
    k = np.where(gridMask%8==2)
    innerPointsIndex = [ (i,j) for i,j in zip(k[0],k[1])]
    pointImagePair = {}
    for ipoint in innerPointsIndex:
        i,j=ipoint
        point = (xgrid[j],ygrid[i])
        images=[]
        for splines in curves:
            for s in splines:
                scalarProd = ((s[0] - point[0])*s[0].deriv()+(s[1] - point[1])*s[1].deriv())
                r = [ v.real for v in scalarProd.roots() if abs(v.imag)<1e-15 and 0<=v.real<=1]
                if len(r)>0:
                    for x in r:
                        images.append(np.array((s[0](x),s[1](x),x)))
        if len(images)>0:
            dist = [ math.hypot(im[0]-point[0],im[1]-point[1]) for im in images]
            minDistanceIndex = dist.index(min(dist))
            normal = images[minDistanceIndex][:2]-point
            normal = normal/np.sqrt((normal[0]**2+normal[1]**2))
            pointImagePair[ ipoint ] = [point, images[minDistanceIndex][:2], 2*images[minDistanceIndex][:2]-point, normal]
    return pointImagePair
