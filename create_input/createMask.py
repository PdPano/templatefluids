import maskUtils
import numpy as np
import constants
import scipy.ndimage.morphology as morph
import flagsAndShifts
import splineslib as slib
import derFlags

def createMaskAndPairs(curves,boundarySpecs,xgrid,ygrid):
    bodyMask=createBodyMask(curves,xgrid,ygrid)
    boundaryMask=createBoundaryMask(boundarySpecs,xgrid,ygrid)
    mask=composeMasks(bodyMask,boundaryMask)
    imagePointsPair = slib.getPointImagePair(mask,curves,xgrid,ygrid)
    removeImagesOutsideDomain(imagePointsPair,mask,xgrid,ygrid)
    derFlags.setDerFlags(mask)
    return mask,imagePointsPair

def createMaskAndDiscontinuities(curves,boundarySpecs,xgrid,ygrid,params):
    bodyMask=createBodyMask(curves,xgrid,ygrid)
    boundaryMask=createBoundaryMask(boundarySpecs,xgrid,ygrid)
    mask=composeMasks(bodyMask,boundaryMask)
    discontinuities,discontinuityPos=createDiscontinuitiesList(curves,xgrid,ygrid,params)
    return mask,discontinuities,discontinuityPos

def createDiscontinuitiesList(curves,xgrid,ygrid,params):
    gam = params['gamma']
    mach = params['mach']
    T_inf = params['T_inf']
    E_inf = T_inf*1.0/(gam*(gam-1)*mach**2)
    posX,derX=maskUtils.getXCrossings(curves,xgrid,ygrid)
    posY,derY=maskUtils.getYCrossings(curves,xgrid,ygrid)
    discontinuityList=[]
    dx=xgrid[1]-xgrid[0]
    dy=ygrid[1]-ygrid[0]
    xmin=xgrid.min()
    ymin=ygrid.min()
    xgrid_trick = np.array(list(xgrid)+[xgrid[-1]+dx])
    ygrid_trick = np.array(list(ygrid)+[ygrid[-1]+dy])
    for point in posX:
        try:
            i = next(v[0] for v in enumerate(ygrid_trick) if v[1]>point[1])-1
            j = next(v[0] for v in enumerate(xgrid_trick) if v[1]>point[0])-1
        except:
            continue
        fracI=(point[1]-ymin)/dy - i
        if(0<=i<len(ygrid)-1 and 0<=j<len(xgrid)):
            discontinuityList.append( ('y',i,j,fracI,0.0,
                                    1.0,0.0,0.0,E_inf,#Point left
                                    1.0,0.0,0.0,E_inf #Point right
                                    ) )
    for point in posY:
        try:
            i = next(v[0] for v in enumerate(ygrid_trick) if v[1]>point[1])-1
            j = next(v[0] for v in enumerate(xgrid_trick) if v[1]>point[0])-1
        except:
            continue
        fracJ=(point[0]-xmin)/dx - j
        if(0<=i<len(ygrid) and 0<=j<len(xgrid)-1):
            discontinuityList.append( ('x',i,j,fracJ,0.0,
                                    1.0,0.0,0.0,E_inf,#Point left
                                    1.0,0.0,0.0,E_inf #Point right
                                    ) )
    discontinuityPos = posX+posY
    return discontinuityList,discontinuityPos


def removeImagesOutsideDomain(imagePointsPair,gridMask,xgrid,ygrid):
    nPointsI,nPointsJ=gridMask.shape
    dx=xgrid[1]-xgrid[0]
    dy=ygrid[1]-ygrid[0]
    xmin=xgrid.min()
    ymin=ygrid.min()
    origKeys=imagePointsPair.keys()
    for ind in sorted(origKeys):
        i = int(np.floor((imagePointsPair[ind][2][1]-ymin)/dy))
        j = int(np.floor((imagePointsPair[ind][2][0]-xmin)/dx))
        neighb = [(i,j), (i+1,j),(i+1,j+1),(i,j+1)]
        if((0<=np.array(neighb)[:,0]).all() and (0<=np.array(neighb)[:,1]).all() and (np.array(neighb)[:,0]<nPointsI).all() and (np.array(neighb)[:,1]<nPointsJ).all()):
            continue
        else:
            gridMask[ind]=1
            imagePointsPair.pop(ind)


def composeMasks(bodyMask,boundaryMask):
    mask=np.zeros_like(bodyMask).astype(int)
    nPointsI,nPointsJ=bodyMask.shape
    for i in range(nPointsI):
        for j in range(nPointsJ):
            bodyType = flagsAndShifts.pointType(bodyMask[i,j])
            boundaryType = flagsAndShifts.pointType(boundaryMask[i,j])
            if(bodyType == 1 or boundaryType == 1): #Solid point
                mask[i,j]=1
            elif(bodyType == 2): #Ghost point
                mask[i,j]=2
            elif(boundaryType==3):#Boundary
                mask[i,j]=boundaryMask[i,j]
    return mask


def createBodyMask(curves,xgrid,ygrid):
    mask=np.zeros((len(ygrid),len(xgrid))).astype(int)
    posX,derX=maskUtils.getXCrossings(curves,xgrid,ygrid)
    posY,derY=maskUtils.getYCrossings(curves,xgrid,ygrid)

    maskBodiesFromCrossings(posX,derX,posY,derY,xgrid,ygrid,mask)
    fillMask(mask)
    return mask

def createBoundaryMask(boundarySpecs,xgrid,ygrid):
    mask=np.zeros((len(ygrid),len(xgrid))).astype(int)

    nPointsI,nPointsJ = mask.shape

    dx = xgrid[1]-xgrid[0]
    dy = ygrid[1]-ygrid[0]

    xmin=xgrid[0]
    ymin=ygrid[0]

    for boundary in boundarySpecs:
        p1,p2,btype,bdict=boundary
        isbx,isby,bxtype,bytype,bflag,bshift,ishift,jshift=constants.boundaryDict[btype]
        bflag=0
        equalX = abs(p1[0]-p2[0])<0.1*dx
        equalY = abs(p1[1]-p2[1])<0.1*dy
        if(equalX):
            j=int((p1[0]-xmin)/dx+0.5) #closest grid line
            i1=int((p1[1]-ymin)/dy+0.5)
            i2=int((p2[1]-ymin)/dy+0.5)
            imin,imax = (i1,i2) if (i1<i2) else (i2,i1)
            for i in range(imin,imax+1):
                if(0<=i<nPointsI and 0<=j<nPointsJ): #Closest grid line receives boundary
                    mask[i,j] = mask[i,j] | (3+(bflag<<bshift))
                if(0<=i+ishift<nPointsI and 0<=j+jshift<nPointsJ):
                    if(mask[i+ishift,j+jshift]==0): #Other side is marked as solid
                        mask[i+ishift,j+jshift]=1
        elif(equalY):
            i=int((p1[1]-ymin)/dy+0.5) #closest grid line
            j1=int((p1[0]-xmin)/dx+0.5)
            j2=int((p2[0]-xmin)/dx+0.5)
            jmin,jmax = (j1,j2) if (j1<j2) else (j2,j1)
            for j in range(jmin,jmax+1):
                if(0<=i<nPointsI and 0<=j<nPointsJ): #Closest grid line receives boundary
                    mask[i,j] = mask[i,j] | (3+(bflag<<bshift))
                if(0<=i+ishift<nPointsI and 0<=j+jshift<nPointsJ):
                    if(mask[i+ishift,j+jshift]==0): #Other side is marked as solid
                        mask[i+ishift,j+jshift]=1
        else:
            print("Warning: not using ",boundary)

    dilationMask=~(mask%4==3)
    image = mask==1
    inner = morph.binary_dilation(image,mask=dilationMask,iterations=2*max(mask.shape)).astype(dilationMask.dtype)
    mask[inner]=1
    return mask


def maskBodiesFromCrossings(posX,derX,posY,derY,xgrid,ygrid,mask):
    nPointsI,nPointsJ = mask.shape

    dx = xgrid[1]-xgrid[0]
    dy = ygrid[1]-ygrid[0]

    xmin=xgrid[0]
    ymin=ygrid[0]

    xmax=xgrid[-1]
    ymax=ygrid[-1]

    for pos,der in zip(posX,derX):
        ifloat = (pos[1]-ymin)/dy
        i=int(np.floor(ifloat))
        jfloat = (pos[0]-xmin)/dx+0.1
        j=int(np.floor(jfloat))

        if(der>0):
            if(0<=i<nPointsI and 0<=j<nPointsJ):
                mask[i,j]+=    1 #outside
            if(0<=i+1<nPointsI and 0<=j<nPointsJ):
                mask[i+1,j]+= -1 #inside
        else:
            if(0<=i<nPointsI and 0<=j<nPointsJ):
                mask[i,j]+=   -1 #inside
            if(0<=i+1<nPointsI and 0<=j<nPointsJ):
                mask[i+1,j]+=  1 #outside

    for pos,der in zip(posY,derY):
        ifloat = (pos[1]-ymin)/dy+0.1
        i=int(np.floor(ifloat))
        jfloat = (pos[0]-xmin)/dx
        j=int(np.floor(jfloat))

        if(der>0):
            if(0<=i<nPointsI and 0<=j<nPointsJ):
                mask[i,j]+=   -1 #inside
            if(0<=i<nPointsI and 0<=j+1<nPointsJ):
                mask[i,j+1]+=  1 #outside
        else:
            if(0<=i<nPointsI and 0<=j<nPointsJ):
                mask[i,j]+=    1 #outside
            if(0<=i<nPointsI and 0<=j+1<nPointsJ):
                mask[i,j+1]+= -1 #inside

def fillMask(mask):
    nPointsI,nPointsJ=mask.shape
    dilationMask=~(mask>0)
    image = mask<0
    inner = morph.binary_dilation(image,mask=dilationMask,iterations=2*max(mask.shape)).astype(dilationMask.dtype)
    mask[inner]=1
    mask[~inner]=0
    for i in range(nPointsI):
        for j in range(nPointsJ):
            if(mask[i,j]):
                keep=True
                for p in [(1,1),(1,0),(1,-1),(0,1),(0,-1),(-1,1),(-1,0),(-1,-1)]:
                    if(0<=i+p[0]<nPointsI and 0<=j+p[1]<nPointsJ):
                        if(keep and 1<=mask[i+p[0],j+p[1]]<=2):
                            keep=True
                        else:
                            keep=False
                if(not keep):
                    mask[i,j]=2
