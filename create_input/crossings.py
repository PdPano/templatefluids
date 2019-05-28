import scipy.interpolate as interp
import scipy.ndimage.morphology as morph
import numpy as np


def getYcrossings(curveX, curveY, xgrid, ygrid):
    splineY, u = interp.splprep((curveX, curveY), s=0, k=3, per=True)
    ty, cy, ky = splineY
    xvalsy, yvalsy = cy

    posY = [[], []]
    derY = [[], []]
    for i, shifty in enumerate(ygrid):
        newcy = ((np.array(xvalsy), np.array(yvalsy)-shifty))
        shiftedSplineY = (ty, newcy, ky)

        rootsY = interp.sproot(shiftedSplineY)

        if(len(rootsY[1])):
            posYtemp = interp.splev(rootsY[1], splineY)

            posY[0] = posY[0]+list(posYtemp[0])
            posY[1] = posY[1]+list(posYtemp[1])

            derYtemp = interp.splev(rootsY[1], splineY, der=1)

            for ind in range(len(derYtemp[0])):
                derYtemp[0][ind], derYtemp[1][ind] = (
                    np.array(
                            (derYtemp[0][ind], derYtemp[1][ind]))
                    / np.sqrt((derYtemp[0][ind]**2+derYtemp[1][ind]**2)))

            derY[0] = derY[0]+list(derYtemp[0])
            derY[1] = derY[1]+list(derYtemp[1])
    return np.array(posY).T, np.array(derY).T


def getXcrossings(curveX, curveY, xgrid, ygrid):
    splineX, u = interp.splprep((curveY, curveX), s=0, k=3, per=True)

    tx, cx, kx = splineX
    xvalsx, yvalsx = cx

    posX = [[], []]
    derX = [[], []]

    for shiftx in xgrid:
        newcx = ((np.array(xvalsx), np.array(yvalsx)-shiftx))
        shiftedSplineX = (tx, newcx, kx)

        rootsX = interp.sproot(shiftedSplineX)

        if(len(rootsX[1])):
            posXtemp = interp.splev(rootsX[1], splineX)

            posX[0] = posX[0]+list(posXtemp[1])
            posX[1] = posX[1]+list(posXtemp[0])

            derXtemp = interp.splev(rootsX[1], splineX, der=1)
            for ind in range(len(derXtemp[0])):
                derXtemp[0][ind], derXtemp[1][ind] = (
                    np.array(
                            (derXtemp[0][ind], derXtemp[1][ind]))
                    / np.sqrt((derXtemp[0][ind]**2+derXtemp[1][ind]**2))
                )
            derX[0] = derX[0]+list(derXtemp[1])
            derX[1] = derX[1]+list(derXtemp[0])
    return np.array(posX).T, np.array(derX).T


def createMask(curves, xgrid, ygrid):
    mask = np.zeros((len(ygrid), len(xgrid)))
    posX = []
    posY = []
    derX = []
    derY = []
    for curve in curves:
        posXtemp, derXtemp = getXcrossings(
            curve[:, 0], curve[:, 1], xgrid, ygrid)
        posYtemp, derYtemp = getYcrossings(
            curve[:, 0], curve[:, 1], xgrid, ygrid)
        posX = np.array(list(posX)+list(posXtemp))
        posY = np.array(list(posY)+list(posYtemp))
        derX = np.array(list(derX)+list(derXtemp))
        derY = np.array(list(derY)+list(derYtemp))
    maskFromCrossings(posX, derX, posY, derY, xgrid, ygrid, mask)
    return mask


def maskFromCrossings(posX, derX, posY, derY, xgrid, ygrid, mask):
    nPointsI = len(ygrid)
    nPointsJ = len(xgrid)

    dx = xgrid[1]-xgrid[0]
    dy = ygrid[1]-ygrid[0]

    xmin = xgrid[0]
    ymin = ygrid[0]

    for pos, der in zip(posX, derX):
        ifloat = (pos[1]-ymin)/dy
        i = int(np.floor(ifloat))
        jfloat = (pos[0]-xmin)/dx+0.1
        j = int(np.floor(jfloat))

        if(der[0] > 0):
            if(0 <= i < nPointsI and 0 <= j < nPointsJ):
                mask[i, j] += 1  # outside
            if(0 <= i+1 < nPointsI and 0 <= j < nPointsJ):
                mask[i+1, j] += -1  # inside
        else:
            if(0 <= i < nPointsI and 0 <= j < nPointsJ):
                mask[i, j] += -1  # inside
            if(0 <= i+1 < nPointsI and 0 <= j < nPointsJ):
                mask[i+1, j] += 1  # outside

    for pos, der in zip(posY, derY):
        ifloat = (pos[1]-ymin)/dy+0.1
        i = int(np.floor(ifloat))
        jfloat = (pos[0]-xmin)/dx
        j = int(np.floor(jfloat))

        if(der[1] > 0):
            if(0 <= i < nPointsI and 0 <= j < nPointsJ):
                mask[i, j] += -1  # inside
            if(0 <= i < nPointsI and 0 <= j+1 < nPointsJ):
                mask[i, j+1] += 1  # outside
        else:
            if(0 <= i < nPointsI and 0 <= j < nPointsJ):
                mask[i, j] += 1  # outside
            if(0 <= i < nPointsI and 0 <= j+1 < nPointsJ):
                mask[i, j+1] += -1  # inside
    dilationMask = ~(mask < 0)
    image = mask > 0
    inner = morph.binary_dilation(image, mask=dilationMask, iterations=max(
        mask.shape)).astype(dilationMask.dtype)
    mask[inner] = 0
    mask[~inner] = 1
    for i in range(nPointsI):
        for j in range(nPointsJ):
            if(mask[i, j]):
                keep = True
                for p in [(1, 1), (1, 0), (1, -1),
                          (0, 1), (0, -1), (-1, 1),
                          (-1, 0), (-1, -1)]:
                    if(0 <= i+p[0] < nPointsI and 0 <= j+p[1] < nPointsJ):
                        if(keep and mask[i+p[0], j+p[1]]):
                            keep = True
                        else:
                            keep = False
                if(not keep):
                    mask[i, j] = 2
