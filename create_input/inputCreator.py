"""main part"""
import numpy as np
import splineslib as slib
import createMask
import plotFunctions
import outputFunctions
import createBoundaryDictionary

pointsParams = np.arange(0, 2*np.pi, 0.1)
points = np.array([(1.202*np.cos(t)+4.8,
                    2.007512502*np.sin(t)-0.0) for t in pointsParams])
points2 = np.array([(1.202*np.cos(t)+1.8,
                     1.007512502*np.sin(t)-2.0) for t in pointsParams])
points = np.array([(0.552*np.cos(t)+10.0, 0.552*np.sin(t)+0.0)
                   for t in pointsParams])

curvePoints = [points, points2]
curvePoints = [points]
#curvePoints = []
curves = slib.splinesFromPoints(curvePoints)
pointsFromSplines = slib.getPointsFromSplines(curves, 10)

outputFunctions.printSplinesToFile(pointsFromSplines)

nx = 401
ny = 1

# xgrid = np.linspace(-4,4,nx)
# ygrid = np.linspace(-2.5,2.5,ny)
# xgrid = np.linspace(-1,1,nx)
# ygrid = np.linspace(-1,1,ny)
xmin, xmax = (0.0, 10.0)
ymin, ymax = (0, 1)
xgrid = np.linspace(xmin, xmax, nx)
ygrid = np.linspace(ymin, ymax, ny)


boundarySpecs = [[(-1, -1), (-1, 1), 'SubsonicInletLeft',
                  {'u': 'parabola',
                   'timeFunctionType': 'linearRamp',
                   'timeParameterA': 2.0}],
                 [(-1, 1), (1, 1), 'AdiabaticNoSlipWallTop', {}],
                 [(1, 1), (1, -0.5), 'AdiabaticNoSlipWallRight', {}],
                 [(0.5, -0.5), (1, -0.5), 'AdiabaticNoSlipWallBottom', {}],
                 [(0.5, -1), (0.5, -0.5), 'SubsonicOutletRight', {}],
                 [(-1, -1), (0.5, -1), 'AdiabaticNoSlipWallBottom', {}]]

boundarySpecs = [[(-1, -1), (-1, 1), 'SubsonicInletLeft',
                  {'u': 'parabola',
                   'timeFunctionType': 'linearRamp',
                   'timeParameterA': 10.0}],
                 [(-1, 1), (1, 1), 'AdiabaticNoSlipWallTop', {}],
                 [(1, 1), (1, -1), 'SubsonicOutletRight', {}],
                 [(-1, -1), (1, -1), 'AdiabaticNoSlipWallBottom', {}]]

boundarySpecs = [[(xmin, ymin), (xmin, ymax), 'SubsonicInletLeft',
                 {'u': 1,
                  'timeFunctionType': 'linearRamp',
                  'timeParameterA': 100.0}],
                 [(xmin, ymax), (xmax, ymax), 'SupersonicOutletTop', {}],
                 [(xmax, ymax), (xmax, ymin), 'SubsonicOutletRight', {}],
                 [(xmin, ymin), (xmax, ymin), 'SupersonicOutletBottom', {}]]

# boundarySpecs=[[(xmin,ymin),(xmin,ymax),'AdiabaticNoSlipWallLeft',{}],
#               [(xmin,ymax),(xmax,ymax),'AdiabaticNoSlipWallTop',{}],
#               [(xmax,ymax),(xmax,ymin),'AdiabaticNoSlipWallRight',{}],
#               [(xmin,ymin),(xmax,ymin),'AdiabaticNoSlipWallBottom',{}]]

boundarySpecs = [[(xmin, ymin), (xmin, ymax), 'SupersonicInletLeft',
                 {'u': 5/4.,
                  'rho': 8/3.,
                  'T':3*np.sqrt(3)/4,
                  'timeFunctionType': 'linearRamp',
                  'timeParameterA': 0.0}],
                 ]

# params={'gamma':1.4,'mach':3.2,'T_inf':1.0,'init':'explosion','explosion_pos':(int(2*nx/3),int(3*ny/4))}
params = {'gamma': 1.4, 'mach': 1.0, 'T_inf': 1.0}
mask, discontinuities, discontinuityPos = (
        createMask.createMaskAndDiscontinuities(curves, boundarySpecs,
                                                xgrid, ygrid, params))
mask2, discontinuities2, discontinuityPos2 = (
        createMask.createMaskAndDiscontinuities([], boundarySpecs,
                                                xgrid, ygrid, params))
outputFunctions.printImmersedDiscontinuities(
        discontinuities, fname="./results/immersedInterface2.dat")
gridMask, imagePointsPair = createMask.createMaskAndPairs(curves,
                                                          boundarySpecs,
                                                          xgrid, ygrid)
gridMask2, imagePointsPair2 = createMask.createMaskAndPairs([],
                                                            boundarySpecs,
                                                            xgrid, ygrid)
print(gridMask)
boundaryDictionary = createBoundaryDictionary.processBoundarySpecs(
                                                                 boundarySpecs,
                                                                 gridMask,
                                                                 xgrid, ygrid,
                                                                 params=params)
outputFunctions.printImmersedBoundary(imagePointsPair, gridMask, xgrid, ygrid)
outputFunctions.printFlagMap(gridMask, xgrid, ygrid)
outputFunctions.printFlagMap(gridMask2, xgrid, ygrid,
                             fname="./results/gridInfo2.dat")
outputFunctions.printBoundaryConfiguration(boundaryDictionary)
outputFunctions.printInitialCondition(gridMask, params=params)

plotFunctions.plotEverything(gridMask, imagePointsPair, pointsFromSplines,
                             xgrid, ygrid, discontinuityPos)
