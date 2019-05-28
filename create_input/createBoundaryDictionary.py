import flagsAndShifts
import constants
import numpy as np
import collections

def processBoundarySpecs(boundarySpecs,gridMask,xgrid,ygrid,params={'gamma':1.4,'mach':1.0,'T_inf':1.0}):
    nPointsI,nPointsJ = gridMask.shape

    dx = xgrid[1]-xgrid[0]
    dy = ygrid[1]-ygrid[0]

    xmin=xgrid[0]
    ymin=ygrid[0]

    T_inf = params['T_inf']
    machNum = params['mach']

    boundaryDictionary = collections.defaultdict(lambda:{
                                                         'isbx':0,
                                                         'isby':0,
                                                         'bxtype':0,
                                                         'bytype':0,
                                                         'rho':1.0,
                                                         'u':0.0,
                                                         'v':0.0,
                                                         'T':T_inf,
                                                         'p':T_inf*1.0/(params['gamma']*machNum**2),
                                                         'timeFunctionType':0,
                                                         'timeParameterA':1.0,
                                                         'timeParameterB':0.0})
    for boundary in boundarySpecs:
        p1,p2,btype,bdict=boundary
        isbx,isby,bxtype,bytype,bflag,bshift,ishift,jshift=constants.boundaryDict[btype]
        equalX = abs(p1[0]-p2[0])<0.1*dx
        equalY = abs(p1[1]-p2[1])<0.1*dy
        if(equalX):
            j=int((p1[0]-xmin)/dx+0.5) #closest grid line
            i1=int((p1[1]-ymin)/dy+0.5)
            i2=int((p2[1]-ymin)/dy+0.5)
            stepx = 1 if (i1<i2) else -1
            for i in range(i1,i2+stepx,stepx):
                ind = i*nPointsJ+j
                pmin=min(i1,i2)*dy+ymin
                pmax=max(i1,i2)*dy+ymin
                pos =i*dy+ymin
                if(flagsAndShifts.pointType(gridMask[i,j])==3):
                    boundaryDictionary[ind]=fillBoundary(isbx,isby,bxtype,bytype,btype,bdict,pmin,pmax,pos,params,boundaryDictionary[ind])
        elif(equalY):
            i=int((p1[1]-ymin)/dy+0.5) #closest grid line
            j1=int((p1[0]-xmin)/dx+0.5)
            j2=int((p2[0]-xmin)/dx+0.5)
            stepy = 1 if (j1<j2) else -1
            for j in range(j1,j2+stepy,stepy):
                ind = i*nPointsJ+j
                pmin=min(j1,j2)*dx+xmin
                pmax=max(j1,j2)*dx+xmin
                pos =j*dx+xmin
                if(flagsAndShifts.pointType(gridMask[i,j])==3):
                    boundaryDictionary[ind]=fillBoundary(isbx,isby,bxtype,bytype,btype,bdict,pmin,pmax,pos,params,boundaryDictionary[ind])
    return boundaryDictionary

def fillBoundary(isbx,isby,bxtype,bytype,btype,bdict,pmin,pmax,pos,params,ret):
    ret['isbx']=ret['isbx']+isbx
    ret['isby']=ret['isby']+isby
    ret['bxtype']=ret['bxtype']+bxtype
    ret['bytype']=ret['bytype']+bytype
    #Parameters that only 'Inlet' can set
    if('Inlet' in btype):
        if('rho' in bdict):
            ret['rho']=bdict['rho']
        if('T' in bdict):
            ret['T']=bdict['T']
        if('timeFunctionType' in bdict):
            ret['timeFunctionType']=convertTimeFunction(bdict['timeFunctionType'])
        if('timeParameterA' in bdict):
            ret['timeParameterA']=bdict['timeParameterA']
        if('timeParameterB' in bdict):
            ret['timeParameterB']=bdict['timeParameterB']

    #Parameters that only 'Outlet' can set
    if('Outlet' in btype):
        if('p' in bdict):
            ret['p']=bdict['p']

    #Parameters that both 'Wall' and 'Inlet' can set
    if ('Inlet' in btype) or ('Wall' in btype):
        if('u' in bdict):
            ret['u']=getValueFromProfile(bdict['u'],pmin,pmax,pos)
        if('v' in bdict):
            ret['v']=getValueFromProfile(bdict['v'],pmin,pmax,pos)
    return ret

def getValueFromProfile(profile,pmin,pmax,pos,profileParam=0.1):
    if (type(profile) is float) or (type(profile) is int):
        return float(profile)
    if (type(profile) is str):
        if 'exponential' in profile:
            return 1.0*( (1-np.exp(-(pos-pmin)/profileParam))
                        +(1-np.exp(+(pos-pmax)/profileParam))
                        -(1-np.exp(-(pmax-pmin)/profileParam)))
        elif 'parabola' in profile:
            return -4.0/(pmax-pmin)**2 * (pos-pmax) * (pos-pmin)
        else:
            print("Profile",profile,"not found")
            return 0.0

    print("createBoundaryDictionary.py/getValueFromProfile: Should not have been here!")
    return 0.0

def convertTimeFunction(timeFunctionName):
    if timeFunctionName in constants.convertTimeFunctionNameToFlag:
        return constants.convertTimeFunctionNameToFlag[timeFunctionName]
    else:
        print("Couldn't find",timeFunctionName,"in conversion table!")
        return -1
