def getXCrossings(curves,xgrid,ygrid):
    posX=[]
    derX=[]
    for splines in curves:
        for s in splines:
            for x in xgrid:
                r = [ v.real for v in (s[0]-x).roots() if abs(v.imag)<1e-15 and 0<=v.real<=1]
                if len(r)>0:
                    for t in r:
                        posX.append((x,s[1](t)))
                        derX.append((s[0].deriv())(t))
    return posX,derX

def getYCrossings(curves,xgrid,ygrid):
    posY=[]
    derY=[]
    for splines in curves:
        for s in splines:
            for y in ygrid:
                r = [ v.real for v in (s[1]-y).roots() if abs(v.imag)<1e-15 and 0<=v.real<=1]
                if len(r)>0:
                    for t in r:
                        posY.append((s[0](t),y))
                        derY.append((s[1].deriv())(t))
    return posY,derY
