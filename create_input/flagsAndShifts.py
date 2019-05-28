derXshift=2
derYshift=10
boundaryXshift=18
boundaryYshift=22

derXsize=1<<8
derYsize=1<<8

boundaryXsize=1<<4
boundaryYsize=1<<4



def pointType(flag):
    return flag%4

def derXtype(flag):
    return (flag>>derXshift)%derXsize

def derYtype(flag):
    return (flag>>derYshift)%derYsize

def boundaryXtype(flag):
    return (flag>>boundaryXshift)%boundaryXsize

def boundaryYtype(flag):
    return (flag>>boundaryYshift)%boundaryYsize
