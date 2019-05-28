import flagsAndShifts


def setDerFlags(mask):
    setDerXFlags(mask)
    setDerYFlags(mask)


def setDerXFlags(mask):
    imax,jmax = mask.shape
    for i in range(imax):
        for j in range(jmax):
            if(flagsAndShifts.pointType(mask[i,j])==0 or flagsAndShifts.pointType(mask[i,j])==3):
                left=0
                flag=0
                while(left<16 and flag==0):
                    if(0<=j-left):
                        if(flagsAndShifts.pointType(mask[i,j-left])==1):
                            flag=1
                        elif(flagsAndShifts.pointType(mask[i,j-left])==2):
                            left=left+1
                            flag=1
                        else:
                            left=left+1
                    else:
                        flag=1
                left=max(0,left-1)
                flag=0
                right=0
                while(right<16 and flag==0):
                    if(j+right<jmax):
                        if(flagsAndShifts.pointType(mask[i,j+right])==1):
                            flag=1
                        elif(flagsAndShifts.pointType(mask[i,j+right])==2):
                            right=right+1
                            flag=1
                        else:
                            right=right+1
                    else:
                        flag=1
                right=max(0,right-1)
                mask[i,j]=mask[i,j]+((left+16*right)<<flagsAndShifts.derXshift)

def setDerYFlags(mask):
    imax,jmax = mask.shape
    for i in range(imax):
        for j in range(jmax):
            if(flagsAndShifts.pointType(mask[i,j])==0 or flagsAndShifts.pointType(mask[i,j])==3):
                down=0
                flag=0
                while(down<16 and flag==0):
                    if(0<=i-down):
                        if(flagsAndShifts.pointType(mask[i-down,j])==1):
                            flag=1
                        elif(flagsAndShifts.pointType(mask[i-down,j])==2):
                            down=down+1
                            flag=1
                        else:
                            down=down+1
                    else:
                        flag=1
                down=max(0,down-1)
                up=0
                flag=0
                while(up<16 and flag==0):
                    if(i+up<imax):
                        if(flagsAndShifts.pointType(mask[i+up,j])==1):
                            flag=1
                        elif(flagsAndShifts.pointType(mask[i+up,j])==2):
                            up=up+1
                            flag=1
                        else:
                            up=up+1
                    else:
                        flag=1
                up=max(0,up-1)
                mask[i,j]=mask[i,j]+((down+16*up)<<flagsAndShifts.derYshift)
