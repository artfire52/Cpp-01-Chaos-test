#!/usr/bin/env python3
#./gaussian_map.py chaotic 3.7 5000 0.5
#./gaussian_map.py regular 2.8 5000 0.5
import sys
def logistic(r,n,x0): #un+1=r(un)*(1-un)
    res=[x0]
    for i in range(n):
       res.append(res[-1]*(1-res[-1])*r) 
    return res

def ouputfile(filename,r,n,x0):
    f=open(filename,"w")
    for i,l in enumerate(logistic(r,n,x0)):
        f.write(str(i)+" "+str(l)+"\n") 
    f.close()

if __name__== "__main__":
    argv=sys.argv
    if(len(argv))!=5:
        exit()
    filename=argv[1]
    r=eval(argv[2])
    n=int(eval(argv[3]))
    x0=float(eval(argv[4]))
    ouputfile(filename,r,n,x0)
