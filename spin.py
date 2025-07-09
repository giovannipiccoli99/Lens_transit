import numpy as np
import matplotlib.pyplot as pl

def rotX(f):
    cs,sn = np.cos(f), np.sin(f)
    matr = np.eye(3)
    matr[1,1] = matr[2,2] = cs
    matr[1,2] = -sn
    matr[2,1] = sn
    return matr

def rotY(f):
    cs,sn = np.cos(f), np.sin(f)
    matr = np.eye(3)
    matr[0,0] = matr[2,2] = cs
    matr[2,0] = -sn
    matr[0,2] = sn
    return matr


def rotZ(f):
    cs,sn = np.cos(f), np.sin(f)
    matr = np.eye(3)
    matr[0,0] = matr[1,1] = cs
    matr[0,1] = -sn
    matr[1,0] = sn
    return matr

def merid(lat):
    return np.array([0,np.cos(lat),np.sin(lat)]).T

wye = np.array([0,1,0]).T

lon = np.linspace(0,2*np.pi,50)

def frame(ay,az):
    pl.clf()
    for lat in np.linspace(-np.pi/2,np.pi/2,12):
        p = np.array([rotZ(az) @ rotY(ay) @ rotZ(l) @ rotX(lat) @ wye for l in lon])
        if lat > 0:
            col = 'green'
        else:
            col = 'cyan'
        pl.plot(p[:,1],p[:,2],color=col)
    pl.gca().set_aspect(1)
    pl.pause(.1)
    
az = 0.5
for ay in np.linspace(0,2*np.pi,256):
    frame(ay,az)
    
ay = 0.5
for az in np.linspace(0,2*np.pi,256):
    frame(ay,az)
