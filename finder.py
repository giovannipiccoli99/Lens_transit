import numpy as np

def find(lensmap,src,zc,wd,toplev=True):
    eps = 1e-4
    if toplev:
        N = 1000
    else:
        N = 10
    x = np.linspace(zc.real,zc.real+wd,N+1)
    y = np.linspace(zc.imag,zc.imag+wd,N+1)
    x,y = np.meshgrid(x,y)
    z = x + 1j*y
    wd /= N
    s = lensmap(z) - src
    a = s[:-1,:-1].real * s[1:,1:].real < 0
    b = s[1:,:-1].real * s[:-1,1:].real < 0
    c = s[:-1,:-1].imag * s[1:,1:].imag < 0
    d = s[1:,:-1].imag * s[:-1,1:].imag < 0
    f = z[:-1,:-1][(a+b)*(c+d)]
    #print('candidates at wd = ',wd)
    #print(f)
    lyst = []
    if wd < eps:
        for p in f:
            lyst.append(p)
    else:
        for p in f:
            lyst += find(lensmap,src,p,wd,toplev=False)
    if not toplev:
        return lyst
    tlyst = []
    for p in lyst:
        for q in tlyst:
            if abs(p-q) < 2*eps:
                break
        else:
            tlyst.append(p)
    return tlyst

if __name__ == '__main__':

    import matplotlib.pyplot as pl
    
    def lensmap(z):
        return z - z/abs(z+1e-20) + z.real/2
        
    zc = -2 - 2j
    wd = 4
    
    for f in np.linspace(-1,1,100):
        src = f * (1+1j)
        lyst = find(lensmap,src,zc,wd,True)
        pl.clf()
        pl.plot(src.real,src.imag,'o',color='cyan')
        for im in lyst:
            pl.plot(im.real,im.imag,'o',color='red')
        pl.xlim(zc.real,(zc+wd).real)
        pl.ylim(zc.imag,zc.imag+wd)
        pl.gca().set_aspect(1)
        pl.pause(0.1)
    
