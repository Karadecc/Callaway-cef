import numpy as np
from frame import load
def color_mask(im, kind):
    r,g,b=im[...,0],im[...,1],im[...,2]
    if kind=='black':   return (r<90)&(g<90)&(b<90)
    if kind=='blue':    return (b>170)&(r<90)&(g<90)
    if kind=='red':     return (r>170)&(g<110)&(b<90)
    if kind=='green':   return (g>80)&(g<160)&(r<110)&(b<90)&(g>r+20)
    if kind=='magenta': return (r>170)&(b>170)&(g<120)
    raise ValueError(kind)
def column_curve(mask, x0,x1,y0,y1, exclude=()):
    m=mask.copy()
    m[:, :x0]=0; m[:, x1:]=0; m[:y0]=0; m[y1:]=0
    for (a,b,c,d) in exclude: m[c:d, a:b]=0
    xs,ys=[],[]
    for c in range(x0,x1):
        rows=np.where(m[:,c])[0]
        if len(rows)>=2: xs.append(c); ys.append(np.median(rows))
    return np.array(xs,float),np.array(ys,float),m
def calib(p1,v1,p2,v2,log=True):
    """pixel->value map through two reference ticks"""
    f=(lambda v: np.log10(v)) if log else (lambda v: v)
    a=(f(v2)-f(v1))/(p2-p1)
    def conv(p):
        u=f(v1)+a*(np.asarray(p)-p1)
        return 10**u if log else u
    return conv
def logbin(x,y,n):
    """average consecutive points in n log-spaced bins of x (reduces column oversampling)"""
    e=np.logspace(np.log10(x.min()),np.log10(x.max())*1.0000001,n+1)
    X,Y=[],[]
    for a,b in zip(e[:-1],e[1:]):
        s=(x>=a)&(x<b)
        if s.sum(): X.append(10**np.mean(np.log10(x[s]))); Y.append(10**np.mean(np.log10(y[s])))
    return np.array(X),np.array(Y)

def clean(xs, ys, win=9, tol=6.0):
    """drop columns whose y deviates by more than tol px from the running median"""
    from scipy.ndimage import median_filter
    med=median_filter(ys, size=win, mode='nearest')
    keep=np.abs(ys-med)<=tol
    return xs[keep], ys[keep]
