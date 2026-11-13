from frame import *
def ticks_bottom(im, bottom, x0, x1, inside=True):
    d=dark(im); out=[]
    for c in range(x0+3,x1-2):
        n=0; r=bottom-1 if inside else bottom+1; step=-1 if inside else 1
        while d[r,c]: n+=1; r+=step
        if n>=4: out.append((c,n))
    return out
def ticks_left(im, left, y0, y1, inside=True):
    d=dark(im); out=[]
    for r in range(y0+3,y1-2):
        n=0; c=left+1 if inside else left-1; step=1 if inside else -1
        while d[r,c]: n+=1; c+=step
        if n>=4: out.append((r,n))
    return out
def group(t):
    g=[];cur=[t[0]]
    for a in t[1:]:
        if a[0]-cur[-1][0]<=1: cur.append(a)
        else: g.append(cur); cur=[a]
    g.append(cur)
    return [(round(np.mean([x[0] for x in c]),1), max(x[1] for x in c)) for c in g]
