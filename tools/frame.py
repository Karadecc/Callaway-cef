from PIL import Image; import numpy as np, sys
def load(f): return np.array(Image.open(f).convert('RGB')).astype(int)
def dark(im,th=90): return (im.max(axis=2)<th)
def frame(im):
    d=dark(im); H,W=d.shape
    rows=[(r,d[r].sum()) for r in range(H)]; cols=[(c,d[:,c].sum()) for c in range(W)]
    rows=sorted(rows,key=lambda x:-x[1])[:8]; cols=sorted(cols,key=lambda x:-x[1])[:8]
    return sorted(rows),sorted(cols)
if __name__=='__main__':
    for f in sys.argv[1:]:
        im=load(f); r,c=frame(im); print(f); print(' rows',r); print(' cols',c)
