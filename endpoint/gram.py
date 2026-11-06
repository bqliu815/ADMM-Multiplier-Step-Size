"""Rational Fourier factors ensure PSD by construction; interval supports only.

Each factor (u,v) contributes H=(u,v)(u,v)^*. The plus-sine cross
convention is realized by real coordinates of u exp(i theta phase).
No numerical eigenvalue is used in the acceptance criterion.
"""
import time
from fractions import Fraction as F
import intervals as iv

def add(*terms):
    r={}
    for a,s in terms:
        for k,v in a.items():r[k]=r.get(k,0)+s*v
    return {k:v for k,v in r.items() if v}

def mode_integers(cert):
    n=cert['n']; d=cert['denominator']
    assert type(d) is int and d>0
    modes={}
    for f in cert['factors']:
        j=f['mode']; a,b,c,e=f['coordinates']
        assert type(j) is int and 1<=j<=n//2
        assert all(type(v) is int for v in [a,b,c,e])
        if j==n//2: assert b==e==0
        p,q,re,im=modes.get(j,(0,0,0,0))
        modes[j]=(p+a*a+b*b,q+c*c+e*e,re+a*c+b*e,im+b*c-a*e)
    return modes

def rows(cert,table=None):
    n=cert['n'];u=n//2+1;d=cert['denominator']; modes=mode_integers(cert)
    table=iv.roots(n) if table is None else table
    arr=[(0,0) for _ in range(2*u+n)]
    for k in range(n):
        aa=bb=cc=(0,0)
        for j,(p,q,re,im) in modes.items():
            co,si=table[j*k%n]; mult=1 if j==n//2 else 2
            if k<u:
                aa=iv.plus(aa,iv.scale(co,mult*p))
                bb=iv.plus(bb,iv.scale(co,mult*q))
            cc=iv.plus(cc,iv.plus(iv.scale(co,mult*re),iv.scale(si,mult*im)))
        if k<u: arr[k]=aa;arr[u+k]=bb
        arr[2*u+k]=cc
    return arr,n*d*d*iv.S

def entry(n,row,i,j):
    u=n//2+1
    if i<n and j<n:
        d=(j-i)%n; return row[min(d,n-d)]
    if i>=n and j>=n:
        d=(j-i)%n; return row[u+min(d,n-d)]
    if i<n: return row[2*u+(j-n-i)%n]
    return row[2*u+(i-n-j)%n]

def audit(cert,table=None):
    start=time.monotonic();n=cert['n'];g=F(cert['gamma'])
    row,den=rows(cert,table)
    def ip(a,b):
        ans=(0,0)
        for i,ai in a.items():
            for j,bj in b.items():ans=iv.plus(ans,iv.scale(entry(n,row,i,j),ai*bj))
        return ans
    eta=[{i:1} for i in range(n)];w=[{n+i:1} for i in range(n)]
    r=[add((w[i],1),(w[i-1],-1)) for i in range(n)]
    x=[add((eta[i],1),(r[i],1)) for i in range(n)]
    a0=add((eta[-1],1),(x[0],-1));a1=add((w[-1],-1))
    supports=[]
    for block,h0,h1,current,points in [('X',a0,a1,x[0],x),('Y',r[0],w[-1],eta[0],eta)]:
        for j in ['origin']+list(range(1,n)):
            z={} if j=='origin' else points[j]; delta=add((current,1),(z,-1))
            b0,b1=ip(h0,delta),ip(h1,delta)
            val=iv.plus(iv.scale(b0,g.denominator),iv.scale(b1,g.numerator))
            supports.append({'block':block,'target':j,
                'intercept_bounds':[str(F(k,den)) for k in b0],
                'slope_bounds':[str(F(k,den)) for k in b1],
                'lower':str(F(val[0],den*g.denominator)),
                'upper':str(F(val[1],den*g.denominator))})
    lo=min(supports,key=lambda x:F(x['lower']));hi=min(supports,key=lambda x:F(x['upper']))
    res=ip(r[0],r[0]);positive=all(F(v['lower'])>0 for v in supports) and res[0]>0
    return {'status':'PASS' if positive else 'REJECT','n':n,'gamma':str(g),
        'psd_by_exact_factors':True,'factor_count':len(cert['factors']),
        'reduced_support_count':len(supports),'full_support_count':2*n*n,
        'minimum_support_lower':lo['lower'],'minimum_upper_support':hi['upper'],
        'worst_lower':{'block':lo['block'],'target':lo['target']},
        'negative_upper_count':sum(F(v['upper'])<0 for v in supports),
        'residual_bounds':[str(F(v,den)) for v in res],
        'supports':supports,'cyclicity':'Exact algebraic real Fourier representation',
        'seconds':time.monotonic()-start}

