"""Fixed-point outward-rounded complex roots of unity, by half-angle recursion."""
from math import isqrt
S=10**50
def ceildiv(a,b):return -((-a)//b)
def plus(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def minus(a,b):return plus(a,neg(b))
def times(a,b):
    v=[x*y for x in a for y in b]
    return min(v)//S,ceildiv(max(v),S)
def half(a):return a[0]//2,ceildiv(a[1],2)
def root(a):
    assert a[0]>=0
    l=isqrt(a[0]*S);h=isqrt(a[1]*S)
    return l,h+(h*h<a[1]*S)
def scale(a,k):return (a[0]*k,a[1]*k) if k>=0 else (a[1]*k,a[0]*k)
def roots(n):
    assert n>=4 and n&(n-1)==0
    c=(0,0);s=(S,S);m=4
    while m<n:
        nc=root(half(plus((S,S),c)))
        s=root(half(minus((S,S),c)));c=nc;m*=2
    quadrant=[((S,S),(0,0))]
    for k in range(n//4):
        a,b=quadrant[-1]
        quadrant.append((minus(times(a,c),times(b,s)),plus(times(a,s),times(b,c))))
    assert quadrant[-1][0][0]<=0<=quadrant[-1][0][1]
    assert quadrant[-1][1][0]<=S<=quadrant[-1][1][1]
    table=[]
    for k in range(n):
        q,r=divmod(k,n//4);a,b=quadrant[r]
        table.append([(a,b),(neg(b),a),(neg(a),neg(b)),(b,neg(a))][q])
    return table
