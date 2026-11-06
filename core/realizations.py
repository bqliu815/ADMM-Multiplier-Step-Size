#!/usr/bin/env python3
"""Reconstruct strict, smooth, and growing ADMM cycles using rational arithmetic.
This independently implemented reconstruction does not import verify.py.
Requires Python 3.9+; do not use optimization flags.
"""
from fractions import Fraction as R
from pathlib import Path
import argparse, hashlib, json

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path)
args = parser.parse_args()
if not __debug__:
    raise SystemExit("Run Python without -O: assertions perform exact checks.")
SOURCE=HERE/'cycle48.json'
data=json.loads(SOURCE.read_bytes())
N=data['N'];assert N==48
gamma,gamma0,sigma=map(R,[data['gamma_hat'],data['gamma0'],data['sigma']])
epsilon=R(1,100);kappa=R(1000000000001,1000000000000)
p,q,cross=list(data['p']),list(data['q']),list(data['c'])
assert len(p)==len(q)==25 and len(cross)==48
p+=p[-2:0:-1];q+=q[-2:0:-1]
K = [[(p[(j-i)%N] if i<N and j<N else
       q[(j-i)%N] if i>=N and j>=N else
       cross[(j-i)%N] if i<N else cross[(i-j)%N])
      for j in range(2*N)] for i in range(2*N)]
J = [[2*K[i][j] - int(i==j) for j in range(2*N)] for i in range(2*N)]
assert J == list(map(list,zip(*J)))
# Independent fraction-free Schur update, checking all 96 leading minors.
work = [r[:] for r in J]
previous = 1
minors=[]
for k in range(2*N):
    pivot=work[k][k]
    assert pivot > 0
    minors.append(pivot)
    for i in range(k+1,2*N):
        for j in range(k+1,2*N):
            quotient, remainder = divmod(pivot*work[i][j]-work[i][k]*work[k][j],previous)
            assert not remainder
            work[i][j]=quotient
    previous=pivot

h=[-sigma*((-1)**i) for i in range(N)]+[gamma0*((-1)**i) for i in range(N)]

def vec(*terms):
    out={}
    for scale,part in terms:
        for j,a in part.items(): out[j]=out.get(j,R(0))+scale*a
    return {j:a for j,a in out.items() if a}

def e(i): return {i:R(1)}
def dot(a,b,aux_weight=epsilon):
    x={i:v for i,v in a.items() if i<2*N}; y={i:v for i,v in b.items() if i<2*N}
    base=sum((v*w*J[i][j] for i,v in x.items() for j,w in y.items()),R(0))
    base+=data['rank_one_coefficient']*sum((h[i]*v for i,v in x.items()),R(0))*sum((h[i]*v for i,v in y.items()),R(0))
    return base+aux_weight*sum((v*b.get(i,R(0)) for i,v in a.items() if i>=2*N),R(0))

def diff(a,b): return vec((1,a),(-1,b))
def scaled(a,c): return vec((c,a))

def build(aux=True,growth=R(1)):
    etas=[]; ws=[]
    for i in range(N):
        eta=e(i); w=scaled(e(N+i),1/gamma0)
        if aux:
            # (2I-P)^(-1) = sum(P^t/2^(t+1))/(1-2^(-N)).
            for t in range(N):
                coeff=R(1,2**(t+1))/(1-R(1,2**N))
                eta=vec((1,eta),(-2*coeff,e(2*N+(i-t)%N)),((2-gamma)*coeff,e(2*N+(i-t-1)%N)))
            w=vec((1,w),(1,e(2*N+i)))
        etas.append(eta); ws.append(w)
    nus=[scaled(w,gamma) for w in ws]
    residuals=[scaled(vec((1,nus[i]),(-1/growth,nus[i-1])),1/gamma) for i in range(N)]
    xis=[vec((1,etas[i]),(1,residuals[i])) for i in range(N)]
    aa=[vec((1/growth,etas[i-1]),(-1/growth,nus[i-1]),(-1,xis[i])) for i in range(N)]
    cc=[vec((1/growth,nus[i-1]),(1,residuals[i])) for i in range(N)]
    return etas,nus,residuals,xis,aa,cc

def support_summary(z,normals):
    values=[dot(normals[0],z[0])]+[dot(normals[0],diff(z[0],z[j])) for j in range(1,N)]
    # Orthogonal cyclic shifts reduce N^2 comparisons to these 48 representatives.
    return values

# Verify coefficient and metric cyclic invariance; this is exact for every phase.
shift=lambda i: (i//N)*N+(i+1)%N
assert all(J[shift(i)][shift(j)]==J[i][j] for i in range(2*N) for j in range(2*N))
assert all(h[shift(i)]==-h[i] for i in range(2*N))

def check_cycle(data,growth):
    eta,nu,r,x,a,c=data
    for i in range(N):
        assert vec((1,x[i]),(1,a[i]),(-1/growth,eta[i-1]),(1/growth,nu[i-1]))=={}
        assert vec((1,eta[i]),(1,c[i]),(-1,x[i]),(-1/growth,nu[i-1]))=={}
        assert vec((1,nu[i]),(-1/growth,nu[i-1]),(-gamma,r[i]))=={}
        for family in data:
            assert {shift(j):v for j,v in family[i].items()}==family[(i+1)%N]
    assert dot(r[0],r[0])>0

base=build(False)
check_cycle(base,R(1))
baseX=support_summary(base[3],base[4]); baseY=support_summary(base[0],base[5])
assert min(baseX)==0 and baseX[38]==0 and sum(v==0 for v in baseX)==1
assert min(baseY)>2 and min(v for v in baseX if v)>R(1,10)
assert all(dot(diff(base[1][i],base[1][j]),diff(base[1][i],base[1][j]))>0 for i in range(N) for j in range(i))

basewprev=scaled(base[1][-1],1/gamma)
slopesX=[dot(scaled(basewprev,-1),base[3][0])]+[dot(scaled(basewprev,-1),diff(base[3][0],base[3][j])) for j in range(1,N)]
slopesY=[dot(basewprev,base[0][0])]+[dot(basewprev,diff(base[0][0],base[0][j])) for j in range(1,N)]
assert min(slopesX)>30000 and min(slopesY)>20000
assert slopesX[38]==R(287451196000,1186923)

strict=build(True)
check_cycle(strict,R(1))
sX=support_summary(strict[3],strict[4]); sY=support_summary(strict[0],strict[5])
assert min(sX)>0 and min(sY)>0
# Check the auxiliary first-block identity a_i^circle=xi_i^circle exactly.
assert all({j:v for j,v in strict[3][i].items() if j>=2*N}=={j:v for j,v in strict[4][i].items() if j>=2*N} for i in range(N))

# Global smooth interpolation is checked through its finite sufficient inequalities.
def smooth(z,normals):
    alpha=dot(normals[0],z[0]); theta=alpha/2
    supports=support_summary(z,normals)
    delta=min([theta]+supports[1:])
    D=max([dot(z[0],z[0])]+[dot(diff(z[0],z[j]),diff(z[0],z[j])) for j in range(1,N)])
    mu=delta/(2*D)
    pp=[vec((1,normals[i]),(-mu,z[i])) for i in range(N)]
    maxp=max([dot(pp[0],pp[0])]+[dot(diff(pp[0],pp[j]),diff(pp[0],pp[j])) for j in range(1,N)])
    K0=max(R(1),maxp/delta)
    # Two origin orientations and all phase displacements, plus self equalities.
    comparisons=[(theta,z[0],pp[0]),(alpha-theta,z[0],pp[0])]
    comparisons += [(supports[j],diff(z[0],z[j]),diff(pp[0],pp[j])) for j in range(1,N)]
    margins=[]
    for gap,dz,dp in comparisons:
        lhs=gap-mu*dot(dz,dz)/2
        rhs=dot(dp,dp)/(2*K0)
        assert lhs>=3*delta/4 and lhs>=rhs
        margins.append(lhs-rhs)
    return {'mu':str(mu),'K0':str(K0),'L':str(K0+mu),'theta':str(theta),'minimum_interpolation_margin':str(min(margins))}

smoothX=smooth(strict[3],strict[4]); smoothY=smooth(strict[0],strict[5])
grow=build(True,kappa)
check_cycle(grow,kappa)
gX=support_summary(grow[3],grow[4]); gY=support_summary(grow[0],grow[5])
assert min(gX)>0 and min(gY)>0
# Active normal comparisons at a vertex, essential for one fixed homogeneous objective.
for z,a in [(grow[3],grow[4]),(grow[0],grow[5])]:
    alpha=dot(a[0],z[0]); assert alpha>0
    assert all(dot(a[j],z[0])<alpha for j in range(1,N))

# Universal identity: degree at most two in gamma and homogeneous quadratic in four states.
# Polarization at every state basis vector and basis sum recovers every state coefficient.
# Three distinct rational gamma values then determine each degree-two coefficient polynomial.
def energy_error(v,g):
    y,w,yold,wold=v
    r=w-wold; x=y+r; a=yold-g*wold-x; c=g*wold+r; dy=y-yold
    lhs=g*g*(w*w-wold*wold)+g*(y*y-yold*yold)
    rhs=-2*g*(a*x+c*y)-g*dy*dy-2*g*dy*r-g*(2-g)*r*r
    return lhs-rhs
basis=[[R(int(i==j)) for i in range(4)] for j in range(4)]
for gamma_probe in [R(0),R(1),R(2)]:
    for i in range(4):
        assert energy_error(basis[i],gamma_probe)==0
        for j in range(i):
            assert energy_error([x+y for x,y in zip(basis[i],basis[j])],gamma_probe)==0
result={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
 'author_verifier_imported':False,'gamma_search':False,
 'independent_positive_definite_integer_minors':len(minors),
 'base_active_displacement':-10,'distinct_multiplier_states':48,
 'minimum_base_X_slope':str(min(slopesX)),'minimum_base_Y_slope':str(min(slopesY)),
 'fixed_strictification_weight':str(epsilon),'strict_X_margin':str(min(sX)),'strict_Y_margin':str(min(sY)),
 'smooth_X':smoothX,'smooth_Y':smoothY,
 'fixed_growth_factor':str(kappa),'growth_X_margin':str(min(gX)),'growth_Y_margin':str(min(gY)),
 'growth_residual_squared':str(dot(grow[2][0],grow[2][0])),
 'growth_objective_coefficient':str((dot(grow[4][0],grow[3][0])+dot(grow[5][0],grow[0][0]))/2),
 'all_coefficient_recurrences_and_phase_symmetry_exact':True,
 'universal_energy_identity_symbolically_zero':True,
 'scope':'finite arithmetic and polynomial identities; global conjugacy and infinite replay are analytic arguments in the article',
 'status':'PASS'}
if args.output: args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'], 'source_sha256':result['source_sha256'], 'strictification_weight':str(epsilon), 'growth_factor':str(kappa), 'strict_supports':True, 'smooth_interpolation':True, 'growth_recursions':True},indent=2))
