#!/usr/bin/env python3
"""Check four finite-history certificates and the period-48 ADMM counterexample.
Only the Python standard library and the two accompanying JSON inputs are used.
The continuation and realization arguments are proved in the article.
"""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import argparse
import json

def matrix(rows, cols=None):
    return [[F(0) for _ in range(rows if cols is None else cols)] for _ in range(rows)]


def transpose(A):
    return list(map(list, zip(*A)))


def multiply(A, B):
    return [[sum((x*y for x, y in zip(row, col)), F(0)) for col in zip(*B)] for row in A]


def subtract(A, B):
    return [[x-y for x, y in zip(a, b)] for a, b in zip(A, B)]


def vector_add(a, b, scale=F(1)):
    return [x+scale*y for x, y in zip(a, b)]


def sym_add(A, x, y, scale):
    for i in range(len(A)):
        for j in range(len(A)):
            A[i][j] += scale*(x[i]*y[j]+y[i]*x[j])/2


def constant_columns(size):
    C = matrix(size, 2)
    for i in range(size):
        C[i][i % 2] = F(1)
    return C


def all_pivots(A):
    n = len(A)
    assert A == transpose(A)
    L = matrix(n)
    d, dets = [], []
    product = F(1)
    for i in range(n):
        p = A[i][i]-sum((L[i][k]*L[i][k]*d[k] for k in range(i)), F(0))
        assert p > F(1, 10**7)
        d.append(p)
        product *= p
        dets.append(product)
        L[i][i] = F(1)
        for j in range(i+1, n):
            L[j][i] = (A[j][i]-sum((L[j][k]*L[i][k]*d[k] for k in range(i)), F(0)))/p
    assert multiply(multiply(L, [[d[i] if i == j else F(0) for j in range(n)] for i in range(n)]), transpose(L)) == A
    return d, dets


def operators(data, gamma, mu):
    ell = data['memory'] if 'memory' in data else len(data['P_left'])//2-1
    m, n = 2*(ell+1), 2*(ell+2)
    y, w = [], []
    for j in range(ell+2):
        yy, ww = [F(0)]*n, [F(0)]*n
        yy[2*j], ww[2*j+1] = F(1), F(1)
        y.append(yy)
        w.append(ww)
    r = [vector_add(w[j], w[j+1], F(-1)) for j in range(ell+1)]
    x = [vector_add(y[j], r[j]) for j in range(ell+1)]
    a = [vector_add(vector_add(y[j+1], x[j], F(-1)), w[j+1], -gamma) for j in range(ell+1)]
    c = [vector_add(r[j], w[j+1], gamma) for j in range(ell+1)]
    weights = {(row['block'], row['lag']): F(row['value']) for row in data['weights']}
    Q = matrix(n)
    sym_add(Q, a[0], x[0], F(data['omega']))
    sym_add(Q, c[0], y[0], F(data['omega']))
    sym_add(Q, r[0], r[0], F(data['epsilon_scaled']))
    for block, positions, normals in [('f', x, a), ('g', y, c)]:
        for d in range(1, ell+1):
            sym_add(Q, normals[0], vector_add(positions[0], positions[d], F(-1)), weights.get((block, -d), F(0)))
            sym_add(Q, normals[d], vector_add(positions[d], positions[0], F(-1)), weights.get((block, d), F(0)))
    M = matrix(m)
    for block, positions, normals in [('f', x, a), ('g', y, c)]:
        z = {-1: [F(0)]*m, **{j: positions[j][:m] for j in range(ell)}}
        g = {-1: [F(0)]*m, **{j: normals[j][:m] for j in range(ell)}}
        for (b, i, j), value in mu.items():
            if b == block:
                # G_ij = phi_i-phi_j - <g_j,z_i-z_j>.
                sym_add(M, g[j], vector_add(z[i], z[j], F(-1)), -value)
    return Q, M


def reconstruct_P(A, Q):
    """Use D C=0 to recover the last two rows and columns of P."""
    n, m = len(Q), len(Q)-2
    R = multiply(Q, constant_columns(n))
    X = matrix(m, 2)
    for j in range(m//2):
        for a in range(2):
            for b in range(2):
                X[2*j+a][b] = -sum((R[2*t+a][b] for t in range(j+1)), F(0))
    C0 = constant_columns(m-2)
    B = subtract(X[:-2], multiply(A, C0))
    E = subtract(X[-2:], multiply(transpose(B), C0))
    assert E == transpose(E)
    P = [A[i]+B[i] for i in range(m-2)] + [transpose(B)[i]+E[i] for i in range(2)]
    assert multiply(P, constant_columns(m)) == X
    assert X[-2:] == R[-2:]
    return P


def positive_certificates(cases):
    reports = []
    for item in cases:
        ell = item["memory"]; g0,g1 = map(F, item["interval"])
        omega,eps,zeta = map(F,[item["omega"],item["epsilon"],item["zeta"]])
        dl,dm = item["lag_denominator"],item["flow_denominator"]
        assert 0 < g0 < g1 and min(omega,eps,zeta,dl,dm) > 0
        weights = {}
        assert len(item["lag_numerators"]) == ell
        for d,row in enumerate(item["lag_numerators"],1):
            assert len(row)==4
            for count,(b,s) in zip(row,[("f",-1),("f",1),("g",-1),("g",1)]):
                assert type(count) is int and count>=0
                weights[b,s*d]=F(count,dl)
        h={b:[sum((weights[b,-d]-weights[b,d] for d in range(i+1,ell+1)),F(0))
              for i in range(ell)] for b in ["f","g"]}
        mu={}
        for b,name in [("f","F"),("g","J")]:
            labels=[(b,i,j) for i in range(-1,ell) for j in range(ell) if i!=j]
            packed=item["flow_numerators"][name]
            assert len(packed)==len(labels)==ell**2
            mu.update({k:F(v,dm) for k,v in zip(labels,packed)})
            for i in range(ell):
                mu[b,i,-1]=h[b][i]+sum((mu[b,j,i] for j in range(-1,ell) if i!=j),F(0))-sum((mu[b,i,j] for j in range(ell) if i!=j),F(0))
        assert len(mu)==2*ell*(ell+1) and min(mu.values())>=0
        data={"memory":ell,"omega":str(omega),"epsilon_scaled":str(eps),
              "weights":[{"block":b,"lag":k,"value":str(v)} for (b,k),v in weights.items()]}
        endpoints=[]
        assert len(item["endpoint_principal_blocks"]) == 2
        for gamma,block in zip([g0,g1],item["endpoint_principal_blocks"]):
            denominator=block["denominator"];packed=block["upper_triangle"]
            assert denominator>0 and len(packed)==ell*(2*ell+1)
            A=matrix(2*ell);offset=0
            for i in range(2*ell):
                for j in range(i,2*ell):
                    A[i][j]=A[j][i]=F(packed[offset],denominator);offset+=1
            Q,M=operators(data,gamma,mu)
            P=reconstruct_P(A,Q)
            m,n = 2*(ell+1),2*(ell+2)
            D = [[-value for value in row] for row in Q]
            for i in range(m):
                for j in range(m):
                    D[i][j] -= P[i][j]
                    D[i+2][j+2] += P[i][j]
            assert multiply(D,constant_columns(n)) == matrix(n,2)
            B = subtract(P,M)
            B[0][0] -= zeta
            B[1][1] -= zeta
            dp,_ = all_pivots([row[:m] for row in D[:m]])
            bp,_ = all_pivots(B)
            endpoints.append({'gamma':str(gamma),'pivot_count':len(dp)+len(bp),'all_pivots_strictly_above':'1/10000000'})
        reports.append({"interval":list(map(str,[g0,g1])),"memory":ell,"endpoints":endpoints})
    assert len(reports)==4 and reports[0]['interval'][0]=='3/2' and reports[-1]['interval'][1]=='377/200'
    assert all(a['interval'][1]==b['interval'][0] for a,b in zip(reports,reports[1:]))
    return reports


def sparse_add(a, b, scale=F(1)):
    result = dict(a)
    for key, value in b.items():
        result[key] = result.get(key,F(0))+scale*value
    return {k:v for k,v in result.items() if v}


def sparse_dot(a, G, b):
    return sum((u*G[i][j]*v for i,u in a.items() for j,v in b.items()),F(0))


def bareiss_positive(A):
    """Unpivoted integer elimination; its pivots are leading determinants."""
    B = [list(row) for row in A]
    n, previous = len(B), 1
    assert B == transpose(B)
    determinants = []
    for k in range(n):
        pivot = B[k][k]
        assert isinstance(pivot,int) and pivot >= 1
        determinants.append(pivot)
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator = pivot*B[i][j]-B[i][k]*B[k][j]
                quotient, remainder = divmod(numerator,previous)
                assert remainder == 0
                B[i][j] = quotient
        previous = pivot
    return determinants


def negative_certificate(data):
    n = data["N"]; assert n == 48
    p,q,c=list(data["p"]),list(data["q"]),list(data["c"])
    assert len(p)==len(q)==25 and len(c)==48
    assert all(type(v) is int for v in p+q+c)
    count=len(p)+len(q)+len(c)
    p+=p[-2:0:-1];q+=q[-2:0:-1]
    K = [[0]*(2*n) for _ in range(2*n)]
    for i in range(n):
        for j in range(n):
            K[i][j] = p[(j-i)%n]
            K[n+i][n+j] = q[(j-i)%n]
            K[i][n+j] = K[n+j][i] = c[(j-i)%n]
    J0 = [[2*K[i][j]-(i==j) for j in range(2*n)] for i in range(2*n)]
    determinants = bareiss_positive(J0)
    g0,sigma,gamma=map(F,[data["gamma0"],data["sigma"],data["gamma_hat"]])
    coefficient=data["rank_one_coefficient"];assert coefficient>=0
    assert F(377,200)<gamma<2
    alt = [(-1)**i for i in range(n)]
    h = [-sigma*x for x in alt]+[g0*x for x in alt]
    H = [[F(J0[i][j])+coefficient*h[i]*h[j] for j in range(2*n)] for i in range(2*n)]
    shift = [(i//n)*n+(i+1)%n for i in range(2*n)]
    assert all(H[shift[i]][shift[j]]==H[i][j] for i in range(2*n) for j in range(2*n))
    eta = [{i:F(1)} for i in range(n)]
    w = [{n+i:1/g0} for i in range(n)]
    rho = [sparse_add(w[i],w[i-1],F(-1)) for i in range(n)]
    xi = [sparse_add(eta[i],rho[i]) for i in range(n)]
    a0 = sparse_add(eta[-1],xi[0],F(-1))
    a1 = {i:-v for i,v in w[-1].items()}
    c0,c1 = rho[0],w[-1]
    supports = {'X':[],'Y':[]}
    for block,positions,intercept,slope in [('X',xi,a0,a1),('Y',eta,c0,c1)]:
        for j in [-1]+list(range(1,n)):
            difference = positions[0] if j==-1 else sparse_add(positions[0],positions[j],F(-1))
            b = sparse_dot(intercept,H,difference)
            a = sparse_dot(slope,H,difference)
            supports[block].append({'comparison':j,'slope':a,'intercept':b,'at_endpoint':a*gamma+b})
    assert len(supports['X'])==len(supports['Y'])==48
    assert min(s['slope'] for s in supports['X'])>30000
    assert min(s['slope'] for s in supports['Y'])>20000
    zeros = [s for s in supports['X'] if s['at_endpoint']==0]
    assert len(zeros)==1 and zeros[0]['comparison']==38
    assert min(s['at_endpoint'] for s in supports['X'] if s['comparison']!=38)>F(1,10)
    assert min(s['at_endpoint'] for s in supports['Y'])>2
    active_slope=F(data["active_slope"])
    assert zeros[0]['slope']==active_slope and zeros[0]['intercept']==-active_slope*gamma
    even_difference = sparse_add(xi[0],xi[38],F(-1))
    assert sum((v*h[i] for i,v in even_difference.items()),F(0))==0
    residuals = [sparse_dot(r,H,r) for r in rho]
    assert residuals[0]>0 and len(set(residuals))==1
    for i in range(n):
        assert sparse_add(rho[i],sparse_add(xi[i],eta[i],F(-1)),F(-1))=={}
    return {'independent_integers':count,'positive_integer_leading_determinants':len(determinants),
            'all_leading_determinants_at_least':1,'cyclic_gram_identity_exact':True,
            'gamma_0':str(g0),'sigma':str(sigma),'endpoint':str(gamma),
            'representative_supports':96,'supports_covered_by_exact_cyclic_symmetry':4608,
            'minimum_X_slope':str(min(s['slope'] for s in supports['X'])),
            'minimum_Y_slope':str(min(s['slope'] for s in supports['Y'])),
            'minimum_inactive_X_endpoint':str(min(s['at_endpoint'] for s in supports['X'] if s['comparison']!=38)),
            'minimum_Y_endpoint':str(min(s['at_endpoint'] for s in supports['Y'])),
            'only_X_zero_displacement':-10,'active_support_matches_display':True,
            'active_support_independent_of_rank_one_weight':True,
            'constant_positive_residual_squared':str(residuals[0]),
            'exact_period':48,'period_reason':'positive-definite Gram makes the multiplier vertices distinct'}


if __name__=="__main__":
    if not __debug__:
        raise SystemExit("Run Python without -O: assertions perform exact checks.")
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    pos=(here/"positive.json").read_bytes()
    neg=(here/"cycle48.json").read_bytes()
    cases=positive_certificates(json.loads(pos)["cases"])
    cycle=negative_certificate(json.loads(neg))
    out={"status":"PASS","input_sha256":{"positive.json":sha256(pos).hexdigest(),"cycle48.json":sha256(neg).hexdigest()},
         "positive_cases":cases,"positive_pivots":sum(e["pivot_count"] for c in cases for e in c["endpoints"]),
         "negative_certificate":cycle}
    text=json.dumps(out,indent=2)+"\n"
    if args.output: args.output.write_text(text)
    print(text,end="")
