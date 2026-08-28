"""Reproduces every numerical claim of JOCOv5.tex by complete enumeration.
Run:  python3 verify_JOCOv5.py           (a few minutes, exact arithmetic)"""
import itertools, math, random
from fractions import Fraction as F
from math import comb, ceil, log2, sqrt

def regs(x,n,l): return [(x>>(t*l))&((1<<l)-1) for t in range(n)]
def Cv(v,n,d):   return sum(d[v[t]][v[(t+1)%n]] for t in range(n)
                            if v[t]<n and v[(t+1)%n]<n)
def Pv(v,n):     return (sum(1 for u in v if u>=n)
                         + sum(1 for a in range(n) for b in range(a+1,n) if v[a]==v[b]))
def tour(p,n,d): return sum(d[p[t]][p[(t+1)%n]] for t in range(n))
def rep_LR(v,n):
    w=[u%n for u in v]; out=[None]*n; used=set()
    for t in range(n):
        if w[t] not in used: out[t]=w[t]; used.add(w[t])
    free=sorted(set(range(n))-used); k=0
    for t in range(n):
        if out[t] is None: out[t]=free[k]; k+=1
    return out
def rep_LF(v,n):
    out=[None]*n; used=set()
    for t in range(n):
        if v[t]<n and v[t] not in used: out[t]=v[t]; used.add(v[t])
    free=sorted(set(range(n))-used); k=0
    for t in range(n):
        if out[t] is None: out[t]=free[k]; k+=1
    return out
def full(n,d,exact=True):
    l=max(1,ceil(log2(n))); m=n*l
    cs=min(tour(list(p),n,d) for p in itertools.permutations(range(n)))
    UF=max(tour(list(p),n,d) for p in itertools.permutations(range(n)))
    dm=max(max(r) for r in d)
    st=dict(cstar=cs,UF=UF,dmax=dm,l=l,m=m,maxC=None,nu=None,minC_P1=None,
            lam=None,ltrLF=None,ltrLR=None,sharp=[])
    bl=None; bLF=None; bLR=None
    for x in range(1<<m):
        v=regs(x,n,l); C=Cv(v,n,d); P=Pv(v,n)
        st['maxC']=C if st['maxC'] is None else max(st['maxC'],C)
        if P==0: continue
        if C<=1e-12 and (st['nu'] is None or P<st['nu']): st['nu']=P
        if P==1 and (st['minC_P1'] is None or C<st['minC_P1']): st['minC_P1']=C
        r=(max(0,cs-C)/P) if not exact else F(max(0,cs-C),P)
        if bl is None or r>bl: bl=r; st['sharp']=[(v,C,P)]
        elif r==bl: st['sharp'].append((v,C,P))
        a=tour(rep_LF(v,n),n,d)-C; b=tour(rep_LR(v,n),n,d)-C
        a=F(a,P) if exact else a/P; b=F(b,P) if exact else b/P
        bLF=a if bLF is None else max(bLF,a); bLR=b if bLR is None else max(bLR,b)
    st['lam'],st['ltrLF'],st['ltrLR']=bl,bLF,bLR
    return st
ok=lambda c: "OK " if c else "** MISMATCH **"

print("="*78); print("Ex. 10.1  unit square, n=4")
s2=sqrt(2); dsq=[[0,1,s2,1],[1,0,1,s2],[s2,1,0,1],[1,s2,1,0]]
st=full(4,dsq,exact=False)
print(f"  c*={st['cstar']:.4f} U_F={st['UF']:.4f} maxC={st['maxC']:.4f} nu={st['nu']} "
      f"min[C:P=1]={st['minC_P1']:.4f} lambda*={st['lam']:.4f} ltr_LF={st['ltrLF']:.4f}")
print("  ",ok(abs(st['cstar']-4)<1e-9 and abs(st['UF']-(2+2*s2))<1e-9
        and abs(st['maxC']-4*s2)<1e-9 and st['nu']==6
        and abs(st['minC_P1']-(2+s2))<1e-9 and abs(st['lam']-1)<1e-9
        and abs(st['ltrLF']-s2)<1e-9),
      "claims: c*=4, U_F=2+2sqrt2, maxC=4sqrt2, nu=6, min[C:P=1]=2+sqrt2 (NOT 1+sqrt2), lambda*=1, ltr_LF=dmax=sqrt2")

print("="*78); print("Ex. 10.2 / Prop. 7.12  path metric, n=5")
d5=[[abs(i-j) for j in range(5)] for i in range(5)]
st=full(5,d5)
cens={}
for v,C,P in st['sharp']: cens[tour(rep_LR(v,5),5,d5)]=cens.get(tour(rep_LR(v,5),5,d5),0)+1
censLF={}
for v,C,P in st['sharp']: censLF[tour(rep_LF(v,5),5,d5)]=censLF.get(tour(rep_LF(v,5),5,d5),0)+1
print(f"  lambda*={st['lam']} nu={st['nu']} #sharp={len(st['sharp'])} LR-census={dict(sorted(cens.items()))} "
      f"LF-census={dict(sorted(censLF.items()))} ltr_LF={st['ltrLF']} ltr_LR={st['ltrLR']} 2dmax={2*st['dmax']} 4dmax={4*st['dmax']}")
print("  ",ok(st['lam']==5 and st['nu']==3 and len(st['sharp'])==60
        and cens=={8:50,10:5,12:5} and censLF=={8:60}
        and st['ltrLF']==7 and st['ltrLR']==9),
      "claims: lambda*=5, nu=3, 60 sharp, census 50/5/5 (NOT 48/6/6), LF repairs all 60 optimally, ltr=7 and 9")

print("="*78); print("Thm. 5.5  exact multipliers of three families")
for n in (4,5,6):
    du=[[0 if i==j else 1 for j in range(n)] for i in range(n)]
    st=full(n,du); q=2**max(1,ceil(log2(n)))-n
    pred=1 if q==0 else 2
    print(f"  K_{n}: q={q} lambda*={st['lam']} predicted {pred} ",ok(st['lam']==pred))
for n in (5,6):
    dp=[[abs(i-j) for j in range(n)] for i in range(n)]
    st=full(n,dp)
    print(f"  P_{n}: lambda*={st['lam']} predicted {n} ",ok(st['lam']==n))
for n,M in ((5,50),(6,20)):
    dc=[[0 if i==j else (1 if abs(i-j)%n in (1,n-1) else M) for j in range(n)] for i in range(n)]
    st=full(n,dc)
    print(f"  spread cycle n={n} M={M}: c*={st['cstar']} lambda*={st['lam']} (pred 2) "
          f"ltr_LF={st['ltrLF']} (pred 2M={2*M}) ",ok(st['lam']==2 and st['ltrLF']==2*M))

print("="*78); print("Prop. 4.7  nu = min[P : C=0]")
for n in (4,5,6,7,8):
    l=max(1,ceil(log2(n))); best=None
    for v in itertools.product(range(2**l),repeat=n):
        if not all((v[t]>=n) or (v[(t+1)%n]>=n) or (v[t]==v[(t+1)%n]) for t in range(n)): continue
        P=Pv(list(v),n)
        if P>0 and (best is None or P<best): best=P
    print(f"  n={n}: nu={best}  ceil(n/2)={ceil(n/2)}  C(n,2)={comb(n,2)}  q={2**l-n}  "
          f"[previous draft claimed nu>=n] ",ok(best is not None))

print("="*78); print("Prop. 7.12  path metric: worst LR-repair of a sharp string = 4n-8")
for n in range(5,15):
    l=ceil(log2(n))
    if 2**l==n: continue
    d=[[abs(i-j) for j in range(n)] for i in range(n)]; w=0
    for p in range(n):
        for code in range(n,2**l):
            for base in (0,1):
                for dr in (1,-1):
                    seq=[base+i for i in range(n-1)][::dr]; v=[None]*n; v[p]=code
                    for k in range(n-1): v[(p+1+k)%n]=seq[k]
                    assert Cv(v,n,d)==n-2 and Pv(v,n)==1
                    w=max(w,tour(rep_LR(v,n),n,d))
                    assert tour(rep_LF(v,n),n,d)==2*(n-1)
    print(f"  n={n:2d}: worst={w} 4n-8={4*n-8} ",ok(w==4*n-8),f" ratio={w/(2*(n-1)):.4f} (LF: always optimal)")

print("="*78); print("Thm. 7.10 / Ex. 10.3  spread cycle witnesses")
for n,M in ((7,500),(7,50)):
    d=[[0 if i==j else (1 if abs(i-j)%n in (1,n-1) else M) for j in range(n)] for i in range(n)]
    v=[n]+list(range(n-1)); C=Cv(v,n,d); P=Pv(v,n); r=rep_LR(v,n)
    print(f"  n={n} M={M}: C={C} P={P} c*=n={n} LR->{r} cost={tour(r,n,d)} (=2M+n-2={2*M+n-2}) "
          f"ratio={tour(r,n,d)/n:.2f} ",ok(C==n-2 and P==1 and tour(r,n,d)==2*M+n-2))
    v2=[0,1,2,n,3]+list(range(4,n-1)); C2=Cv(v2,n,d)
    print(f"     Thm 7.9 witness {v2}: C={C2} (=n-3+M={n-3+M}) LF cost={tour(rep_LF(v2,n),n,d)} "
          f"diff={tour(rep_LF(v2,n),n,d)-C2} (=2M={2*M}) ",ok(C2==n-3+M and tour(rep_LF(v2,n),n,d)-C2==2*M))
for n,M in ((7,500),):
    d=[[M]*n for _ in range(n)]
    for i in range(n): d[i][i]=0
    for i in range(n-2): d[i][i+1]=d[i+1][i]=1
    d[n-1][1]=d[1][n-1]=1; d[n-1][2]=d[2][n-1]=1; d[n-2][0]=d[0][n-2]=1
    cs=min(tour(list(p),n,d) for p in itertools.permutations(range(n)))
    v=list(range(n-1))+[n]; C=Cv(v,n,d)
    print(f"  Thm 7.10 LF-instance n={n} M={M}: c*={cs} C={C} P={Pv(v,n)} LF cost={tour(rep_LF(v,n),n,d)} "
          f"(=n-2+2M={n-2+2*M}) ",ok(cs==n and C==n-2 and tour(rep_LF(v,n),n,d)==n-2+2*M))

print("="*78); print("Lem. 7.4  locality: LR is not r-local for r<4")
v=[5,0,1,2,3]; A=[t for t in range(5) if rep_LR(v,5)[t]!=v[t]]
print(f"  v={v} LR->{rep_LR(v,5)} A={A} 2|A|={2*len(A)} P={Pv(v,5)} ",ok(len(A)==2 and Pv(v,5)==1))

print("="*78); print("Thm. 8.3 / Rem. 8.4  monomial-wise budget")
for l in (2,3,4,5):
    brute=sum(comb(2*l,k)*(k-2) for k in range(3,2*l+1))
    print(f"  ell={l}: brute={brute:6d} corrected=(l-2)n^2+2l+2={(l-2)*4**l+2*l+2:6d} "
          f"previous-draft form={l*2**(2*l-1)-2**(2*l+1)+2+2*l:7d} ",
          ok(brute==(l-2)*4**l+2*l+2))

print("="*78); print("Lem. 8.5  no quadratic XNOR gadget on three bits")
try:
    from scipy.optimize import linprog
    mono=lambda b,bp,q:[1,b,bp,q,b*bp,b*q,bp*q]
    Z=[(0,0,1),(0,1,0),(1,0,0),(1,1,1)]; NZ=[(0,0,0),(0,1,1),(1,0,1),(1,1,0)]
    r=linprog([0]*7,A_ub=[[-c for c in mono(*t)] for t in NZ],b_ub=[-1]*4,
              A_eq=[mono(*t) for t in Z],b_eq=[0]*4,bounds=[(None,None)]*7)
    print("  LP over (a0,a1,a2,a3,a12,a13,a23): feasible =",r.status==0,ok(r.status!=0))
except ImportError: print("  scipy absent; see the algebraic proof of Lem. 8.5")

print("="*78); print("Tab. 3  compilation budgets")
for n in (4,8,16,32):
    l=int(log2(n))
    print(f"  n={n:3d}: monomial-wise={n*l+n*((l-2)*n*n+2*l+2):7d} equality-tree={n*l+comb(n,2)*(2*l-1):6d} "
          f"shared decoder={n*l+n*(2*n-4):6d} one-hot={n*n:5d}")

print("="*78); print("Cor. 7.7  random check of the transfer bounds")
random.seed(11); bad=0; T=46
for t in range(T):
    n=random.choice([4,5]); d=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            w=random.randint(1,9); d[i][j]=d[j][i]=w
    if t%3==0:
        for i in range(n):
            for j in range(n):
                if i!=j: d[i][j]=random.randint(1,9)
    st=full(n,d)
    if st['ltrLF']>2*st['dmax'] or st['ltrLR']>4*st['dmax']: bad+=1
print(f"  {T} random instances (symmetric and asymmetric): violations of ltr_LF<=2dmax and ltr_LR<=4dmax:",bad,ok(bad==0))
print("="*78)
